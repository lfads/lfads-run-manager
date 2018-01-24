# Goals
# -----
#
# What we want is to launch and monitor multiple shell scripts simultaneously.
# We have a list of tasks. Each task has a certain amount of memory that it needs in GPU memory.
#
# Launch tensorboard with all runs queued.
#
# At the start, we loop through the queue and determine if the next task can run on any GPU. If so, we launch it on that GPU. On each tick (1 second), we perform this basic operation where we launch the first item in the queue that can be run on a GPU. Print a message when a new process is launched on which GPU.
#
# While this is happening, monitor each of the output processes that are in play and dump output with appropriate identifying queues.

import subprocess
import os
import csv
import re
import numpy as np
import time
import shlex
from multiprocessing import Process, Queue, cpu_count, Lock
from Queue import Empty
import sys, traceback

# needed for synchronization to ensure tmux isn't called by multiple processes simultaneously
mutex = Lock()

class GpuStatus(object):
    def __init__(self, index, name, memfree, memtotal, uuid, num_tasks=0):
        self.index = index
        self.name = name
        self.memfree = memfree
        self.memtotal = memtotal
        self.uuid = uuid
        self.num_tasks = num_tasks

    def __repr__(self):
        return 'GpuStatus {}:"{}" free {}/{} MiB, {} tasks'.format(self.index, self.name, self.memfree, self.memtotal, self.num_tasks)

    def incr_num_tasks(self):
        self.num_tasks += 1

    def decr_num_tasks(self):
        if self.num_tasks > 0:
            self.num_tasks -= 1

class Task(object):
    index = None
    name = ''
    command = '' # shell command
    outfile = None # output from stdout and stderr
    donefile = None # created when task completed
    memory_req = 0 # memory needed on GPU
    running_on_gpu = None # which GPU index running on
    has_finished = False
    has_failed = False
    skipped_donefile_exists = False

    process = None
    popen = None

    def is_running(self):
        return self.running_on_gpu is not None and not self.has_finished

    def __init__(self, index, name='', command='', memory_req=0,
                 outfile=None, donefile=None, tmux_session=None):
        self.index = index
        self.name = name
        self.command = command
        self.memory_req = memory_req
        self.outfile = outfile
        self.donefile = donefile
        if tmux_session is None:
            self.tmux_session = tmux_session

    def __repr__(self):
        if self.has_finished:
            if self.skipped_donefile_exists:
                status = 'skipped, donefile exists'
            else:
                status = 'finished'
        elif self.is_running():
            status = 'running on GPU {}'.format(self.running_on_gpu)
        else:
            status = 'not running'
        return 'Task {} {}: mem req {}, {}'.format(self.index, self.name, self.memory_req, status)

    def mark_finished_if_donefile_exists(self):
        if self.donefile:
            if os.path.exists(self.donefile):
                self.has_finished = True
                self.skipped_donefile_exists = True

    def delete_donefile(self):
        if self.donefile and os.path.exists(donefile):
            os.remove(self.donefile)

def query_gpu_status():
    """
    Calls nvidia_smi to poll GPU free memory

    Returns:
        gpu_status (list of dicts) : fields memfree (MiB), memtotal (MiB), name
    """
    nvidia_out = subprocess.check_output(['nvidia-smi',
                                          '--query-gpu=name,memory.free,memory.total,gpu_uuid',
                                          '--format=csv,noheader,nounits'])
    nvidia_with_header = "name,memfree,memtotal,uuid\n" + nvidia_out


    # parse output of nvidia-smi query
    reader = csv.DictReader(nvidia_with_header.splitlines())
    gpu_status = []
    for idx, info in enumerate(reader):
        gpu_status.append(GpuStatus(index=idx, name=info['name'].strip(),
                                    memfree=float(info['memfree']),
                                    memtotal=float(info['memtotal']),
                                    uuid=info['uuid'].strip()))

    # now figure out how many python processes are running on each gpu
    nvidia_out = subprocess.check_output(['nvidia-smi',
                                          '--query-compute-apps=process_name,gpu_uuid',
                                          '--format=csv,noheader,nounits'])
    nvidia_with_header = "pname,uuid\n" + nvidia_out

    # parse output of nvidia-smi query
    reader = csv.DictReader(nvidia_with_header.splitlines())

    # increment num tasks on python tasks by matching up GPU uuids
    for task in reader:
        if task['pname'] == 'python':
            task_gpu_uuid = task['uuid'].strip()
            [s.incr_num_tasks() for s in gpu_status if s.uuid == task_gpu_uuid]

    return gpu_status

def find_gpu_ready_for_task(gpu_status, task):
    """
    Returns the index of thegpu in gpu_status with memfree >= memory_req
    with the fewest running tasks

    Args:
        gpu_status (list of GpuStatus)
        task (list of Task)

    Returns:
        gpu.index (int)

    """

    # of the set of GPUs with sufficient memory
    gpu_eligible = [gpu for gpu in gpu_status if gpu.memfree >= task.memory_req]
    if not gpu_eligible:
        return None

    # and the fewest running tasks. This also ensures that every GPU will be
    # used at least once before tasks are doubled up
    gpu = sorted(gpu_eligible, key = lambda gpu: gpu.num_tasks)[0]
    return gpu.index

def find_gpu_by_index(gpu_list, index):
    return [gpu for gpu in gpu_list if gpu.index == index][0]

def pick_next_task(gpu_status, tasks):
    """
    Returns the index of the next task that isn't yet running and for which
    there is a GPU that can run it

    Args:
        gpu_status : info about gpu memory usage
        tasks (list of Tasks) : list of tasks to consider

    Returns:
        task_index
        gpu_index

    """
    for task_index, task in enumerate(tasks):
        if not task.has_finished and not task.is_running():
            gpu_index = find_gpu_ready_for_task(gpu_status, task)
            if gpu_index is not None:
                return (task_index, gpu_index)

    return (None, None)

def check_all_tasks_completed_or_running(tasks):
    return all([task.has_finished or task.is_running() for task in tasks])

def check_num_tasks_running(tasks):
    return sum([1 if task.is_running() else 0 for task in tasks])

def check_all_tasks_complete(tasks):
    return all([task.has_finished for task in tasks])

def build_task_list(task_specs):
    return [Task(index=idx, **spec) for (idx, spec) in enumerate(task_specs)]

def generate_tmux_command(name, command):
    """
    Generates a tmux command that will run command in tmux session name
    and wait for it to complete before returning. The bash trap ensures that the
    popen subprocess will complete even if the command in tmux is Ctrl-C'ed
    """
#     return "tmux new-session -ds {name} 'trap \"tmux wait-for -S {name}_done\" SIGHUP SIGTERM SIGINT; {command}; tmux wait-for -S {name}_done' && tmux wait-for {name}_done".format(name=name, command=command)
    return "tmux new-session -ds {name} 'export PATH={path}; {command}'".format(path=os.environ['PATH'], name=name, command=command)

def generate_tee_command(command, outfile, append=True):
    """Appends tee redirection for stdout and stderr to a command.
    Command may be a list of strings, in which case each subcommand will be tee'd
    and the set will be joined by semicolons"""
    if append:
        tee = 'tee -a'
    else:
        tee = 'tee'

    if type(command) is str:
        return '{} 2>&1 | {} {}'.format(command, tee, outfile)
    else:
        commands = ['{} 2>&1 | {} {}'.format(cmd, tee, outfile) for cmd in command]
        return '; '.join(commands)

def follow_file(file):
    """Generator that reads lines appended to the end of a file"""

    file.seek(0,2) # begin yielding at end of file in case we're using tee -a
    while True:
        line = file.readline()
        if not line:
            yield '' # need to yield '' so we don't close the generator
        else:
            yield line

def get_tail_file(filepath, nlines=3):
    return subprocess.check_output(shlex.split('tail -n {} {}'.format(nlines, filepath)))

def touch_file(outfile):
    # mkdir -p
    if not os.path.exists(os.path.dirname(outfile)):
        try:
            os.makedirs(os.path.dirname(outfile))
        except OSError as exc: # Guard against drace condition
            if exc.errno != errno.EEXIST:
                raise

    # touch outfile to ensure it exists
    try:
        temp = open(outfile, 'a')
    except IOError as exc:
        print('Error opening out file: {}'.format(exc))

    temp.close()

def get_list_tmux_sessions():
    with mutex:
        return subprocess.check_output(shlex.split("tmux list-sessions -F '#{session_name}'")).splitlines()

def check_tmux_session_exists(session):
    return session in get_list_tmux_sessions()

def get_list_tmux_sessions_name_starts_with(prefix):
    return filter(lambda sess: sess.startswith(prefix), get_list_tmux_sessions())

class TaskCompletedMessage(object):
    def __init__(self, task_index, success, tail):
        self.task_index = task_index
        self.success = success
        self.tail = tail

class TaskStartedMessage(object):
    def __init__(self, task_index, pid, tmux_session):
        self.task_index = task_index
        self.pid = pid
        self.tmux_session = tmux_session

class TaskExceptionMessage(object):
    def __init__(self, task_index):
        self.task_index = task_index
        T, V, TB = sys.exc_info()
        self.message = ''.join(traceback.format_exception(T,V,TB))

class ChildProcessError(Exception) :
    pass

def process_launch_task_in_tmux(queue, task, gpu_index, filter_output=True):
    """Run command in tmux and monitor output"""
    def print_task(x):
        print('Task {}: {}'.format(task.name, x.rstrip('\n').strip()))

    def print_relevant_output_return_success_status(outlines):
        """Prints out relevant lines of output (LFADS specific)
        and return True if the graceful exit line is printed"""
        task_success = False
        line = next(outlines)
        while line is not None and line != '':
            ## These strings are LFADS specific
            if "Stopping optimization" in line:
                # LFADS trained successfully
                task_success = True
                print_task(line)

            # determine whether to print output
            elif not filter_output or "learning rate" in line:
                print_task(line)

            line = next(outlines)

        return task_success

    try:
        # touch outfile to ensure it exists
        touch_file(task.outfile)

        # postpend tee redirection appending to outfile
        # generate_tee_command can take a list of commands too and will combine with semicolons
        #append = True
        #command = generate_tee_command(task.command, task.outfile, append)

        # prepend gpu specification
        command = 'export CUDA_VISIBLE_DEVICES={}; {}'.format(gpu_index, task.command)

        # wrap in tmux session and wait for completion
        tmux_command = generate_tmux_command(task.name, command)
        #print(tmux_command)

        # if the tmux session is already running, throw an error, usually means
        # this model is already training
        if check_tmux_session_exists(task.name):
            print('Queue: Task {} is already running. Monitoring'.format(task.name))
            #raise ChildProcessError('Tmux session {} already exists.'.format(task.name))
            pid = None
            task.running_on_gpu = None

        else:
            # need to launch this task
            with mutex:
                #print(tmux_command)
                subp = subprocess.Popen(tmux_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        shell=True, universal_newlines=True);

            time.sleep(0.1)
            if not check_tmux_session_exists(task.name):
                # usually this means that the command immediately failed
                raise ChildProcessError('Tmux session immediately terminated running "{}" '.format(tmux_command))
            pid = subp.pid

        # notify master process we're starting
        queue.put(TaskStartedMessage(task.index, pid, task.name))

        # monitor outfile for output in real time and print when it comes in
        task_success = False
        with open(task.outfile, 'r') as outfile:
            outlines = follow_file(outfile)

#             while subp.poll() is None and check_tmux_session_exists(task.name):
            while check_tmux_session_exists(task.name):
                task_success = task_success or print_relevant_output_return_success_status(outlines)
                time.sleep(1)

            # check one last time after termination
            task_success = task_success or print_relevant_output_return_success_status(outlines)

        # Mark it finished so we don't need to redo next time
        if task_success:
            touch_file(task.donefile)

        # notify master process that this task is done
        queue.put(TaskCompletedMessage(task.index, task_success, get_tail_file(task.outfile, 10)))

    except Exception as e:
        # something went wrong, pass this along to main
        queue.put(TaskExceptionMessage(task.index))


def print_task_status_summary(tasks):
    num_finished = num_running = num_skipped = num_failed = 0

    for task in tasks:
        if task.skipped_donefile_exists:
            num_skipped += 1
        elif task.has_failed:
            num_failed += 1
        elif task.has_finished:
            num_finished += 1
        elif task.is_running():
            num_running += 1

    print('Queue: {0} skipped, {1} finished, {2} failed, {3} running'
          .format(num_skipped, num_finished, num_failed, num_running))

def run_command_in_tmux_session_no_monitoring(session_name, command):
    """Run command in tmux session in a detached way, no monitoring of output"""

    # wrap in tmux session and wait for completion
    tmux_command = generate_tmux_command(session_name, command)

    # if the tmux session is already running, throw an error, usually means
    # this model is already training
    if check_tmux_session_exists(session_name):
        print('Tmux session {} already exists'.format(session_name))
        raise ChildProcessError('Tmux session {} already exists.'.format(session_name))

    with mutex:
        subp = subprocess.Popen(tmux_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                shell=True, universal_newlines=True);
    return subp

def get_open_port():
    """Use bind(0) to get a free open port"""
    import socket
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.bind(("",0))
    port = s.getsockname()[1]
    s.close()
    return port

def launch_tensorboard_in_tmux(session_name, tensorboard_script, port):
    command = 'bash {} --port={}'.format(tensorboard_script, port)
    print(command)
    return run_command_in_tmux_session_no_monitoring(session_name, command)

def run_lfads_queue(queue_name, tensorboard_script_path, task_specs,
                    gpu_list=None, one_task_per_gpu=True, max_tasks_simultaneously=None, ignore_donefile=False):

    WAIT_TIME = 0.2

    if 'TMUX' in os.environ:
        print('Warning: tmux sessions will be nested inside the current session')
        del os.environ['TMUX']

    tasks = build_task_list(task_specs)

    gpu_status = query_gpu_status()

    if gpu_list:
        gpu_status = [gpu_status[i] for i in gpu_list]
    num_gpus = len(gpu_status)

    # compute number of tasks we can do simultaneously
    # factoring number of GPUs. If num_gpus == max_tasks_simultaneously, the
    # load balancer will implicitly place one task on each gpu
    num_cpus = cpu_count()
    if one_task_per_gpu:
        # setting max_tasks_simultaneously <= num_gpus ensures that no gpu will
        # ever have more than one task due to the scheduling algorithm
        if max_tasks_simultaneously is None:
            max_tasks_simultaneously = num_gpus
        else:
            max_tasks_simultaneously = min(max_tasks_simultaneously, num_gpus)
    elif max_tasks_simultaneously is None:
            max_tasks_simultaneously = num_cpus-1

    def print_status(x):
        print('Queue: ' + x.rstrip('\n'))

    # is tensorboard running in tmux?
    tensorboard_session_prefix = '{}_tensorboard'.format(queue_name)
    running_tensorboard_sessions = get_list_tmux_sessions_name_starts_with(tensorboard_session_prefix)

    if running_tensorboard_sessions:
        # tensorboard already running
        m = re.search('port(?P<port>\d+)', running_tensorboard_sessions[0])
        port = m.group('port') if m is not None else None
        print_status('TensorBoard already running on port {} in tmux session {}'.format(port, running_tensorboard_sessions[0]))
    else:
        # launch the tensorboard on an open port in a tmux session (if not already open)
        port = get_open_port()
        tensorboard_session = '{}_port{}'.format(tensorboard_session_prefix, port)
        print_status('Launching TensorBoard on port {} in tmux session {}'.format(port, tensorboard_session))
        launch_tensorboard_in_tmux(tensorboard_session, tensorboard_script_path, port)

    print_status('Initializing with {} GPUs and {} CPUs, max {} simultaneous tasks'
                 .format(len(gpu_status), num_cpus, max_tasks_simultaneously))

    # check for tasks already completed
    if not ignore_donefile:
        for task in tasks:
            task.mark_finished_if_donefile_exists()
            if task.skipped_donefile_exists:
                print('Task {}: skipping, task already completed'.format(task.name))

    # communication queue for each process
    message_queue = Queue(100)

    while not check_all_tasks_complete(tasks):

        # check queue for new messages
        do_status_summary =- False
        while message_queue.qsize() > 0:
            try:
                msg = message_queue.get_nowait()

                if type(msg) is TaskStartedMessage:
                    task = tasks[msg.task_index]

                    if msg.pid is not None:
                        # None means the task was already running previously, so don't print anything
                        print('Task {}: started in tmux session {} on GPU {} with PID {}'.format(task.name, msg.tmux_session, task.running_on_gpu, msg.pid))

                        # deduct from gpu memory
                        gpu = find_gpu_by_index(gpu_status, task.running_on_gpu)
                        gpu.memfree -= task.memory_req
                        gpu.incr_num_tasks()

                    sys.stdout.flush()

                elif type(msg) is TaskCompletedMessage:
                    task = tasks[msg.task_index]
                    if msg.success:
                        print('Task {}: completed successfully'.format(task.name))
                    else:
                        task.has_failed = True
                        if len(msg.tail) > 0:
                            print('Task {}: TERMINATED UNEXPECTEDLY. Final output:'.format(task.name))
                            print(msg.tail)
                        else:
                            print('Task {}: TERMINATED UNEXPECTEDLY with no output'.format(task.name))

                    task.has_finished = True
                    # return to available gpu memory
                    gpu = find_gpu_by_index(gpu_status, task.running_on_gpu)
                    gpu.memfree += task.memory_req
                    gpu.decr_num_tasks()
                    do_status_summary = True

                    sys.stdout.flush()

                elif type(msg) is TaskExceptionMessage:
                    task = tasks[msg.task_index]
                    task.has_finished = True
                    task.has_failed = True
                    print('Task {}: INTERNAL ERROR. Exception was:'.format(task.name))
                    print(msg.message)
                    do_status_summary = True

                    if task.running_on_gpu is not None:
                        gpu = find_gpu_by_index(gpu_status, task.running_on_gpu)
                        gpu.memfree += task.memory_req
                        gpu.decr_num_tasks()

                    sys.stdout.flush()

                else:
                    print('Unknown message {}'.format(msg))

            except Empty:
                pass

        # check again since tasks have now been marked complete
        if check_all_tasks_complete(tasks):
            break

        if do_status_summary:
            print_task_status_summary(tasks)

        # only run a certain number of tasks at the same time to avoid inefficient use of the CPUs
        if check_num_tasks_running(tasks) >= max_tasks_simultaneously:
            #print_status('Waiting for free CPU to become available')
            time.sleep(WAIT_TIME)
            continue;

        if check_all_tasks_completed_or_running(tasks):
            #print_status('All tasks launched or finished, waiting for last batch to complete')
            time.sleep(WAIT_TIME)
            continue;

        # find next task for which there is sufficient GPU memory
        (task_index, gpu_index) = pick_next_task(gpu_status, tasks)

        if task_index is None:
            #print_status('Waiting for GPU memory to become available')
            time.sleep(WAIT_TIME)
            continue;

        task = tasks[task_index]
        #print('Task {}: launching on gpu {}'.format(task.name, gpu_index))
        sys.stdout.flush()

        # mark task as running
        task.running_on_gpu = gpu_index

        # launch a process to monitor the task, also receive messages via Queue
        p = Process(target=process_launch_task_in_tmux, args=(message_queue, task, gpu_index, True))
        task.process = p
        p.start()

        time.sleep(WAIT_TIME)

    print_status('All tasks completed.')
    print_task_status_summary(tasks)

    # wait for all the subprocesses to complete, should be quick since all tasks are reported done now
    for task in tasks:
        if task.process is not None:
            task.process.join()

    message_queue.close()

    return tasks


__all__ = [run_lfads_queue, query_gpu_status, GpuStatus, Task, get_list_tmux_sessions, check_tmux_session_exists]
