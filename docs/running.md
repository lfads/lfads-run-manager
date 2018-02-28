# Running LFADS

To train the LFADS model using Python+Tensorflow, you need to generate shell scripts that will actually call `run_lfads.py` and do the work of training the model. `lfads-run-manager` provides two ways to go about this.

!!! warning "Add the `run_lfads.py` folder to your shell PATH"
    Be sure that the LFADS python source folder is on your shell path, such that running `which run_lfads.py` prints the directory where the Python+Tensorflow code LFADS is located. If not, you'll need to run something like `export PATH=$PATH:/path/to/models/research/lfads` and consider adding this to your `.bashrc` file.

    If Matlab is able to determine the location of `run_lfads.py` (meaning that it's own inherited `PATH` was set correctly), it will prepend an `export PATH=...` statement to each generated shell script for you. If not, you can try calling `setenv('PATH', '...')` from within Matlab to add `run_lfads.py` to the path. before generating the shell scripts.

    Alternatively, you can hard-code the location to `run_lfads.py` by passing along the fully specified path to each of the `writeShellScript...` methods as `'path_run_lfads_py', '/path/to/run_lfads.py'`

!!! tip "Virtualenv support"
    Each of the methods below supports a `:::matlab 'virtualenv', 'environmentName'` parameter-value argument. If specified, a `source activate environmentName` will be prepended to each script that calls Python for you. This is needed when Tensorflow is installed inside a virtual environment.

## Launching each run individually from shell scripts

It is possible to run each model individually, but you'll probably prefer to [queue everything at once](#lfads-queue-automatically-queueing-many-runs).

### Training the model
The first is to manually generate shell scripts for each run and then run them yourself. First, for each run `i`, you will call:

```matlab
rc.runs(i).writeShellScriptLFADSTrain('cuda_visible_devices', 0, 'display', 0);
```

Here, you should specify options that will be written into the shell script, the key ones being:

* `cuda_visible_devices` - which GPU index to run this model on, e.g. `0`. Use the `nvidia-smi` to enumerate the available GPUs on your system
* `display` - the X display to use, e.g. `0`, which will set `DISPLAY` to `:0`. The python code generates plots during training that will appear in TensorBoard. Generating these plots requires a display. When running in a remote server, you'll need to specify this, and possibly to launch an X server using something like `tightvnc` or `vncserver`.
* `appendPosteriorMeanSample` - `true` or `false` specifying whether to chain the posterior mean sampling operation after the training is finished. The default is `false`, but if you set this to `true`, you won't need to call `writeShellScriptPosteriorMeanSample` below.
* `appendWriteModelParams` - `true` or `false` specifying whether to chain the posterior mean sampling operation after the training is finished. The default is `false`, but if you set this to `true`, you won't need to call `writeShellScriptWriteModelParams` below.

This will generate an `lfads_train.sh` in the corresponding run's folder. For the first run in our example, this is at
```
~/lorenz_example/runs/exampleRun/param_Qr2PeG/single_dataset001/lfads_train.sh
```

The script essentially launches Python to run `run_lfads.py` with the specific parameters you've indicated in `RunParams` and pointing at the corresponding datasets, which were saved earlier when we called `rc.prepareForLFADS`.

```bash
#!/bin/bash

path_to_run_lfads=$(which run_lfads.py)
if [ ! -n "$path_to_run_lfads" ]; then
    echo "Error: run_lfads.py not found on PATH. Ensure you add LFADS to your system PATH."
    exit 1
fi

DISPLAY=:0 CUDA_VISIBLE_DEVICES=0 python $(which run_lfads.py) --data_dir=/home/djoshea/lorenz_example/runs/exampleSingleSession/param_YOs74u/single_dataset001/lfadsInput --data_filename_stem=lfads --lfads_save_dir=/home/djoshea/lorenz_example/runs/exampleSingleSession/param_YOs74u/single_dataset001/lfadsOutput --cell_clip_value=5.000000 --factors_dim=8 --ic_enc_dim=64 --ci_enc_dim=128 --gen_dim=64 --keep_prob=0.950000 --learning_rate_decay_factor=0.980000 --device=/gpu:0 --co_dim=0 --do_causal_controller=false --do_feed_factors_to_controller=true --feedback_factors_or_rates=factors --controller_input_lag=1 --do_train_readin=true --l2_gen_scale=500.000000 --l2_con_scale=500.000000 --batch_size=150 --kl_increase_steps=900 --l2_increase_steps=900 --ic_dim=64 --con_dim=128 --learning_rate_stop=0.001000 --temporal_spike_jitter_width=0 --allow_gpu_growth=true --kl_ic_weight=1.000000 --kl_co_weight=1.000000 --inject_ext_input_to_gen=false
```

Running the `lfads_train.sh` script will launch the Tensorflow training which will take some time. You likely want to launch this in a `tmux` session if running remotely.

### Sampling the posterior means
Next, generate the `lfads_posterior_mean_sample.sh` script to sample the posterior means, which can be launched after training has completed. If you set `appendPosteriorMeanSample` to `true` in `writeShellScriptLFADSTrain`, you can skip this step.

```matlab
rc.runs(i).writeShellScriptLFADSPosteriorMeanSample('cuda_visible_devices', 0);
```

### Writing the model parameters

Lastly, we want to export the trained model parameters to disk as an HD5 file. We do this by generating the shell script using

```matlab
rc.runs(i).writeShellScriptWriteModelParams('cuda_visible_devices', 0);
```

If you set `appendWriteModelParams` to `true` in `writeShellScriptLFADSTrain`, you can skip this step. These results will be written to a file called `lfadsOutput/model_params`, though these results can be loaded into Matlab using `run.loadModelTrainedParams()`.


### Launching Tensorboard
You can monitor the progress of each run by generating a script that launches TensorBoard.
```matlab
rc.writeTensorboardShellScript();
```

This will create `launch_tensorboard.sh` which will launch Tensorboard which can then be visited at `http://localhost:PORT`.

## LFADS Queue: Automatically queueing many runs

Manually running each of these shell scripts in sequence can be tedious, especially if you don't have enough GPUs or CPUs to run them all in parallel and individual runs take hours or days to complete. To make this part of the process more complete, you can alternatively use the Python task queueing system which will take care of training all the LFADS models for you.

!!! warning "Only supported on Linux"
    Unfortunately, this task queueing system is not supported on Mac OS at the moment, primarily because it depends on `nvidia-smi`, though it's theoretically possible with `cuda-smi` with light code changes. However, Tensorflow has discontinued explicit GPU support on Mac OS anyway. This has also never been tested on Windows, as you'd need to get `tmux` working.

First, we'll generate the Python script from Matlab that enumerates all of the runs:

```matlab
rc.writeShellScriptRunQueue('display', 500, 'gpuList', [0 1 2 3]);
```

The first argument `display` specifies the X11 display for plotting as before. `gpuList` enumerates the indices of GPUs that can be used for the runs. This argument is optional if all of the GPUs are viable for Tensorflow on your system.

!!! note "Capping the number of simultaneous runs"
    You can also manually specify `maxTasksSimultaneously` if you wish to cap the number of simultaneous runs. By default this is set to the minimum of the number of CPUs on your system and the available GPU memory. By default, each LFADS task is assumed to use 2000 MB of GPU memory, but you can adjust this by specifying `gpuMemoryRequired`.

This will generate a Python script `run_lfads.py`, which for our example lives here:

```bash
~/lorenz_exajjmple/runs/exampleRun/run_lfadsqueue.py
```

!!! note "lfads-run-manager repo folder will be added to your PYTHONPATH automatically"
    The `run_lfadsqueue.py` script depends on `lfadsqueue.py`, which lives in the root of the `lfads-run-manager` repository. This will be added to your PYTHONPATH environment variable in the `run_lfadsqueue.py` script.

!!! warning "Install and configure `tmux`"
    The LFADS queue launches each LFADS run inside its own `tmux` session to make it easy to monitor the runs as they are running. You'll need to install `tmux`.

    Also, `tmux` is finnicky about environment variables, which are only loaded when the `tmux` server first launches, not when a new session is started. The main one you need is that `run_lfads.py` must be on your `PATH` somewhere. If Matlab is able to determine this location (meaning that it's own inherited `PATH` was set correctly), it will prepend an `export PATH=...` statement to each `lfads_train.sh` script for you. If not, you can try calling `setenv('PATH', '...')` from within Matlab to add `run_lfads.py` to the path. before generating the shell scripts.

    If you're having trouble, you might want to launch a new `tmux` session using:

    ```bash
    tmux new-session
    ```

    Then from inside `tmux`, test that `which run_lfads.py` prints a location and that you are able to launch python and run `import tensorflow as tf` without any issues.

You can then kick everything off by running `python run_lfadsqueue.py` at the command line. It's recommended to do this from inside your own `tmux` session if you're running on a remote server, so you can monitor the task runner.

!!! tip "Python virtual environments"

    If tensorflow is installed in a Python virtual environment, you can have this environment be automatically `source activate`d within the training scripts using:
    ```matlab
    rc.writeShellScriptRunQueue('display', 500, 'gpuList', [0 1 2 3], 'virtualenv', 'tensorflow');
    ```

A few notes on how the system works:

* Output from Python will be `tee`'d into `lfads.out`, so you can check the output during or afterwards either there or in the `tmux` session.
* When a model finishes training and posterior mean sampling, a file called `lfads.done` will be created
* If the task runner detects an `lfads.done` file, it will skip that run. Unless you pass `:::matlab 'rerun', true` to `writeShellScriptRunQueue`, in which case every run will be rerun. This is convenient if you've added additional runs and just want the new ones to run.
* If a run fails, the error will be printed by the task runner and `lfads.done` will not be created
* A running tally of how many runs are currently running, have finished, or have failed will be printed
* You can enter a run's `tmux` session directly to monitor it. The list of sessions can be obtained using `tmux list-sessions`. You can also abort it using `Ctrl-C` and it will be marked as failed by the task runner.

The `run_lfadsqueue.py` script will periodically output updates about how the runs are proceeding:

```bash
(tensorflow) ➜  exampleRun python run_lfadsqueue.py
Warning: tmux sessions will be nested inside the current session
Queue: Launching TensorBoard on port 42561 in tmux session exampleRun_tensorboard_port42561
bash /home/djoshea/lorenz_example/runs/exampleRun/launch_tensorboard.sh --port=42561
Queue: Initializing with 2 GPUs and 12 CPUs, max 4 simultaneous tasks
Task lfads_param_Qr2PeG__single_dataset001: launching on gpu 0
Task lfads_param_Qr2PeG__single_dataset001: started in tmux session lfads_param_Qr2PeG__single_dataset001 on GPU 0 with PID 19498
Task lfads_param_Qr2PeG__single_dataset002: launching on gpu 1
Task lfads_param_Qr2PeG__single_dataset002: started in tmux session lfads_param_Qr2PeG__single_dataset002 on GPU 1 with PID 19527
Task lfads_param_Qr2PeG__single_dataset003: launching on gpu 0
Task lfads_param_Qr2PeG__single_dataset003: started in tmux session lfads_param_Qr2PeG__single_dataset003 on GPU 0 with PID 19551
Task lfads_param_Qr2PeG__all: launching on gpu 1
Task lfads_param_Qr2PeG__all: started in tmux session lfads_param_Qr2PeG__all on GPU 1 with PID 19585
Task lfads_param_Qr2PeG__single_dataset003:      Decreasing learning rate to 0.009800.
Task lfads_param_Qr2PeG__single_dataset001:      Decreasing learning rate to 0.009800.
Task lfads_param_Qr2PeG__single_dataset001:      Decreasing learning rate to 0.009604.
Task lfads_param_Qr2PeG__single_dataset003:      Decreasing learning rate to 0.009604.
Task lfads_param_Qr2PeG__single_dataset003:      Decreasing learning rate to 0.009412.
Task lfads_param_Qr2PeG__single_dataset001:      Decreasing learning rate to 0.009412.
```

As the tasks run, the task queue will print out messages related to decreasing the learning rate, which is one way to measure ongonig progress towards the termination criterion (when the learning rate hits `c_learning_rate_stop`). When a task fails or completes, the queue will print out a running tally.

Note that TensorBoard has automatically been launched on an available port, here on `42561`. You can also directly attach to the tmux sessions whose names are indicated in the script as "Tasks", which can be listed using `tmux list-sessions`.

```bash
➜  exampleRun tmux list-sessions
matlab: 4 windows (created Tue Oct  3 21:51:49 2017) [201x114] (attached)
exampleRun_tensorboard_port42561: 1 windows (created Fri Oct  6 14:43:16 2017) [201x113]
lfads_param_Qr2PeG__all: 1 windows (created Fri Oct  6 14:43:17 2017) [201x113]
lfads_param_Qr2PeG__single_dataset001: 1 windows (created Fri Oct  6 14:43:16 2017) [201x114]
lfads_param_Qr2PeG__single_dataset002: 1 windows (created Fri Oct  6 14:43:16 2017) [201x113]
lfads_param_Qr2PeG__single_dataset003: 1 windows (created Fri Oct  6 14:43:17 2017) [201x113]
```

If you wish to abort ongoing runs, you can either attach to them directly and use `Ctrl-C`, or use `tmux kill-session SESSIONNAME`. When everything has completed, you'll see something like this:

```bash
Task lfads_param_Qr2PeG__all: Stopping optimization based on learning rate criteria.
Task lfads_param_Qr2PeG__all: completed successfully
Queue: All tasks completed.
Queue: 0 skipped, 4 finished, 0 failed, 0 running
```
