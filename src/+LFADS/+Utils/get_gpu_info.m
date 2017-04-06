function [info, index_list] = get_gpu_info()

    fname = tempname;
    cmd = sprintf('nvidia-smi --query-gpu=index,name,utilization.gpu,memory.free,memory.total,display_mode --format=csv,nounits > %s', fname);
    [status, out] = system(cmd);
    if status
        error('Error running nvidia-smi command: %s', out);
    end
    info = readtable(fname);
    delete(fname);
    
    index_list = info.index;
    
end