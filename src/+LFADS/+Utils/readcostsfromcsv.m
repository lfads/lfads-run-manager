function [total_tr, total_val, recon_tr, recon_val, kl_tr, kl_val, l2, kl_w, l2_w] = readcostsfromcsv(csvlog_path)
    % total_tr: total cost (train)
    % total_val: total cost (validation)
    % recon_tr: reconstruction cost (train)
    % recon_val: reconstruction cost (validation)
    % kl_tr: KL cost (train)
    % kl_val: KL cost (validation)
    % l2: l2 cost
    % kl_w: KL weight
    % l2_w: l2 weight
    fid = fopen(csvlog_path);
    lfadslog = textscan(fid , '%s %f %s %f %s %f %f %s %f %f %s %f %f  %s %f %s %f %s %f' , 'Delimiter' , ',');
    fclose(fid);
    total_tr = lfadslog{6};
    total_val = lfadslog{7};
    recon_tr = lfadslog{9};
    recon_val = lfadslog{10};
    kl_tr = lfadslog{12};
    kl_val = lfadslog{13};
    l2 = lfadslog{15};
    kl_w = lfadslog{17};
    l2_w = lfadslog{19};
    


    
function aver= avglast10(vec)
    % Calc average of the last 10 elements
    if numel(vec>9)
        aver = mean(vec(end-9:end));
    else
        aver = mean(vec);
    end
        
