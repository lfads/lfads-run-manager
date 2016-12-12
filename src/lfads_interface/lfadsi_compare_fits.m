function lfads_compare_two_fits(fs, ylims, labels, minEpoch)
% function lfads_compare_two_fits(fs, ylims, labels, minEpoch)

clf;

if ~exist('minEpoch','var')
    minEpoch = 1;
end

for nn = 1:numel(fs)
    if ~exist(fs{nn},'file')
        error(sprintf('Can''t find file: %s', fs{nn}));
    end

    y1 = lfadsi_read_fitlog(fs{nn});
    ep{nn} = str2double(y1(:,2));
    trcost{nn} = str2double(y1(:,9));
    vacost{nn} = str2double(y1(:,10));
    

    subplot(2,1,1)
    plot(ep{nn},trcost{nn});
    hold on;
    axis('tight');
    if exist('ylims','var') && numel(ylims)
        set(gca,'ylim',ylims);
    end
    ylabel('rec cost train');
    xlabel('epoch');

    if nn == numel(fs)
        if exist('labels','var') && numel(labels)
            legend(labels);
        end
    end

    subplot(2,1,2)
    plot(ep{nn},vacost{nn});
    hold on
    axis('tight');
    if exist('ylims','var') && numel(ylims)
        set(gca,'ylim',ylims);
    end
    ylabel('rec cost valid');
    xlabel('epoch');

    disp(sprintf('%g: epochs: %g, train: %g, valid: %g', nn, ep{nn}(end), min(trcost{nn}(minEpoch:end)), min(vacost{nn}(minEpoch:end))));
 
    if nn == numel(fs)
        if exist('labels','var') && numel(labels)
            legend(labels);
        end
    end
   
end



disp(' ');
