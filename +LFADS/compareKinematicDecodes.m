function compareKinematicDecodes(seq, resultSet, varargin)

p = inputParser();
p.addParameter('colormap', distinguishable_colors(numel(resultSet), {'k', 'w'}), @(x) true); % ms
p.addParameter('separateConditions', false, @islogical);
p.parse(varargin{:});


% which condition is each trial
[cnames, ~, cond] = unique({seq.targetDirectionName});
nC = numel(cnames);

pydim = 2;

cmap = TrialDataUtilities.Color.toRGB(p.Results.colormap);
resultSet = wrapCell(resultSet);
h = gobjects(numel(seq), numel(resultSet)+1);
hAxes = gobjects(nC, 1);
clf;

for icond = 1:nC
    ttp = find(cond==icond);
    for itr = 1:numel(ttp)
        ntr = ttp(itr);

        T = resultSet{1}.T;
        
        % make the figure for raw spks decoding
        if p.Results.separateConditions
            hAxes(icond) = subtightplot(pydim, ceil(nC / pydim), icond);
        end
        axis off;
        h(ntr, 1) = plot(T(ntr).X(1,:), T(ntr).X(2,:), '-', 'color', [0.25 0.25 0.25]);
        hold on;
        
        if p.Results.separateConditions
            ylim([-1 1]*160);
            xlim([-1 1]*160);
        end
            
        set(gca,'xtick',[], 'ytick', []);

        for iR = 1:numel(resultSet)
            Tdec = resultSet{iR}.T_dec;
            h(ntr, iR+1) = plot(Tdec(ntr).xk(1,:), Tdec(ntr).xk(2,:), 'Color', cmap(iR, :));
        end     
% 
%         show_legend = false;
%         if show_legend && itr == (pydim-1)*pxdim+1
%             xlabel('x pos');
%             ylabel('y pos');
%             
%             legend({'true cursor', 'raw neural + vkf', ['lfads ' ...
%                                 'facts + vkf']},'location','best');
%         end
        %set(gca,'box','off');
        %set(gca,'xtick',[-150 150], 'ytick', [-150 150]);
    end
end
%equalize_axes(ah);

if p.Results.separateConditions
    linkaxes(hAxes, 'xy');
else
    axis tight;
end
axis off;

TrialDataUtilities.Plotting.setLineOpacity(h, 0.4);

% figure(1)
% set(gcf,'PaperOrientation', 'landscape')
% set(gcf,'paperposition', [0 0 10.5, 10.5/13*7])
% set(gcf,'papersize', [10.5 5.7])
% print('-dpdf', fullfile(outdir1, 'decoded_trajectories_raw'))
% 
% figure(2)
% set(gcf,'PaperOrientation', 'landscape')
% set(gcf,'paperposition', [0 0 10.5, 10.5/13*7])
% set(gcf,'papersize', [10.5 5.7])
% print('-dpdf', fullfile(outdir1, 'decoded_trajectories_lfads'))
