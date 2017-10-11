function showInLegend(h, names)
% showInLegend(h, names)
% show object h in legend with specified by default
% use legend(axh, 'show') to activate default legend
% names is either char (for scalar h), or cellstr

    if nargin >= 2 && ~isempty(names) && ischar(names)
        names = repmat({names}, numel(h), 1);
    end

    for i = 1:numel(h)
        if ~TrialDataUtilities.Plotting.isGraphicsHandle(h(i)), continue; end
        ann = get(h(i), 'Annotation');
        leg = get(ann, 'LegendInformation');
        set(leg, 'IconDisplayStyle', 'on');
        
        if ~isempty(names)
            set(h(i), 'DisplayName', names{i});
        end
    end

end