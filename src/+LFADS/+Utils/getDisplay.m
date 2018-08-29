function display = getDisplay()
    str = getenv('DISPLAY');
    if ~isempty(str)
        [~, rem] = strtok(str, ':');
        if isempty(rem) || numel(rem) < 2
            display = [];
        else
            display = str2double(rem(2:end));
            if isnan(display)
                display = [];
            end
        end
    else
        display = [];
    end
end