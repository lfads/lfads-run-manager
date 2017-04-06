function display = getDisplay()
    str = getenv('DISPLAY');
    if ~isempty(str)
        display = str2double(str(2:end));
    else
        display = [];
    end
end