function outputString = teeify_string(inputString, teefile, append)
    % function outputString = tmuxify_string(inputString, tmux_session_name)

    if append
        outputString = sprintf('%s 2>&1 | tee -a %s', inputString, teefile);
    else 
        outputString = sprintf('%s 2>&1 | tee %s', inputString, teefile);
    end
  
end