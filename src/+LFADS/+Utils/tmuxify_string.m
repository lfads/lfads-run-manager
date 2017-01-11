function outputString = tmuxify_string(inputString, tmux_session_name)
% function outputString = tmuxify_string(inputString, tmux_session_name)

outputString = sprintf(['CMD="%s"\nTMUX= tmux new-session -d -s "%s" "$CMD; read -n 1"\n'], ...
                       inputString, tmux_session_name);

%outputString = sprintf(['CMD="%s"\nTMUX= tmux new-session -d -s "%s" "echo ''$CMD'';$CMD; read -n 1"'], ...
%                       inputString, tmux_session_name);
