function outputString = tmuxify_string(inputString, tmux_session_name, varargin)
% function outputString = tmuxify_string(inputString, tmux_session_name)

p = inputParser;
p.addParameter('keepSessionAlive', true, @islogical);
p.parse(varargin{:});

if p.Results.keepSessionAlive
    outputString = sprintf(['CMD="%s"\nTMUX= tmux new-session -d -s "%s" "$CMD; read -n 1"\n'], ...
                           inputString, tmux_session_name);
else 
    outputString = sprintf('CMD="%s"\nTMUX= tmux new-session -d -s "%s" "$CMD"\n', ...
                           inputString, tmux_session_name);
end
                       
%outputString = sprintf(['CMD="%s"\nTMUX= tmux new-session -d -s "%s" "echo ''$CMD'';$CMD; read -n 1"'], ...
%                       inputString, tmux_session_name);
