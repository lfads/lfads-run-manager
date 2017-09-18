function [ outPath, projPath] = mkProjectPath(runsRootPath, projectName, runName, dateTimeVec)
% function [ outPath] = mkProjectPath(runsRootPath, projectName, dateTimeVec)
% create a path based on the project name with added date/time 
% runRootPath, top path to where the lfads runs are stored
% projectName, name of the project for this lfads run
% dataTimeVec (optional), specific datetime vector [YY MM DD HH MM], if not
% provided current date and time are used
    if ~exist('dateTimeVec', 'var')
        dateTimeVec = datevec(now);
    end
    
    dateTimeVec = dateTimeVec(1:5);   % discard seconds
    projPath = sprintf('%02.0f%02.0f%02.0f_%02.0f%02.0f_%s', mod(dateTimeVec(1),100), ...
        dateTimeVec(2), dateTimeVec(3), dateTimeVec(4), dateTimeVec(5), runName);
    outPath = fullfile(runsRootPath, projectName, projPath);
    

end

