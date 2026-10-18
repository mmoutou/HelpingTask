function where = where2findHT1(showWhere)
% function where2findHT1 provides structure 'where' with 
% various basics so that other code can look up e.g. 
% where.home or where.code etc. See ADJUST lines and modify!
    
try showWhere; catch showWhere = 1; end  % by default, print where we'll work as set here.

cwd = cd;         % just a record of where this is being run from.
fs = filesep();   % the character that separates folders from subfolders in the filesystem.
if ispc           % paths depend on the operating system, ie. Windoze vs. Linux/Mac
    homeDir = fullfile(getenv('HOMEDRIVE'), getenv('HOMEPATH'));
else
    homeDir = getenv('HOME');
end
try
    where.home = [homeDir fs];
    % Where Michael does coding : 
    where.code = [where.home 'localgit' fs 'HelpingTask' fs 'MentalizingSelfOther' fs 'matlab' fs];      
    % General purpose 'sandpit' for work: 
    where.sandpit = [where.home 'Dropbox' fs 'task_code' fs 'IPD' fs 'sandpit4HT' fs 'fitMentSO_wrk' fs];         
    where.results_July26 = [where.home 'Dropbox' fs 'FIL_aux' fs 'MSc_iBSc_PhD_student' fs 'Kwan_(Tyler)' ...
                            fs 'HelpingTaskProject_sharing' fs 'results_July26' fs];  % Results of Tyler's MSc data (Helping Task 1 Study 1)

    where.HT1Study2 = [where.home 'Dropbox' fs 'FIL_aux' fs 'MSc_iBSc_PhD_student' fs 'Kwan_(Tyler)' ...
                            fs 'HelpingTaskProject_sharing' fs 'ProlificHelpingTask' fs];  % Results of Tyler's MSc data (Helping Task 1 Study 1)
    % Can put nice results here for sharing:
    where.HT1St2MentSOres = [where.home 'localgit' fs 'HelpingTask' fs 'MentalizingSelfOther' fs 'results' fs];

    cd(where.sandpit);   % attempt to change to this directory to see if it's OK. 
    cd(cwd);
catch  % ADJUST LINE INDICATED BELOW TO SUIT YOUR OWN COMPUTER(S) :
    where.home = [homeDir fs];
    % ADJUST THE FOLLOWING LINES:
    warning('Attempted but could not cd to ''sandpit'' directory, using current working directory');
    where.sandpit = cwd;   
    cd(where.sandpit);
end


if showWhere
    disp(where);
end

return; % whole function.