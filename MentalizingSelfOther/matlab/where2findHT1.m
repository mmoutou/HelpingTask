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
    where.code = [where.home 'Dropbox' fs 'task_code' fs 'IPD' fs];   % Where Michael does his local coding
    where.sandpit = [where.code 'sandpit' fs];                        % Where Michael puts rough results
    where.results_July26 = [where.home 'Dropbox' fs 'FIL_aux' fs 'MSc_iBSc_PhD_student' fs 'Kwan_(Tyler)' ...
                            fs 'HelpingTaskProject_sharing' fs 'results_July26' fs];  % Results of Tyler's MSc data (Helping Task 1 Study 1)

    where.HT1Study2 = [where.home 'Dropbox' fs 'FIL_aux' fs 'MSc_iBSc_PhD_student' fs 'Kwan_(Tyler)' ...
                            fs 'HelpingTaskProject_sharing' fs 'ProlificHelpingTask' fs];  % Results of Tyler's MSc data (Helping Task 1 Study 1)
    where.HT1St2MentSOres = [where.HT1Study2 'fitMentSO' fs]; % For MentalizingSelfOther (not TrustState HMM) 
                                                              %  model fits  & results.

    cd(where.HT1Study2);   % attempt to change to this directory to see if it's OK. 
    cd(cwd);
catch  % ADJUST LINE INDICATED BELOW TO SUIT YOUR OWN COMPUTER(S) :
    where.home = [homeDir fs];
    % ADJUST THE FOLLOWING LINES:


    cd(where.HT1Study2);   % attempt to change to this directory to see if it's OK. 
    cd(cwd);
end

addpath(where.code);     % working matlab code is here.

if showWhere
    disp(where);
end

return; % whole function.