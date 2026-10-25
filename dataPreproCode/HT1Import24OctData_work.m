

%%  Working on converting from matlab table to p, d, using HT1csv2pd03.m
datDir24Oct = '/home/michael/Dropbox/FIL_aux/MSc_iBSc_PhD_student/Kwan_(Tyler)/HelpingTaskProject_sharing/ProlificHelpingTask/HelpingTask24Oct26/data/';
cwd = cd; cd(datDir24Oct);
taskDatFList = dir('helpTask2025*.csv');    %  raw json data
load('dataDemogr24Oct.mat')                  %  load the Prolific demographic data, dataDemogr24Oct
% 
datFTot = length(taskDatFList);
prolD = {}; 
% **************  KEY CALLING OF HT1csv2pd03: **********************
for fN = 1:datFTot
    [v, d, p] = HT1csv2pd03(taskDatFList(fN).name,datDir24Oct);
    ro = strcmp(dataDemogr24Oct{:,"ParticipantId"},d.code);
    d.prolDemog = dataDemogr24Oct(ro,:);
    p.data_dir = datDir24Oct;
    p.fileName = taskDatFList(fN).name;

    prolD{fN}.p = p; 
    prolD{fN}.d = d;
    prolD{fN}.v = v;
end
% Check the prolific IDs the other way around:
totProlN = size(dataDemogr24Oct,1);
noTaskData = {}; 
for ptN = 1:39 % totProlN
    found = 0;
    prolDemogrID = dataDemogr24Oct{ptN,'ParticipantId'};
    for fN = 1:datFTot
        if strcmp(prolD{fN}.d.code, prolDemogrID)
            found = 1;
        end
    end
    if ~found
        noTaskData{end+1} = prolDemogrID{1}; 
        disp(['Pt in row ',num2str(ptN), 'of dataDemogr* not found:'])
        disp(dataDemogr24Oct(ptN,:));
    end
end
save("allProlificTaskData24Oct.mat","prolD","noTaskData")



%% eof ------------------------------------------------------------------------

