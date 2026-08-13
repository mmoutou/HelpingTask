function [v, d, p] = csv2pdHT1(fileN,data_dir,HTdir)
% [v, d, p] = CSV2PDHT1(FILEN,DATA_DIR,HTDIR) convert Jae's csv to structures p and d
% csv2pdHT1 also initializes the variable v. All these are used in HT1ll1 and the like.
%   fileN is e.g. helpTask2025_PARTICIPANT_SESSION_2026-07-02_16h42.46.468.csv 
%   The rest are optional -- where to find and store the data, and the directory
%   with the matlab code HTdir.
% Demo at the very end of this file.

fs = filesep(); 
homeD = [getenv('HOMEDRIVE') getenv('HOMEPATH') fs];
datT = HT1json2mat(fileN); 

%%  Fill in p ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p.helperNum = max(datT{:,"Block"});
p.code = datT.ProlificId{end}; % THIS IS PROLIFIC ID FROM PAVLOVIA, NOT FROM PROLIFIC DEMOGRAPHICS!
p.PID = datT.PID{end}; 
try
    p.HTdir = HTdir;
catch
    p.HTdir = []; % main subdir with code.
end
if isempty(p.HTdir); p.HTdir = [homeD 'Dropbox' fs 'task_code' fs 'IPD' fs]; end  % main subdir with code.
try
    p.data_dir = data_dir;
catch
    p.data_dir = [];
end
if isempty(p.data_dir)
    p.data_dir = [homeD 'Dropbox' fs 'FIL_aux' fs 'MSc_iBSc_PhD_student' fs 'Kwan_(Tyler)' fs 'HelpingTaskProject_sharing' fs 'sandpit' fs];
    warning(['data_dir not given, set to ' p.data_dir]);
end

% estimate number of blocks (helpers), their trial number and type.
p.helperType = nan(p.helperNum,2);
for blN = 1:p.helperNum
    blInd = find(datT{:,"Block"}==blN);
    blLen = sum(datT{:,"Block"} == blN);
    if  isnan(blLen)
        error(['Failed at block ' num2str(blN)]); 
    else
        p.helperType(blN,2) = blLen;
    end
    if strcmp(datT{blInd(1),'Helper Policy'},"UNCARING")
        p.helperType(blN,1) = 1;
    elseif strcmp(datT{blInd(1),'Helper Policy'},"CARING")
        p.helperType(blN,1) = 2;
    else
        error(['Not ready for helperType= ', datT{blInd(1),'Helper Policy'}])
    end

end

p = prepParHT1(p); 


%%  Fill in d ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d.taskDesign = datT{1,'Task Design'}; 
[v,d,p] = prepDatHT1(p,d); 
TRow = cell(1,p.helperNum);    % To hold the row numbers for each helper in datT
              % as they are all in one big list.
% For clarity, get numbers of columns and check that they really correspond
SpolCol = 4;  if ~strcmp(d.evo_hd{SpolCol}, 'Spol');  error('Wrong Spol column');  end
sWorkCol = 8; if ~strcmp(d.evo_hd{sWorkCol},'sWork'); error('Wrong sWork column'); end
OpolCol = 20; if ~strcmp(d.evo_hd{OpolCol}, 'Opol');  error('Wrong Opol column');  end
oWorkCol= 21; if ~strcmp(d.evo_hd{oWorkCol},'oWork'); error('Wrong oWork column'); end
polDecRTCol=7;if ~strcmp(d.evo_hd{polDecRTCol},'RT'); error('Wrong polDecRT column'); end

retSCol=9;    if ~strcmp(d.evo_hd{retSCol},'retS'); error('Wrong retS column'); end
retOCol=10;   if ~strcmp(d.evo_hd{retOCol},'retO'); error('Wrong retO column'); end

SevalCol = 11; if ~strcmp(d.evo_hd{SevalCol},'Seval'); error('Wrong Seval column'); end
OevalCol = 12; if ~strcmp(d.evo_hd{OevalCol},'Oeval'); error('Wrong Oeval column'); end
RTSatisfCol=26;if ~strcmp(d.evo_hd{RTSatisfCol},'RTSeval'); error('Wrong RTSatisf column'); end
RTTrustCol=29; if ~strcmp(d.evo_hd{RTTrustCol},'RTOeval'); error('Wrong RTTrust column'); end
% end of pedantic checking! Now keep a record of which columns will have useful, real data:
d.HT1DatCol = [1:3,SpolCol,polDecRTCol,sWorkCol,oWorkCol,...
                retSCol,retOCol,SevalCol,OevalCol,RTSatisfCol,RTTrustCol ];

for helperN=1:p.helperNum
    TRow{helperN} = find(datT{:,"Block"} == helperN); 

    % Before-outcomes row, tagged 901 in col 2 of d.evo{helperN}
    d.evo{helperN}(p.polRow{helperN},SpolCol) = datT{TRow{helperN},"Effort Intended"}; 
    d.evo{helperN}(p.polRow{helperN},OpolCol) = NaN;   % The Helper policy should not matter.
    
    % Observations row, on the same row as the exchange (rather than the decisions or evaluations): 
    d.evo{helperN}(p.evoSt{helperN},sWorkCol) = datT{TRow{helperN},"Effort Actual"}; 
    d.evo{helperN}(p.evoSt{helperN},oWorkCol) = datT{TRow{helperN},"Partner Effort"}; 
    d.evo{helperN}(p.evoSt{helperN},retSCol) = datT{TRow{helperN},"Points For Help-Seeker"}; 
    d.evo{helperN}(p.evoSt{helperN},retOCol) = datT{TRow{helperN},"Points For Helper"}; 

    % Post-outcomes row, tagged 902 in col 2 of d.evo{helperN}
    d.evo{helperN}(p.QRow{helperN},SevalCol) = datT{TRow{helperN},"Satisfaction Rating"}; 
    d.evo{helperN}(p.QRow{helperN},OevalCol)= datT{TRow{helperN},"Trust Rating"};         
    d.evo{helperN}(p.QRow{helperN},RTTrustCol) = datT{TRow{helperN},"RT Trust Rating (ms)"}; 
    d.evo{helperN}(p.QRow{helperN},RTSatisfCol)= datT{TRow{helperN},"RT Satisfaction Rating (ms)"};     
end

% Some default values -- deliberately set to unlikely but should-compute values
% Be careful not to take them as if actually elicited!
d.Spref    = [0 1 2 2];     % Peaks of one's own preferred response for level of contrib. of Other
d.prevPri= [0.2 0.25  0.3  0.25];  % Prior over prevalance of location of preference peak of Others.   
d.PartnPr= [0.05 0.05 0.45 0.45] ; % participant's preferences about the type of Other
                                   % they would like to deal with.
d.T = 0.1; 
d.blockLR = 0.99;
% The next 3 are rule-of-thumb values from some of Meera's project -- see e.g.
% /home/michael/Dropbox/BASOR/BASOR_output/RiverCrossing/rc4fmri/main/beh_phys/ic2behan/Gosalia_et_al_result_copies/JointPref/magnanimToApr2024/fitIC2Sf.csv
d.evalRL =  [-1.7346 -2.1972 -1.3863 -0.6931 0.3300 0.6700 0.6700 0 -0.8473];
d.feelm  =  5;  
d.feelu  =  3;
warning('Param for pS set to defaults, e.g. e.g. d.Spref=[0 1 2 2] - THESE ARE NOT experimental data!');

%% End piece - tidy up and write copies ----------------------                
try
    save([p.data_dir p.code '_paramHT1'],'p'); 
    save([p.data_dir p.code '_datHT1'],  'd'); 
catch
    warning(['Could not save p or d to ' p.data_dir p.code '_paramHT1.mat'])
end
try 
    save([p.commDir fs [p.code '_paramHT1'] '.mat'],'p'); 
    save([p.commDir fs [p.code '_datHT1'] '.mat'],  'd'); 
catch
    warning('Could not save params. to p.commDir'); 
end

return; % end of function ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Demo:
% fileN =  'helpTask2025_PARTICIPANT_SESSION_2026-07-02_16h42.46.468.csv';
% fs = filesep; homeD = [getenv('HOMEDRIVE') getenv('HOMEPATH') fs];
% codeD = [homeD 'Dropbox' fs 'task_code' fs 'IPD' fs]; 
% sandpitD = [codeD 'sandpit' fs];
% cd(sandpitD);
% datD = [homeD 'Dropbox' fs 'FIL_aux' fs 'MSc_iBSc_PhD_student' fs 'Kwan_(Tyler)' fs 'HelpingTaskProject_sharing' fs 'sandpit' fs];

% [v, d, p] = csv2pdHT1(fileN)

% ---------------------------------------------------------------


