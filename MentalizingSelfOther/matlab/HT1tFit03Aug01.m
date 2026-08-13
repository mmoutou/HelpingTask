% HT1tFit03Aug01.m: Script, not function, to fit a toDo range of data, decisions only, 
% from 03Aug26 Prolific pts,Helping Task. See subesquent versions for conditional blocks
% to get best-fit Hessian and simulated data.

toDo = 26:30; % 1:2; % 3:8; 9:14; 15:20; 21:25; 26:30;

thisFitVer = [num2str(toDo(1)) 'to' num2str(toDo(end)) 'a'];
thisFitStr = ['HT1tFit' thisFitVer];    % to store outputs etc.
write2disk = 0;          % Don't write to disk if we are not sure it works ...
localDebug = 1;          % For code testing - uses single-point (pseudo)grid.

datDir3Aug = '~/Dropbox/FIL_aux/MSc_iBSc_PhD_student/Kwan_(Tyler)/HelpingTaskProject_sharing/ProlificHelpingTask/HelpingTask03Aug26/data/';
cwd = cd; cd(datDir3Aug);
load('dataDemogr3Aug.mat')   % load the Prolific demographic data, dataDemogr3Aug
load('allProlificTaskData03Aug.mat')  % load prolD, with all the p, d, v (initial) for each of 30 pts
resDir = '~/Dropbox/FIL_aux/MSc_iBSc_PhD_student/Kwan_(Tyler)/HelpingTaskProject_sharing/ProlificHelpingTask/fitMentSO/'; 
cd(resDir); 

psInit0 = []; 
psPr0  = []; 
fit={};  % to hold everything.
D={};    % to hold all the d
P={};    % to hold all the P

for ptN=toDo

  d = prolD{ptN}.d; 
  p = prolD{ptN}.p; 
  [ps, slp, sll, psPr, Hess] = HT1MAP01t( d, p, psPr0, psInit0, localDebug ); 

  d.feelm = NaN; d.feelu=NaN;   % Make it crystal clear these have not been fit.
  d.evalRL = NaN * d.evalRL ;   %   ... ditto.
  % Actually fitted:
  d.Spref    = ps.Spref;
  d.prevPri  = ps.prevPri;
  d.PartnPr  = ps.SPartnPr;
  d.prevp    = ps.prevp;
  d.prevu    = ps.prevu;
  d.Spartp   = ps.SPartp;
  d.Spartu   = ps.SPartu;
  d.T        = ps.T;
  d.blockLR  = ps.blockLR;
  d.tHess    = Hess;

  D{ptN} = d;  P{ptN}= p; 

  eval(['fit{' num2str(ptN) '}.ps=ps;']);
  eval(['fit{' num2str(ptN) '}.slp=slp;']);  
  eval(['fit{' num2str(ptN) '}.sll=sll;']);  
  eval(['fit{' num2str(ptN) '}.psPr=psPr;']);  
  eval(['fit{' num2str(ptN) '}.Hess=Hess;']);  

  if write2disk
    % Save after every fit:
    eval(['save(''' thisFitStr '.mat''',',''P'',''D'',''fit'');']);
  end

end


return;  %% ~~~~~~~~~~~~~~~~~~~ eof ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
