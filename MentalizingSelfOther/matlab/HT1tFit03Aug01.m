% HT1tFit03Aug01.m: Script, not function, to fit a toDo range of data, decisions only, 
% from 03Aug26 Prolific pts,Helping Task. See subesquent versions for conditional blocks
% to get best-fit Hessian and simulated data.

toDo = 1:30;     % all participants. Was 26:30.

thisFitVer = [num2str(toDo(1)) 'to' num2str(toDo(end)) 'a'];
thisFitStr = ['HT1tFit' thisFitVer];    % to store outputs etc.
write2disk = 1;          % Was 0, so it fitted and then threw the results away.
localDebug = 0;          % 0 = the full 864-point grid, i.e. the real MAP fit. Hours.
                         % 1 = a single grid point, minutes, useful only as a smoke test.
                         % Was 1 for the 25 Aug 26 shakedown run. Now the real thing.

% Paths derived from wherever this script lives, so nothing needs editing and there is
% no cd into somebody else's Dropbox. Was two hard-coded ~/Dropbox/... lines plus cd.
fs = filesep();  cwd = cd;
thisScript = which(mfilename);
if isempty(thisScript); datDir3Aug = [cwd fs]; else; datDir3Aug = [fileparts(thisScript) fs]; end
resDir = datDir3Aug;
addpath(genpath(datDir3Aug));

% Demographics are not needed for fitting, so load them only if the file is there:
if exist([datDir3Aug 'dataDemogr3Aug.mat'],'file')
    load([datDir3Aug 'dataDemogr3Aug.mat']);
end

% prolD holds p, d and v for each participant. Look beside the script, then in the
% current folder, then anywhere below:
datFile = '';
for cand = {[datDir3Aug 'allProlificTaskData03Aug.mat'], [cwd fs 'allProlificTaskData03Aug.mat']}
    if exist(cand{1},'file'); datFile = cand{1}; break; end
end
if isempty(datFile)
    hits = dir([datDir3Aug '**' fs 'allProlificTaskData03Aug.mat']);
    if ~isempty(hits); datFile = fullfile(hits(1).folder, hits(1).name); end
end
if isempty(datFile)
    error(['allProlificTaskData03Aug.mat not found under %s . That file holds the ' ...
           'parsed task data and nothing can be fitted without it.'], datDir3Aug);
end
load(datFile);   % provides prolD
disp(['Loaded ' num2str(numel(prolD)) ' participants from ' datFile]);
toDo = toDo(toDo <= numel(prolD));

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
  try d.tFixed = ps.tFixed; catch; d.tFixed = []; end

  D{ptN} = d;  P{ptN}= p; 

  eval(['fit{' num2str(ptN) '}.ps=ps;']);
  eval(['fit{' num2str(ptN) '}.slp=slp;']);  
  eval(['fit{' num2str(ptN) '}.sll=sll;']);  
  eval(['fit{' num2str(ptN) '}.psPr=psPr;']);  
  eval(['fit{' num2str(ptN) '}.Hess=Hess;']);  

  if write2disk
    % Save after every fit, under BOTH naming conventions, because
    % HT1fFit03Aug13bBoth wants d, p and tFit while this script builds D, P and fit.
    dSave = D;  pSave = P;  tFit = fit;
    save([resDir 'HT1tFit1to30a.mat'], 'dSave','pSave','tFit');
    save([resDir thisFitStr '.mat'],   'dSave','pSave','tFit','D','P','fit');
  end

end


disp(' ');
disp(['Pass one done, wrote ' resDir 'HT1tFit1to30a.mat']);
disp('Next: run HT1fFit03Aug13bBoth');
if localDebug
    disp('NB localDebug was 1, so a single grid point was used. Set it to 0 and rerun');
    disp('   for the real MAP fit, once you have seen the whole chain work.');
end

return;  %% ~~~~~~~~~~~~~~~~~~~ eof ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
