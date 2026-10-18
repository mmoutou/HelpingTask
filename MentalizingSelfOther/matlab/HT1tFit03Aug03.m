% HT1tFit03Aug03.m: Script, not function, to fit a toDo range of data, decisions only, 
% from 03Aug26 Prolific pts, Helping Task. 
% version 03, i.e. ...03Aug03, Using weakly informative priors for 1/prevu and 1/SPartu,
%         and optimizer that does not depend on smoothness, Hessian etc.
%         See line thisFitVer = ...
% See subesquent versions for conditional blocks
% to get best-fit measures *and simulated data* .

toDo = 1:30;     % all participants. Was 26:30. localDebug uses only the first two.

thisFitVer = [num2str(toDo(1)) 'to' num2str(toDo(end)) '_03'];
thisFitStr = ['HT1tFit' thisFitVer];    % to store outputs etc.
write2disk = 1;          % Was 0, so it fitted and then threw the results away.
localDebug = 0;          % 0 = the full 864-point grid, i.e. the real MAP fit. Hours.
                         % 1 = a single grid point, minutes, useful only as a smoke test.
                         % Was 1 for the 25 Aug 26 shakedown run. Now the real thing.

if localDebug; toDo = toDo(1:2); end
Nl = 4; setLevN = 1;  % In HT1, have 4 levels of choice for effort and only 1 'context' (per Other)
psPr0.Spref0 = nan(Nl,Nl,setLevN); % Spref rows, Owrk cols, context pages, flat. 
for o=1:Nl
    for c=1:setLevN
        m = o/(1+Nl); u = 1.5*m*(1-m); 
        psPr0.Spref0(:,o,c) = discBetaMU(m,u,4); 
    end
end
psPr0.prevp0 = [1.05, 1.05];  % A and B for betalike for pSucc of noisyBino describing prevPri
psPr0.prevu0 = [1.5, 1]; 

psPr0.SPartp0 =  [1.05, 1.05]; % [45, 10];  for highly constrained diagnostic. Uninf. is: psPr0.SPartp0 = [1.05, 1.05]; 
psPr0.SPartu0 = [1.5 1]; % [2, 0.04] for diagnostic;  % check with x = 0:0.01:10; plot(x,gampdf(x,1.5,0.1)) % for less highly constrained.

psPr0.T0 = [1.5, 0.5];        % A and B for gamlike on T
psPr0.blockLR0 = [1.1, 1.4];  % A and B for betalike for betaLR. Modestly discrourages very high apparent LRs.

% psPr0  = [];  % If empty, HT1MAP02t below defaults to very weak priors. 

fs = filesep();  cwd = cd;
try
    dirs = where2findHT1;    % This is a function wherein you can add where data, outputs etc. are to be in your computer.
    datDir3Aug = dirs.sandpit; 
catch    
   % Paths derived from wherever this script lives 
   thisScript = which(mfilename);
   if isempty(thisScript); datDir3Aug = [cwd fs]; else; datDir3Aug = [fileparts(thisScript) fs]; end
end
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
fit={};  % to hold everything.
D={};    % to hold all the d
P={};    % to hold all the P

for ptN=toDo

  d = prolD{ptN}.d; 
  p = prolD{ptN}.p; 
  [ps, slp, sll, psPr, tFitMeasures] = HT1MAP03t( d, p, psPr0, psInit0, localDebug ); 

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
  d.sllt = sll;
  d.slpt = slp; 
  d.tFitMeasures = tFitMeasures;

  try d.tFixed = ps.tFixed; catch; d.tFixed = []; end

  D{ptN} = d;  P{ptN}= p; 

  eval(['fit{' num2str(ptN) '}.ps=ps;']);
  eval(['fit{' num2str(ptN) '}.slp=slp;']);  
  eval(['fit{' num2str(ptN) '}.sll=sll;']);  
  eval(['fit{' num2str(ptN) '}.psPr=psPr;']);  
  eval(['fit{' num2str(ptN) '}.tFitMeasures=tFitMeasures;']);  

  if write2disk
    % Save after every fit, under BOTH naming conventions, because
    % HT1fFit03Aug13bBoth wants d, p and tFit while this script builds D, P and fit.
    dSave = D;  pSave = P;  tFit = fit;
    save([resDir thisFitStr '.mat'],   'dSave','pSave','tFit','D','P','fit');
  end
  
  disp(['***************  pt ' num2str(ptN) ' finished with: **************']);
  disp('Fit measures:');
  disp(tFitMeasures);
  disp(['sllt: ' num2str(sll)]);
  disp('****************************************************');


end


disp(' ');
disp(['''t''reating each other fitting done, wrote ' resDir thisFitStr '.mat']);
disp('Next: run HT1fFit03Aug*Both . NB No Hessians and SEs for tFits !');
if localDebug
    disp('NB localDebug was 1, so a single grid point was used. Set it to 0 and rerun');
    disp('   for the real MAP fit, once you have seen the whole chain work.');
end


return;  %% ~~~~~~~~~~~~~~~~~~~ eof ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
