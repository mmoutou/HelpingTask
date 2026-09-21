function [ps, slp, sll, psPr, hessian, output] = HT1MAP02t( d, p, psPr, psInit, localdebug )
% [ps, slp, sll, psPr, hessian, output] = HT1MAP02t( d, p, psPr, psInit ), following
%   Orestis Zavlis work in HT1MAP01t to fix prevu and SPartu. Here instead given
%         informative priors to obtain well-conditioned Hessian, via psPr below.
%   Max A Posterori parameters for Helping Task, HT1, To be used by called by HT1fit*
%   Find optimal parameters ps for one participant for contribution decisions (=t) using HT1ll1
%         d is the expt. data, p are the parameters *of the task*, not the pt.,
%   psPr is the (possibly group-derived) prior, psInit is the optional
%         value of parameters ps to start exploration from.
%   Provide an empty psPr in order not to use prior at all. Otherwise an
%         uninformative prior will be used.
%   Relies on pslPrHT1(ps,psPr,p) for log-prior of params
%   slp is the sumlogposterior, ssl the sum log likelihood and psPr a 
%         copy of the priors on params.
%   hessian is the approx hessian at the optimum, 
%   output is some detail of the fit e.g. number of iterations etc.

% for debugging at the level of this function only.
try localdebug;  catch  localdebug = 1;  end

% Initialize ps - with 'neutral' values if not initial val provided: 
try  psInit.Spref; catch psInit.Spref=[]; end
if isempty(psInit.Spref)
  psInit.Spref=repmat([0,1,2,3],p.settingLevN);  % will start from a 'neutral' pref
  psInit.prevp = 0.1;  psInit.prevu=2;  % for noisyBino(0.5,2,4) for prevPri
  psInit.prevPri = noisyBino(psInit.prevp,psInit.prevu,p.Nl); 
  psInit.SPartp = 0.9; psInit.SPartu=1;  % ditto for SPartnPr
  psInit.SPartnPr = noisyBino(psInit.SPartp, psInit.SPartu, p.Nl);
  psInit.T = 0.2;
  psInit.blockLR = 0.15;
end
psbest = psInit; 

% Intialize coarse grid over continuous parameters : 
grid.prevp  = [0.05 0.45 0.85];
grid.prevu  = [-5 -1  1  5];
grid.SPartp = [0.2  0.5  0.8];
grid.SPartu = [0.1]   % [-5 -1  1 5];
grid.T      = [0.15 1.35 4.05];  
grid.blockLR  = [0.1 0.6];

if localdebug  % only one point in the grid
  grid.prevp  = [0.65]; %#ok<*NBRAK2>
  grid.prevu  = [1];
  grid.SPartp = [0.6];
  grid.SPartu = [1];
  grid.T      = [0.15];  
  grid.blockLR  = [0.15];
end

% Look for priors. Note that if psPr of NaN is given, we do 
% not fill it in, but set the log-prior component to zero below:
try  psPr;  catch  psPr=[]; end
if isempty(psPr)
  % Weak commonsense priors :
  psPr.Spref0 = nan(p.Nl,p.Nl,p.settingLevN); % Spref rows, Owrk cols, context pages, flat. 
  for o=1:p.Nl; for c=1:p.settingLevN; psPr.Spref0(:,o,c) = noisyBino(o/(1+p.Nl),50,4); end; end;
  psPr.prevp0 = [1.05, 1.05];  % A and B for betapdf for pSucc of noisyBino describing prevPri
  psPr.prevu0 = [2.0, 1.0];    % A and B for gampdf for U of noisyBino describing prevPri
  psPr.T0 = [1.5, 0.5];        % A and B for gampdf on T
  psPr.SPartp0 = [1.05, 1.05]; 
  psPr.SPartu0 = [2.0, 2.0];
  psPr.blockLR0 = [1.1, 1.4];  % A and B for betapdf for betaLR. Modestly discrourages very high apparent LRs.
else
    try
        if  isnan(psPr)
            % if NaN psPr (prior) is provided, set log-prior component to zero:
            nPrior = 0; 
        end
    catch
    end
end

%% Construct auxiliaries esp repertoire of basic almost-linear preferences:
pks = 0:(p.Nl-1);  % possible values of Spref
% Form some linear preferences of peaks, which will
% form the base for search:
basepatN = p.Nl^2;
basepref = zeros(basepatN,p.Nl);
basecnt=0;
for y1=1:p.Nl
  for y4 = 1:p.Nl
    basecnt=basecnt+1;
    basepref(basecnt,:) = y1 + (y4-y1)*((1:4)-1)/3;
  end
end
basepref = round(basepref) - 1;  % as actually all go 0 to p.Nl-1

gridit=0; 
gridtot = length(grid.prevp)*length(grid.prevu)*length(grid.SPartp)*...
          length(grid.SPartu)*length(grid.T)*length(grid.blockLR);

tstart=tic;

%% Main loop over parameter grid -- NB the gradient descent part follows this!
for i1=1:length(grid.prevp)
  for i2 = 1:length(grid.prevu)
    for i3 = 1:length(grid.SPartp)
      for i4 = 1:length(grid.SPartu)
        for i5 = 1:length(grid.T)
          for i6 = 1:length(grid.blockLR)
              psInit.prevp = grid.prevp(i1);  
              psInit.prevu = grid.prevu(i2);  % for noisyBino(0.5,2,4) for prevPri
              psInit.prevPri = noisyBino(psInit.prevp,psInit.prevu,p.Nl);
              psInit.SPartp = grid.SPartp(i3); 
              psInit.SPartu=  grid.SPartu(i4);  % ditto for SPartnPr
              psInit.SPartnPr = noisyBino(psInit.SPartp, psInit.SPartu, p.Nl);
              psInit.T = grid.T(i5);
              psInit.blockLR = grid.blockLR(i6);
              
              ps = psInit;

disp('Loop over base preference patterns with init. prevp:');
disp(ps);
%% Keeping one context prefs constant, look for other prefs:
for context=1:p.settingLevN
  bestbslp = -inf;  % start afresh for different contexts ...
  for pattn = 1:basepatN
    ps.Spref(context,:) = basepref(pattn,:);
    % Check that prior is provided and calc log prior of params:
    if ~isempty(psPr); lnPrior = pslPrHT1(ps,psPr,p); end
    % Debug line to provide detailed output acc. to setting in p, e.g. in fitHT1a :
    try RunType.detailed = p.detailed; catch RunType.detailed =[]; end;
    newslp = HT1ll1(ps, d, p, RunType) + lnPrior;
    if newslp > bestbslp; bestbpatt = pattn; bestbslp=newslp; end;
  end
  ps.Spref(context,:) = basepref(bestbpatt,:); 
end

%% Now scan around the 'base preference' patterns to see if can be improved:
disp('Now scanning around best base preference pattern so far ...');
%  First for Bobi prefs, just for funsies:
bestslp2 = bestbslp;
for context = fliplr(1:p.settingLevN)
  for Owrk=0:3
    for prf = pks(~(pks==ps.Spref(context,Owrk+1)))
       Sprbak = ps.Spref; % back up best so far
       ps.Spref(context,Owrk+1) = prf;
       % Check that prior is provided and calc log prior of params:
       if ~isempty(psPr); lnPrior = pslPrHT1(ps,psPr,p); end;
       newslp = HT1ll1(ps, d, p) + lnPrior;
       % Now the other way round - restore if no improvement!
       if newslp <= bestslp2
         ps.Spref=Sprbak; 
       else
         bestslp2 = newslp;
       end
    end
  end
end

gridit = gridit + 1;
disp(['grid iteration ' num2str(gridit) ' out of ' num2str(gridtot)]);
disp(['Best slp in this iteration = ' num2str(bestslp2) ' for Spref='])
disp(ps.Spref);
tel = toc(tstart)
 
% Store best ps and slp found so far in grid search
% if an improvement has been made :
if gridit==1   
  gridslp = bestslp2; 
  gridps = ps; 
else
  if bestslp2 > gridslp
    gridps = ps;
    gridslp = bestslp2;
  end
end

% close loops over all continuous param grids:
          end
        end
      end
    end
  end
end
% Set ps to the best one found in grid search:
ps = gridps;

%% Now keeping ps.Spref constant, fit the continuous
%  params:  prevp, prevu, SPartp, SPartu, T, blockLR
restpInit = [ps.prevp, ps.prevu, ps.SPartp, ps.SPartu, ps.T, ps.blockLR];

% Declare the function to be minimized by fmincon. The follwing has to be
% re-declared every time we want to do the fit, not e.g. before this loop :
details=0; 
mLP = @(restp)HT1lp2( restp, ps, d, p, psPr, details);
% boring: have to specify empty 'linear constraints' in order to get to 
% the arguments for the lower and upper bounds, acc. to the doc fmincon example ...
A = []; b = []; Aeq = [];  beq = [];
%    prevp, prevu, SPartp, SPartu, T,    blockLR
lb = [0.01,  -50,   0.01,   0.1,  0.001, 0.01];  % lower bounds for restp
ub = [0.99,   50,   0.99,    50,   100,  0.99];  % upper bounds for same
% was:
% lb = [0.01, 0.1,0.01, 0.1,0.001];  % lower bounds for restp
% ub = [0.99, 50, 0.99, 50, 100  ];  % upper bounds for same
disp('Now running fmincon for prevp, prevu, SPartp, SPartu, T ...');
[restpOpt, mmLL, ~, output, ~, ~, hessian] = ...
    fmincon(mLP, restpInit, A, b, Aeq, beq,lb,ub );

% store and display best (so far ...) :
ps.prevp =restpOpt(1); ps.prevu =restpOpt(2); 
ps.prevPri=noisyBino(ps.prevp,ps.prevu,p.Nl);
ps.SPartp=restpOpt(3); ps.SPartu=restpOpt(4); 
ps.SPartnPr=noisyBino(ps.SPartp,ps.SPartu,p.Nl);
ps.T = restpOpt(5);  % may be deliberately spewed out! 
ps.blockLR = restpOpt(6);

disp(['end of map fitting. Best ps ']);
disp(ps);
slp = -mmLL;           % deliberately spewed out!
disp([' final slp: ' num2str(slp)]);

sll = HT1ll1(ps, d, p);
disp([' final sll: ' num2str(sll)]);

return;

