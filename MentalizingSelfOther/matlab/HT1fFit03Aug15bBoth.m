% script HT1fFit03Aug15bBoth
% 
% To fit person-evaluation ratings, given the interactive behaviour (decision making)
% for the Prolific 03Aug26 Helping Task (HT1) dataset (which has no Other->Self ratings).
%
% Version 15  uses an exponential Returns kernel / autoregressive RL approach, and
%             includes a constant intercept E0 = [E0self,E0other]' as well as 
%             allowing Self Ret -> Oeval to be different than vice versa, and wOS22 
%             not constrained to 0 but freely fitted too.
%
% Improved or fixed Aug 26 (OZ), four things:
%  (a) evalRL(13) EvBlockLR now enters the likelihood, see llfeelHT1b. It used to be
%      fitted and reported without being read, so fmincon returned the starting value
%      untouched and the Hessian had an exactly zero row and column. That is why ALL
%      the standard errors came back NaN, not just this one's. Checked numerically:
%      sweeping it across its whole range left the summed log-likelihood identical to
%      ten decimal places and its finite-difference derivative was exactly zero.
%  (b) lambda2 is now set explicitly below, with a note on what it means, instead of
%      being left at a value that contradicted the model description.
%  (c) priors now come from priParEvalRLHT1, which has 13 columns, so parameter 13 is
%      regularised. priParEvalRL returns only 12.
%  (d) a flatness check runs before fmincon and refuses to fit any parameter that
%      leaves the objective unchanged, so this class of bug cannot recur silently.
% Text-search '(Change below for new model version)' and also check HT1lp02f and HT1ll1
% 
%   Fitting starts from fits of the neuroecon decision-making e.g. of HT1MAP01t 
%   (t for treating each other) 
% To change if fitted to all data or just self data, etc.,
% the crucial variable is approvs2fit below. 
%
%  BIC for this fitting was derived fairly rigorously from iBIC approx.
%      For HT1, the small sample approximation isn't needed, so -ln(2*pi) dropped
%  k=parPtN;    % participant level number of parameters
%  c=gpPtN;     % group-level number of parameters
%  N            % 
%  n            % number of datapoints per participant, summed over all blocks.
%  nN           % 
%  BIC_kc := -2 ln Lopt + k*ln(n/(2*pi) + (c/N)* ln (nN)
%                             = -2 ln Lopt + k*(ln n ) + (c/N)* ln (nN)

clear variables; 
%% ~~~~~~~~~~  Menu-like items ~~~~~~~~~~~~~~~~~~
apprStr = '11';   approvs2fit = [1,1];    % weights to place on the 
    % log lik for self, other, and poss. other->self in summing mslf2. 
    % if '10', [1,0] fit only the self data etc. NOT CHANGED IN MOST MODELS,
    % WHICH DEALS WITH THE SELF AND OTHER EVALUATIONS ONLY, AND 
    % PUT THEM IN THE SAME FOOTING.
% More strings to use when naming outputs:
% (Change below for new model version)
nameStr = {'HT1fFit03Aug','15','SO'};
codeName = [nameStr{1} nameStr{2} nameStr{3}]; % name of this script or function,
    % will be used for key outputs.
% (Change below for new model version) :
par2fit = [2,3,4,5,6,7,8,11,12]; % eta, wEx, sig, wOS11, wOS12, wOS21, wOS22, E0s, E0o.
                               % wEx is IN (it is the self-other channel fusion term, the
                               % narcissism vs borderline discriminator) and EvBlockLR is
                               % OUT (median SE 8.67, 1 of 30 significant). See 
                               % below for description.
par2fitN  = length(par2fit);  % To be used for BIC
groupParN = 0;                % params derived from whole group
codeTesting = 0;  % if set to non-zero, use various debugging settings - see below.
toDo = 1:30;      % which participant from the 3 Aug 26 dataset to do
selfMod = [1,15]; % [1 is the t/DM model, 15] is the evaluations model fitted here.

%% initial directory work - read in or set by hand key directories to work with
cwd = cd;         % just a record of where all this is being run from.
fs = filesep();   % the character that separates folders from subfolders in the filesystem.
% directories to use, and where to to output results - see dirs.outPath below.   
% Paths. This used to call where2findHT1, which hard-codes Michael's Dropbox tree and
% then calls cd on it "to check it's OK", with an empty catch, so on any other machine
% it errored before anything else could happen. Everything here is derived from the
% folder this script sits in instead, and nothing is hard-coded.
thisScript = which(mfilename);
if isempty(thisScript); datdir = [cwd fs]; else; datdir = [fileparts(thisScript) fs]; end
addpath(genpath(datdir));
dirs.HT1St2MentSOres = datdir;   dirs.sandpit = datdir;   dirs.code = datdir;

% for testing, direct outputs to a rough work, 'sandpit' directory:
if codeTesting >= 1
    toDo = 1:2;
    warning('codeTesting is >=1, so few fits done and directed to sandbox directory.');
    outDir = dirs.sandpit; 
    addpath(outDir);
elseif codeTesting == 0
    warning('NO codeTesting, so outputs directed to dirs.HT1St2MentSOres.');
    outDir = dirs.HT1St2MentSOres; 
    addpath(outDir);
else
    error('codeTesting code encountered is not provided for');
end
% HT1tFit03Aug01 saves these as P, D and fit, while the loop below wants p, d and
% tFit, so accept either rather than failing on the names:
tFitFile = [datdir 'HT1tFit1to30a.mat'];
if ~exist(tFitFile,'file')
    error(['%s does not exist. It is an OUTPUT of the decision-making fit, not ' ...
           'something that ships with the code: it holds the pass-one parameters ' ...
           'this script conditions on. Run HT1tFit03Aug01 first.'], tFitFile);
end
S = load(tFitFile);
if     isfield(S,'d');     d = S.d;      elseif isfield(S,'dSave'); d = S.dSave;
elseif isfield(S,'D');     d = S.D;
else;  error('HT1tFit1to30a.mat holds no d, dSave or D'); end
if     isfield(S,'p');     p = S.p;      elseif isfield(S,'pSave'); p = S.pSave;
elseif isfield(S,'P');     p = S.P;
else;  error('HT1tFit1to30a.mat holds no p, pSave or P'); end
if     isfield(S,'tFit'); tFit = S.tFit;  elseif isfield(S,'fit'); tFit = S.fit;
else;  error('HT1tFit1to30a.mat holds neither tFit nor fit'); end
clear S;
    
%% Form key headers for csvs, tables etc.
% Columns headings for the output. 
% Column headings for output. First, retain key ones from decision-making fits.
%  Maximum extend of evaluation / feeling param vector for RLish models. '*' = current
%                 *            *     *      *     *     *                     *       *      *
%          1      2     3      4     5      6     7     8       9      10     11     12     13 
%       lambda   eta   wEx    sig   wOS11 wOS12 wOS21  wOS22  lambda2  lps    E0s    E0o  EvBlockLR

hdEval = {'lambda','eta','wEx','sig','wOS11','wOS12','wOS21','wOS22','lambda2','lps','E0s','E0o','EvBlockLR'};
hdEvalFit = {'slpf','sllSf','sllOf','BICf','sllt','slpt','BICt',};
hdDM   = {'Spref1', 'Spref2','Spref3','Spref4','prevp','prevu',...
          'prevPri1','prevPri2','prevPri3','prevPri4',...
          'SPartp','SPartu','SPartnPr1','SPartnPr2','SPartnPr3','SPartnPr4',...
          'T','blockLR'};
errHd = {'numID','lambdaSE','etaSE','wExSE','sigSE','wOS11SE',...
         'wOS12SE','wOS21SE','wOS22SE','lambda2SE','lpsSE','E0sSE','E0oSE','EvBlockLRSE','fRCondHess'...
         'prevpSE','prevuSE','SPartpSE','SPartuSE','TSE','blockLRSE','tRCondHess','prolificID' };
fErrHd = errHd(2:14); fErrHd = fErrHd(par2fit); % for easy reference of fitted values

hdf = ['numID',hdEval,hdEvalFit,hdDM];
hd2=hdf;      % One of these will act as a backup copy.
totColN = length(hd2);     
fitMeasN  = length(hdEvalFit);  % this many columns for slpf, sllf, BIC ...
feelMeasN = length(hdEval);     % max num of evaluation model params ...      

%% Select the initial values for evalRL params from the vector below. SEE LATER FOR FIXED/DEFAULTS.
%  (Change below for new model version). * = current . 
%                         *              *       *      *     *     *                  *     *      *
% %               1       2       3      4       5      6     7     8       9     10   11   12      13    
% %            lambda    eta     wEx    sig    wOS11 wOS12 wOS21  wOS22  lambda2  lps  E0s  E0o  EvBlockLR
feelpInit = [    -20,    1.9,  -1.5,   -0.56,  0.26, 2.41,  2.41, 0.26,   -20,   -20   0     0     0    ];
% wEx now starts at -1.5, i.e. invlogit ~ 0.18, rather than at -20 which is the
% boundary and a poor place to start a parameter you intend to estimate.   
feelpInit = feelpInit(par2fit);

%% Extract decision making, i.e. _t_ reat-each-other data from fit24f _t_                 
totPtN = length(tFit);
grandTrN = 0; for ptN=1:totPtN grandTrN=grandTrN+sum(p{ptN}.trN); end
warning(['We assume that if there are group-level params, they''re based on ' ...
         num2str(grandTrN) ' trials, grand total. IS THIS RIGHT?']);

% for combined key decision-making and eval. fitted params. For output,
% may also use array fitHT1tBf which will not have the string column of prolific ID, 
% only the numerical PID. Preparet the storage variables:
tabHd = [string(hdf) "prolificID"];
fitHT1TBF = array2table(nan(totPtN,totColN),VariableNames=tabHd(1:end-1));
fitHT1TBF.prolificID = repmat("",totPtN,1);
% For estimates of errors:
errMeasN = length(errHd) - 2;  % It has 2 ID columns too!
errHd = string(errHd);
errHT1TBF = array2table(nan(totPtN,errMeasN+1),VariableNames=errHd(1:end-1));
errHT1TBF.prolificID = repmat("",totPtN,1);
% Now fill in what we already have:
for ptN=1:totPtN
    % fill in two forms of ID:
    fitHT1TBF{ptN,"numID"} = str2num(d{ptN}.PID);
    fitHT1TBF{ptN,"prolificID"} = string(d{ptN}.code);
    errHT1TBF{ptN,"numID"} = str2num(d{ptN}.PID);
    errHT1TBF{ptN,"prolificID"} = string(d{ptN}.code);


    % REM In HT1, the SPartnPr1,... and Spref1,... are not filled in.
    % They were not fitted in previous versions either (were they self-reported?)

    % fill in fitted parameters:
    % Fitted by grid search only:
    fitHT1TBF{ptN,22:25} = d{ptN}.Spref;
    % Fitted by grid search + fmincon :
    fitHT1TBF{ptN,'prevp'} = d{ptN}.prevp;
    fitHT1TBF{ptN,'prevu'} = d{ptN}.prevu;
    fitHT1TBF{ptN,'SPartp'} = d{ptN}.Spartp;
    fitHT1TBF{ptN,'SPartu'} = d{ptN}.Spartu;
    fitHT1TBF{ptN,'T'} = d{ptN}.T;
    fitHT1TBF{ptN,'blockLR'} = d{ptN}.blockLR;
    
    % REM tHess variables are:
    % [ps.prevp, ps.prevu, ps.SPartp, ps.SPartu, ps.T, ps.blockLR]
    se = nan(size(d{ptN}.tHess,1),1);
    % Before attempting to estimate Hessian based errors, see that it is well conditioned:
    rcondH = rcond(d{ptN}.tHess);
    if rcondH > 1e-8 && isfinite(abs(rcondH))
        covMat = d{ptN}.tHess \ eye(size(d{ptN}.tHess));
        se = sqrt(diag(covMat));
        try
            se(d{ptN}.tFixed) = NaN;  % prevu and SPartu were held fixed, so their entries
        end                           % in the embedded Hessian are placeholders, not SEs.
    else
        warning('Inverting the Hessian failed, so leaving SEs as NaN');
    end
    errHT1TBF{ptN,string({'prevpSE','prevuSE','SPartpSE','SPartuSE','TSE','blockLRSE'})} = se'; 
    errHT1TBF{ptN,'tRCondHess'} = rcondH;

end

mslf2 = {};

% REM example of d{}:
    %     code: '62e3d33894be1341b31b0e98'
    %       PID: '0803185301909'
    %      date: [2026 8 6 14 34 32.3393]
    %    evo_hd: {1×42 cell}
    %       evo: {[33×42 double]  [33×42 double]  [30×42 double]}
    %     feelm: NaN
    %    evalRL: [NaN NaN NaN NaN NaN NaN NaN NaN NaN]
    % HT1DatCol: [1 2 3 4 7 8 21 9 10 11 12 26 29]
    %     Spref: [0 0 3 3]
    %   prevPri: [0.2242 0.2057 0.2352 0.3349]
    %   PartnPr: [8.3092e-05 0.0039 0.0860 0.9100]
    %         T: 0.2647
    %   blockLR: 0.2203
    %     feelu: NaN
    % prolDemog: [0×23 table]
    %     prevp: 0.3389
    %     prevu: -4.9993
    %    Spartp: 0.9900
    %     tHess: [6×6 double]
    %    Spartu: 1.4811
   
%%  Construct weakly informative priors, unless priors supplied ~~~~~~
try 
  psPr.evalRL0;  
catch
  [psPr.evalRL0, ~] = priParEvalRLHT1;  % 13 columns. priParEvalRL returns 12, which
end                                     % leaves EvBlockLR with no prior at all.
if size(psPr.evalRL0,2) < 13
    error('psPr.evalRL0 has %d columns, need 13. Use priParEvalRLHT1.', size(psPr.evalRL0,2));
end


%% main feelings fit loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
timeStamp = string(datetime('now', 'Format', 'yyyyMMddHHmmss'));
fitName = [codeName num2str(toDo(1)) 'to' num2str(toDo(end)) char(timeStamp)];   % for storage etc.
cd(outDir); 

for ptN= toDo

      % record the version of the model(s) used for t/decisionmaking and 
      % f/evaluations in both p and d:
      p{ptN}.selfMod = selfMod;       d{ptN}.selfMod = selfMod;
      P = p{ptN};                     D = d{ptN}; 
      errHT1TBF{ptN,'prolificID'} = string(D.code);       % prolific ID
      errHT1TBF{ptN,'numID'} = str2num(D.PID);     % numerical ID

      try              % ...  which we will now check for and if necessary fill in:   
            P.rowPolComb;
      catch
            % was: P.rowPolComb = basicP.rowPolComb; etc., but basicP is never defined
            % anywhere in this script, so if the try ever failed the catch threw its own
            % undefined-variable error on top. Rebuild the deterministic task quantities.
            warning(['P.rowPolComb missing for pt ' num2str(ptN) ', rebuilding via prepParHT1']);
            P = prepParHT1(P);
      end
      % Change below for new model version :
      P.selfMod = selfMod;  % [1,15]. ,1] would be original decision model = 1,
                           %  RLish eval model is 3, 7, 10 etc.
      P.synth = [0,0];     % Explicitly say that we don't simulate feelings data.
      
      % Declare the function to be minimized by fmincon. The follwing has to be
      % re-declared every time we want to do the fit, not e.g. before this loop :
      pS    = tFit{ptN}.ps;
      pS.ID = D.PID;          % This should be the corresp. numerical ID.
      details=0; 
      % The following line contains all the important defaults. The to-be-fitted
      % will ofc be replaced within the likelihood fn, e.g. HT1lp02f
      % (Change below for new model version) :
      % Values for every entry NOT in par2fit. Entries that ARE in par2fit get
      % overwritten inside HT1lp02f, so their value here does not matter.
      % Entry 9 is lambda2: -20 gives invlogit ~ 0, evaluations track RETURNS;
      %                     +20 gives invlogit ~ 1, evaluations track PREDICTION ERRORS.
      % 13b as originally run had +20, i.e. a PE fit, despite model 13 being documented
      % as the returns one. Pick one on purpose and say so in nameStr.
      lambda2Transf = -20;   % <<< THE returns / PE switch. -20 returns, +20 PE.
      % TRANSF lambda  eta wEx sig  wOS11  wOS12 wOS21  wOS22 lambda2 lps E0s E0o EvBlockLR 
      pS.evalRL  =[-20,1.9,-20,-0.56, 0.26, 2.41, 2.41,  0.26, lambda2Transf, -20, 0, 0, -1.4];  
      % Sanity check that the free and fixed index sets partition 1:13:
      parFixed = setdiff(1:13, par2fit);
      if ~isequal(sort([par2fit parFixed]), 1:13); error('par2fit does not partition 1:13'); end
      % Health warning on the returns / PE contrast. lambda, entry 1, is fixed at
      % invlogit(-20) ~ 2e-9, so the return expectation never moves off its initial
      % value R0. The prediction error is then Ret - R0, a constant offset from Ret, and
      % that offset is absorbed exactly by the free intercept E0. Verified numerically:
      % lambda2 = 1 and lambda2 = 0 with E0 shifted by -Wr*R0 give the same log
      % likelihood to eight decimal places. So with lambda fixed, 'returns based' and
      % 'PE based' are the SAME model with a reparameterised intercept, and comparing
      % them is not a model comparison. To make the contrast real, free lambda:
      %     par2fit = [1,2,4,5,6,7,8,11,12,13];
      if ~ismember(1,par2fit) && ismember(11,par2fit) && ismember(12,par2fit)
          warning(['lambda is fixed while E0s and E0o are free, so lambda2 is not ' ...
                   'identifiable and the returns/PE label is cosmetic. See comment above.']);
      end
      
      mLP = @(feelp)HT1lp02f( feelp, pS, D, P, psPr, approvs2fit, details);
      % boring: have to specify empty 'linear constraints' in order to get to 
      % the arguments for the lower and upper bounds, acc. to the doc fmincon example ...
      A = []; b = []; Aeq = [];  beq = [];
      % NB fit will take place in transformed space, hence the values of lower and upper bounds.
      % first row has upper bounds for feelp, second has lower. Maximally,
      %prHd={'lambda','eta','wEx','sig','wOS11','wOS12','wOS21','wOS22','lambda2','lps','E0s','E0o','EvBlockLR'};
      Bounds =[[ 10,    10,    10,   5,   20,     20,    20,     20,      10,     10,   20,   20,     10] ; ... 
               [-10,   -10,   -10,  -5,  -20,    -20,   -20,    -20,     -10,    -10,  -20,  -20,    -10]];   
      % (two commas were missing on those rows: inside brackets MATLAB treats the gap
      %  as a separator anyway, so it parsed, but it is one typo away from not doing.)
      ub = Bounds(1,par2fit);  % Only these variables are optimized in version 03Aug13bBoth
      lb = Bounds(2,par2fit);  % Only these variables are optimized in version 03Aug13bBoth
      
      try
          %% Flatness check. A free parameter that does not move the objective is not
          % merely unidentifiable: it puts an exactly zero row and column in the
          % Hessian, so rcond is 0 and EVERY standard error comes back NaN. This is
          % what evalRL(13) was doing. Cheap to check, so check every time.
          fdH = 1e-4;  fdG = nan(1,par2fitN);
          for kPar = 1:par2fitN
              tp = feelpInit; tp(kPar) = tp(kPar) + fdH;
              tm = feelpInit; tm(kPar) = tm(kPar) - fdH;
              fdG(kPar) = (mLP(tp) - mLP(tm)) / (2*fdH);
          end
          deadPar = hdEval(par2fit(abs(fdG) < 1e-9));
          disp('  d(objective)/d(param) at the start point:');
          for kPar = 1:par2fitN
              disp(sprintf('    %-10s %12.4e', hdEval{par2fit(kPar)}, fdG(kPar))); %#ok<*DSPS>
          end
          if ~isempty(deadPar)
              error(['These free parameters leave the objective unchanged: ' ...
                      strjoin(deadPar,', ') '. Fix the model before fitting.']);
          end

          %% Crucial optimization call ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          % (Change below for new model version)
          disp(['Now fitting pt ',num2str(ptN) ' using fmincon, ' codeName ', for ' strjoin(hdEval(par2fit),',')]);
          [feelpOpt, mmLL,~, output, ~, ~, hessian] = fmincon(mLP, feelpInit, A, b, Aeq, beq,lb,ub );
          %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
          % now run HT1lpf once again to return details and augment with the 
          % fit quality indices :
          mslf2{ptN} = HT1lp02f( feelpOpt, pS, D, P, psPr, [1,1], 1); 
          mslf2{ptN}.fitOut = output;
          mslf2{ptN}.Hess = hessian;
          fRCondHess = rcond(hessian);  % Is Hessian well-conditioned ( > 1e-8)
          errHT1TBF{ptN,'fRCondHess'} = fRCondHess;
          if fRCondHess > 1e-8 && isfinite(abs(fRCondHess))
            covMat = hessian \ eye(size(hessian));
            se = sqrt(diag(covMat))';
            errHT1TBF{ptN,string(fErrHd)} = se;         
          else
            warning('Inverting the Evaluations Fit Hessian failed, so leaving SEs as NaN');
          end

          feelP = pS.evalRL; 
          feelP(par2fit) = feelpOpt;   
          % (Change below for new model version) :
          % NO DUPLICATED ONE IN 10 feelP(7) = feelpOpt(end);  % this is the duplicated one.
          
          % Calculate BIC, with correction for group-wide params.
          %   Note this assumes all pts had same number of trials, which
          %   isn't strictly true, but is approximately true (only first two had slightly different):
          % n is the number of RATINGS entering sllf, which is two per trial, the
          % satisfaction one and the trust one, not one per trial:
          nRatings = 2*sum(P.trN);
          BICnN = -2*sum(mslf2{ptN}.sllf(1:2)) + ...       % deliberately not sllt here ! 
                  + par2fitN * (log(nRatings)) + ...   % was + par2fitN*(log(sum(P.trN))) + ...
                  + (groupParN / totPtN)  * log (grandTrN);

          fitHT1TBF{ptN,string(hdEval)} = feelP;
          % REM hdEvalFit = {'slpf','sllSf','sllOf','BICf','sllt','slpt','BICt',};

          % mpostf is MINUS the sum log posterior while the column is called slpf, so
          % negate it and the column means what its name says. Sign differs from the
          % csv files produced before Aug 26.
          fitHT1TBF{ptN,string(hdEvalFit(1:5))} = [-mslf2{ptN}.mpostf, ...
                 mslf2{ptN}.sllf(1:2), BICnN, mslf2{ptN}.sllt]; 
      catch
          warning(['fmincon MAP fitting failed for pt ' num2str(ptN)]);
          % if no valid data etc., create appropriate empty outputs:
          % Needed in case of invalid data:
          feelP = nan(1,feelMeasN-4) ; % should eq. to the length of feelP
          mslf2{ptN}.mpostf = nan; mslf2{ptN}.sllf= nan(1,2); BICnN = nan;
          mslf2{ptN}.fitOut = [];  mslf2{ptN}.Hess= nan(length(feelpInit));

      end
  
      save([outDir fitName '.mat'], 'mslf2','fitHT1TBF','errHT1TBF','d', 'p'); 
      fitHT1tBf = fitHT1TBF{:,1:end-1};
      mat2csv2Dfl(fitHT1tBf, [outDir fitName '.csv'], 0,1, hdf ); 
      errHT1tBf = errHT1TBF{:,1:end-1}; 
      mat2csv2Dfl(errHT1tBf, [outDir 'err' fitName '.csv'], 0,1, errHd(1:end-1) ); 

end

disp([fitName ' done, with NO group params.']);
nBadHess = sum(~isfinite(errHT1TBF.fRCondHess) | errHT1TBF.fRCondHess <= 1e-8);
disp(['Hessians too ill-conditioned to invert: ' num2str(nBadHess) ' of ' num2str(totPtN) ...
      '. Their standard errors are NaN and should not be used.']);
cd(cwd)

