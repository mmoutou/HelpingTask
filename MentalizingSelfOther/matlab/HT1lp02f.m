function mslf = HT1lp02f( restp, pS, d, p, psPr, approvs2use, details)
%MSLF= HT1LP2F( restp, pS, d, p, psPr, approvs2use, details) wrapper giving minus sum ln post. of emot. eval. params 
%   from 'Helping Task' (HT1) expt., once all the interaction  
%   treat-each-other parameters have been fitted and incorporated in pS.  
%   Text-search (Change below for new model version)
%       Returns mainly the minus-sum-log-posterior (prior being psPr over 
%   all params BUT ONLY LOGLIK+LOGPRI CORRESP. TO EMOTIONAL EVALUATIONS PART RETURNED HERE.)
%     - 'if details', also the sum log likelihood over each other's treatment and
%       feelings ratings separately, and also a copy of the parameters as used, i.e.
%       having replaced the relevant parts of 'baseline' params pS with restp.
%     - restp is the emotional evaluation params e.g. evalRL or [feelm, feelu] to use, which
%       should be appropriate to use this fn with fmincon and similar.
%     - pS contains all the rest of the necessary (treat-each-other & 
%       'default' feeling / eval) params
%     - d contains the data. May be synthetic!
%     - p contains the *experiment*, not the participant parameters. If it has 
%       a field p.synth which isn't [0,0], AND 'details' is nonzero, 
%       then mslf has a field mslf.dsynth with the synthetic data.
%     - p also contains info about the version of the (likelihood) model to use.
%       p.selfMod . 
%     - psPr has the priors over [ps.feelm, ps.feelu]
%     - approvs2use allows to use the data over self, other or 
%       both [1,0], [0,1], [1,1] respectively.

% Set below to 1 to test the code, but see below return statement at
% end of function for how to 
codeTest = 1;

try
   runType.selfMod = p.selfMod;
catch
   runType.selfMod = [1,1];  % This is the original, ',1]' 
           % means with joint preference evaluation / feeling model.
end
try
  approvs2use;
catch
  approvs2use = [1,1];  % weights to place on the log lik for self and other.
  % in summing mslf. if [1,0,0] return only for self etc. (no other->self in HT1)
end
try
  details;
catch
  details = 0; 
end
if details % check for providing synthetic data only if details is true.
    % Otherwise set to [0,0] any previously specified p.synth, so 
    % that we can help signal (with 'details') that synth data isn't needed.
    try 
        p.synth(2);
    catch
        p.synth = [0,0];
    end
end
        
%% consturct new participant-parameter structure, keeping the 
%  crucial ps.Spref well alone :
ps = pS;
try                   % as these may already be formed ...
  ps.prevPri;
  ps.SPartnPr;
catch
  ps.prevPri  = noisyBino(ps.prevp,ps.prevu,p.Nl);
  ps.SPartnPr = noisyBino(ps.SPartp,ps.SPartu,p.Nl);
end
% (Change below for new model version)
% ( here, add elseif for new model version ) :
if runType.selfMod(2) == 1   % set up priors for 1st joint preferences eval. model.
                             % and use restp to set up ps.[appropriate eval. params]
    try 
      psPr.feelm0;
    catch
      clear psPr;
      psPr.feelm0 = nan; psPr.feelu0 = nan;
    end
    if isnan(psPr.feelm0)
      psPr.feelm0 = [1.025, 10.0];    % A and B for gampdf for m of rat2resp - flat max near 1
      psPr.feelu0 = [1.025, 10.0];    % A and B for gampdf for u of rat2resp - flat max near 1
    end
    ps.feelm = restp(1); 
    ps.feelu = restp(2);
    try  rmfield(ps, 'evalRL'); end

elseif runType.selfMod(2) == 2   % set up priors for 1st RL autoregression model.
           % and use restp to set up ps.[appropriate eval. params] - only 5 free params.
    try 
      psPr.evalRL0;
    catch
      clear psPr;
      psPr.evalRL0 = nan;
    end
    if isnan(psPr.evalRL0)  % nan can also be given as an argument to force default.
        [psPr.evalRL0, ~] = priParEvalRL();
         % will cater later for those set deterministically, so prior=1.
    end
% D.RLEval(v.trial+1,2:3) =  ( (1-eta)*[[1-wEx, wEx]; [wEx, 1-wEx]] * D.RLEval(v.trial,2:3)' + ...
%                                 eta * wOS*((1-lambda2)*Ret + lambda2*PE))' ; 
% %               1      ==>2         3         ==>4     ==>5  ==>[6 ==== 7]     8       9     10     11    12   ==>13
% %            lambda      eta       wEx        sig      wOS11 wOS12  wOS21  wOS22  lambda2    lps    E0s   E0o   EvBlockLR
    ps.evalRL([2,4,5,6,13]) = restp; % Unique params, incl first of yoked pairs.
    ps.evalRL(7) = ps.evalRL(6) ;    % Yoked params.

    try  rmfield(ps,{'feelm','feelu'}); end

elseif  runType.selfMod(2) == 3   % set up priors for 2nd RL autoregression model, PE
                        % based, and use restp to set up ps.[appropriate eval. params]
    try 
      psPr.evalRL0;
    catch
      clear psPr;
      psPr.evalRL0 = nan;
    end
    if isnan(psPr.evalRL0)  % nan can also be given as an argument to force default.
        [psPr.evalRL0, ~] = priParEvalRL();
         % will cater later for those set deterministically, so prior=1.
    end
% D.RLEval(v.trial+1,2:3) =  ( (1-eta)*[[1-wEx, wEx]; [wEx, 1-wEx]] * D.RLEval(v.trial,2:3)' + ...
%                                 eta * wOS*((1-lambda2)*Ret + lambda2*PE))' ; 
% %            ==>1      ==>2         3        ==>4   ==>5   ==>[6 ==== 7]     8       9      10     11     12  ==>13
% %            lambda      eta       wEx       sig    wOS11  wOS12   wOS21  wOS22   lambda2   lps    E0s    E0o   EvBlockLR
  ps.evalRL([1,2,4,5,6,13]) = restp; % Unique params, incl first of yoked pairs.
  ps.evalRL(7) = ps.evalRL(6) ;      % Yoked params.

  try  rmfield(ps,{'feelm','feelu'}); end

elseif  runType.selfMod(2) == 4   % set up priors for 2nd RL autoregression model, PE
                        % based, and use restp to set up ps.[appropriate eval. params]
                        % Has yoked wOS21, and I think it has no wOS22
    try 
      psPr.evalRL0;
    catch
      clear psPr;
      psPr.evalRL0 = nan;
    end
    if isnan(psPr.evalRL0)  % nan can also be given as an argument to force default.
        [psPr.evalRL0, ~] = priParEvalRL();
         % will cater later for those set deterministically, so prior=1.
    end
% D.RLEval(v.trial+1,2:3) =  ( (1-eta)*[[1-wEx, wEx]; [wEx, 1-wEx]] * D.RLEval(v.trial,2:3)' + ...
%                                 eta * wOS*((1-lambda2)*Ret + lambda2*PE))' ; 
% %              1       ==>2      ==>3        ==>4   ==>5   ==>[6 ==== 7]     8       9       10         11     12  ==>13
% %            lambda      eta       wEx         sig    wOS11  wOS12  wOS21  wOS22   lambda2   lps        E0s    E0o   EvBlockLR
  ps.evalRL([2,3,4,5,6,13]) = restp; % Unique params, incl first of yoked pairs.
  ps.evalRL(7) = ps.evalRL(6) ;      % Yoked params.
 
  try  rmfield(ps,{'feelm','feelu'}); end

elseif  runType.selfMod(2) == 5   % set up priors for first RL autoregression model return AND PE
            % based, with yoked wOS12  wOS21 , and use restp to set up ps.[appropriate eval. params]
    try 
      psPr.evalRL0;
    catch
      clear psPr;
      psPr.evalRL0 = nan;
    end
    if isnan(psPr.evalRL0)  % nan can also be given as an argument to force default.
        [psPr.evalRL0, ~] = priParEvalRL();
         % will cater later for those set deterministically, so prior=1.
    end
% D.RLEval(v.trial+1,2:3) =  ( (1-eta)*[[1-wEx, wEx]; [wEx, 1-wEx]] * D.RLEval(v.trial,2:3)' + ...
%                                 eta * wOS*((1-lambda2)*Ret + lambda2*PE))' ; 
% %           ==>1       ==>2         3        ==>4   ==>5   ==>[6 ==== 7]     8    ==>9       10         11     12  ==>13
% %            lambda      eta       wEx         sig    wOS11  wOS12  wOS21  wOS22   lambda2    lps       E0s    E0o   EvBlockLR
    ps.evalRL([1,2,4,5,6,9,13]) = restp;    % Unique params, incl. first of yoked pairs.
    ps.evalRL(7) = ps.evalRL(6);            % Yoked params
    try  rmfield(ps,{'feelm','feelu'}); end
    
elseif  runType.selfMod(2) == 6   % set up priors for 2nd RL autoregression model, PE 
              % based, and use restp to set up ps.[appropriate eval. params]
              % I think it has no wOS22 yet.
    try 
      psPr.evalRL0;
    catch
      clear psPr;
      psPr.evalRL0 = nan;
    end
    if isnan(psPr.evalRL0)  % nan can also be given as an argument to force default.
        [psPr.evalRL0, ~] = priParEvalRL();
         % will cater later for those set deterministically, so prior=1.
    end
% D.RLEval(v.trial+1,2:3) =  ( (1-eta)*[[1-wEx, wEx]; [wEx, 1-wEx]] * D.RLEval(v.trial,2:3)' + ...
%                                 eta * wOS*((1-lambda2)*Ret + lambda2*PE))' ; 
% %              1       ==>2                  ==>4   ==>5    ==>6   ==>7      8        9       10      11     12  ==>13
% %            lambda      eta       wEx         sig    wOS11  wOS12  wOS21  wOS22   lambda2    lps     E0s    E0o   EvBlockLR
  ps.evalRL([2,4,5,6,7,13]) = restp;  % No yoked params in version 6. 
    try  rmfield(ps,{'feelm','feelu'}); end

elseif  runType.selfMod(2) == 7   % set up priors for RLEval07, the crucial 'full wOS' model.
    % This has params eta, sig, and wOS11 to 22 (full wOS). 
    % See below in this block for how this is uses ps.evalRL (the params) and further down for
    % ps.evalRL0 (the priors on these).
    try 
      psPr.evalRL0;
    catch
      clear psPr;
      psPr.evalRL0 = nan;
    end
    if isnan(psPr.evalRL0)  % nan can also be given as an argument to force default.
        [psPr.evalRL0, ~] = priParEvalRL();
         % will cater later for those being set deterministically, so prior=1.
    end
% D.RLEval(v.trial+1,2:3) =  ( (1-eta)*[[1-wEx, wEx]; [wEx, 1-wEx]] * D.RLEval(v.trial,2:3)' + ...
%                                 eta * wOS*((1-lambda2)*Ret + lambda2*PE))' ; 
% %              1       ==>2                  ==>4   ==>5    ==>6     ==>7   ==>8      9       10      11     12  ==>13
% %            lambda      eta       wEx         sig    wOS11  wOS12  wOS21  wOS22   lambda2    lps     E0s    E0o   EvBlockLR
  ps.evalRL([2,4,5,6,7,8,13]) = restp;  % No yoked param in this version.
    try  rmfield(ps,{'feelm','feelu'}); end

elseif  runType.selfMod(2) == 8   % set up priors for model 8 RL autoregression model, PE
                        % based, and use restp to set up ps.[appropriate eval. params]
    try 
      psPr.evalRL0;
    catch
      clear psPr;
      psPr.evalRL0 = nan;
    end
    if isnan(psPr.evalRL0)  % nan can also be given as an argument to force default.
        [psPr.evalRL0, ~] = priParEvalRL();
         % will cater later for those set deterministically, so prior=1.
    end
% D.RLEval(v.trial+1,2:3) =  ( (1-eta)*[[1-wEx, wEx]; [wEx, 1-wEx]] * D.RLEval(v.trial,2:3)' + ...
%                                 eta * wOS*((1-lambda2)*Ret + lambda2*PE))' ; 
% %              1       ==>2                  ==>4   ==>5    ==>6   ==>7      8       9     ==>10      11     12  ==>13
% %            lambda      eta       wEx         sig    wOS11  wOS12  wOS21  wOS22   lambda2    lps     E0s    E0o   EvBlockLR
  ps.evalRL([2,4,5,6,7,10,13]) = restp;  % No yoked param in this version. 
    try  rmfield(ps,{'feelm','feelu'}); end

elseif  runType.selfMod(2) == 9   % set up priors for model 9 RL autoregression model, PE
         % based, and use restp to set up ps.[appropriate eval. params]
    try 
      psPr.evalRL0;
    catch
      clear psPr;
      psPr.evalRL0 = nan;
    end
    if isnan(psPr.evalRL0)  % nan can also be given as an argument to force default.
        [psPr.evalRL0, ~] = priParEvalRL();
         % will cater later for those set deterministically, so prior=1.
    end
% D.RLEval(v.trial+1,2:3) =  ( (1-eta)*[[1-wEx, wEx]; [wEx, 1-wEx]] * D.RLEval(v.trial,2:3)' + ...
%                                 eta * wOS*((1-lambda2)*Ret + lambda2*PE))' ; 
% %              1       ==>2                  ==>4   ==>5    ==>6     ==>7   ==>8       9    ==>10     11     12  ==>13
% %            lambda      eta       wEx         sig    wOS11  wOS12  wOS21  wOS22   lambda2    lps     E0s    E0o   EvBlockLR
  ps.evalRL([2,4,5,6,7,8,10,13]) = restp;  % No repeated param in this version. 
    try  rmfield(ps,{'feelm','feelu'}); end

elseif  runType.selfMod(2) == 10  % set up priors for RLEval10 model, which is PE based, but
    % otherwise like full wOS (ie. RLEval07 ), and introduces E0 = [E0s E0o]'. 
    % See below in this block for how this is uses ps.evalRL (the params) and further down for
    % ps.evalRL0 (the priors on these).
    try 
      psPr.evalRL0;
    catch
      clear psPr;
      psPr.evalRL0 = nan;
    end
    if isnan(psPr.evalRL0)  % nan can also be given as an argument to force default.
        [psPr.evalRL0, ~] = priParEvalRL();
         % will cater later for those set deterministically, so prior=1.
    end
% D.RLEval(v.trial+1,2:3) =  [E0s E0o]' + ...
%                            ( (1-eta)*[[1-wEx, wEx]; [wEx, 1-wEx]] * D.RLEval(v.trial,2:3)' + ...
%                                 eta * wOS*((1-lambda2)*Ret + lambda2*PE))' ; 
% %  ==>1   ==>2         ==>4   ==>5    ==>6    ==>7  ==>8   9 fix->max    10  ==>11  ==>12    ==>13
% % lambda   eta    wEx   sig    wOS11  wOS12  wOS21  wOS22   lambda2     lps    E0s    E0o      EvBlockLR
  ps.evalRL([1,2,4,5,6,7,8,11,12,13]) = restp;  % No yoked params. 
    try  rmfield(ps,{'feelm','feelu'}); end

elseif  runType.selfMod(2) == 11  % set up priors for RLEval11 model, which is PE based but
    % otherwise like full wOS (ie. RLEval07 ), and has yoked E0 = [E0s E0o=E0s]'. 
    % See below in this block for how this is uses ps.evalRL (the params) and further down for
    % ps.evalRL0 (the priors on these).
    try 
      psPr.evalRL0;
    catch
      clear psPr;
      psPr.evalRL0 = nan;
    end
    if isnan(psPr.evalRL0)  % nan can also be given as an argument to force default.
        [psPr.evalRL0, ~] = priParEvalRL();
         % will cater later for those set deterministically, so prior=1.
    end
% D.RLEval(v.trial+1,2:3) =  [E0s E0o]' + ...
%                            ( (1-eta)*[[1-wEx, wEx]; [wEx, 1-wEx]] * D.RLEval(v.trial,2:3)' + ...
%                                 eta * wOS*((1-lambda2)*Ret + lambda2*PE))' ; 
% %  ==>1   ==>2         ==>4   ==>5    ==>6    ==>7  ==>8   9 fix->max    10   ==>11 = 12   ==>13
% % lambda   eta    wEx   sig    wOS11  wOS12  wOS21  wOS22   lambda2     lps    E0s     E0o   EvBlockLR
  ps.evalRL([1,2,4,5,6,7,8,11,13]) = restp;  % version 11 has yoked E0o param 12, not set here but below:
  ps.evalRL(12) = ps.evalRL(11);          % yoked param - fix E0o := E0s. 
    try  rmfield(ps,{'feelm','feelu'}); end

elseif  runType.selfMod(2) == 13  % set up priors for RLEval13 model, which is Returns based but
    % otherwise like full wOS (ie. RLEval07 ), and has yoked E0 = [E0s E0o=E0s]'. 
    % See below in this block for how this is uses ps.evalRL (the params) and further down for
    % ps.evalRL0 (the priors on these).
    try 
      psPr.evalRL0;
    catch
      clear psPr;
      psPr.evalRL0 = nan;
    end
    if isnan(psPr.evalRL0)  % nan can also be given as an argument to force default.
        [psPr.evalRL0, ~] = priParEvalRLHT1();   % 13 cols, i.e. incl. EvBlockLR.
         % was priParEvalRL(), which returns only 12 and so leaves parameter 13 with
         % no prior at all. will cater later for those set deterministically, so prior=1.
    end
% D.RLEval(v.trial+1,2:3) =  [E0s E0o]' + ...
%                            ( (1-eta)*[[1-wEx, wEx]; [wEx, 1-wEx]] * D.RLEval(v.trial,2:3)' + ...
%                                 eta * wOS*((1-lambda2)*Ret + lambda2*PE))' ; 
% %    1   ==>2         ==>4   ==>5    ==>6    ==>7  ==>8   9 SET BY CALLER  10  ==>11  ==>12  ==>13
% % lambda   eta    wEx   sig    wOS11  wOS12  wOS21  wOS22   lambda2     lps    E0s    E0o     EvBlockLR
  ps.evalRL([2,4,5,6,7,8,11,12,13]) = restp;  % No yoked params. 
  % NB lambda2, entry 9, is NOT set here, it is inherited from whatever the calling
  % script left in pS.evalRL. That is how model 13 came to be documented as the
  % returns-based one (comment '9 fix->min') while HT1fFit03Aug13bBoth set it to +20,
  % i.e. invlogit(20) = 1, and ran it prediction-error based. Set it deliberately in
  % the calling script and say which you meant.
  try  rmfield(ps,{'feelm','feelu'}); end

elseif  runType.selfMod(2) == 15  % RLEval15, the preregistration model.
    % Derived from 13, with three changes, each following from the 26 Aug 26 full-grid fits:
    %  (a) wEx (entry 3) is FREED. It is the cross-talk between the self and other
    %      evaluation channels, i.e. whether feeling badly about oneself drags the view of
    %      the partner down with it. It was fixed at 0 in every earlier version, so a fused
    %      self-other structure could not be distinguished from an insulated one at all.
    %  (b) EvBlockLR (entry 13) is FIXED, not fitted. Its median standard error was 8.67 on
    %      the transformed scale against 0.13 for sig, and 1 participant of 30 had an
    %      estimate distinguishable from zero. Three partners give two between-partner
    %      transitions, which cannot support a rate. It stays wired into the likelihood
    %      (llfeelHT1b), it is simply no longer free.
    %  (c) lambda2 is set here rather than inherited from the calling script.
    % par2fit = [2,3,4,5,6,7,8,11,12] : nine free against 66 ratings per participant.
    try 
      psPr.evalRL0;
    catch
      clear psPr;
      psPr.evalRL0 = nan;
    end
    if isnan(psPr.evalRL0)
        [psPr.evalRL0, ~] = priParEvalRLHT1();   % 13 cols, incl. EvBlockLR
    end
% D.RLEval(v.trial+1,2:3) =  [E0s E0o]' + ...
%                            ( (1-eta)*[[1-wEx, wEx]; [wEx, 1-wEx]] * D.RLEval(v.trial,2:3)' + ...
%                                 eta * wOS*((1-lambda2)*Ret + lambda2*PE))' ; 
% %    1   ==>2   ==>3  ==>4   ==>5    ==>6   ==>7   ==>8   9 set here  10   ==>11  ==>12  13 fixed
% % lambda   eta    wEx   sig   wOS11  wOS12  wOS21  wOS22   lambda2    lps    E0s    E0o    EvBlockLR
  ps.evalRL([2,3,4,5,6,7,8,11,12]) = restp;  % No yoked params.
  ps.evalRL(9)  = -20;   % invlogit(-20) ~ 0 : evaluations track RETURNS. Stated on purpose.
  ps.evalRL(13) =   0;   % invlogit(0) = 0.5 : half of the previous partner's evaluative
                         % state carries into the next, for everyone. Fixed, not estimated.
  try  rmfield(ps,{'feelm','feelu'}); end

elseif  runType.selfMod(2) == 14  % set up priors for RLEval14 model, which is Returns based but
    % otherwise like full wOS (ie. RLEval07 ), and includes E0 = [E0s E0o=E0s]'. 
    % See below in this block for how this is uses ps.evalRL (the params) and 
    % ps.evalRL0 (the priors on these).
    try 
      psPr.evalRL0;
    catch
      clear psPr;
      psPr.evalRL0 = nan;
    end
    if isnan(psPr.evalRL0)  % nan can also be given as an argument to force default.
        [psPr.evalRL0, ~] = priParEvalRL();
         % will cater later for those set deterministically, so prior=1.
    end
% D.RLEval(v.trial+1,2:3) =  [E0s E0o]' + ...
%                            ( (1-eta)*[[1-wEx, wEx]; [wEx, 1-wEx]] * D.RLEval(v.trial,2:3)' + ...
%                                 eta * wOS*((1-lambda2)*Ret + lambda2*PE))' ; 
% %    1   ==>2         ==>4   ==>5    ==>6    ==>7  ==>8   9 fix->min    10   ==>11 == 12   ==>13
% % lambda   eta    wEx   sig    wOS11  wOS12  wOS21  wOS22   lambda2     lps    E0s    E0o   EvBlockLR
  ps.evalRL([2,4,5,6,7,8,11,13]) = restp;  % Will repeat param E0o := E0s. 
  ps.evalRL(12) = ps.evalRL(11);           % yoked param(s) set here.
  try  rmfield(ps,{'feelm','feelu'}); end
    
    
else % if model not covered ...
%% End of blocks for model alternatives ------------------------    
    
    error('runType.selfMod(2)  not catered for.')
end

%% Only now, and if appropriate, calc log prior of params.
%  This is HT1lp02f, but we deliberately use the HT1 functions for the lnPri and
%  the log-lik below, as the same feelings / evaluations model is used at this level.
if isempty(psPr)
  lnPri = 0;
else
  [~, LPri] = pslPrHT1f(ps, psPr);
  % Correct items deterministically set, hence having trivial prior (certain)
  % (Change below for new model version)
  % ( here add 'elseif' for new model version) :
  if runType.selfMod(2) == 2 
    LPri.evalRL([1,3,7:10,11:12]) = 0;  % 2 has par2fit = [2,4,5,6,13];
  elseif  runType.selfMod(2) == 3 
    LPri.evalRL([3,7:10,11:12]) = 0; 
  elseif  runType.selfMod(2) == 4 
    LPri.evalRL([1,7:10,11:12]) = 0; 
  elseif  runType.selfMod(2) == 5 
    LPri.evalRL([3,7:8,10,11:12]) = 0;  % 5 has par2fit = [1,2,4,5,6,9,13];
  elseif  runType.selfMod(2) == 6  
    LPri.evalRL([1,3,8:10,11:12]) = 0;  % 6 has par2fit = [2,4,5,6,7,13]
  elseif  runType.selfMod(2) == 7 
    LPri.evalRL([1,3,9:10,11:12]) = 0;  % 7 has par2fit = [2,4,5,6,7,8,13];
  elseif  runType.selfMod(2) == 8 
    LPri.evalRL([1,3,8,9,11:12]) = 0;   % 8 has par2fit = [2,4,5,6,7,10,13];
  elseif  runType.selfMod(2) == 9 
    LPri.evalRL([1,3,9,11:12]) = 0;     % 9 has par2fit = [2,4,5:8,10,13]; 
  elseif  runType.selfMod(2) == 10 
    LPri.evalRL([3,9,10]) = 0;          % 10 has par2fit = [1,2,4,5:8,11,12,13]; 
  elseif  runType.selfMod(2) == 11 
    LPri.evalRL([3,9,10,12]) = 0;       % 11 has par2fit = [1,2,4,5:8,11,13]; 
  elseif  runType.selfMod(2) == 13 
    LPri.evalRL([1,3,9,10]) = 0;        % 13 has par2fit = [2,4,5:8,11,12,13]; 
  elseif  runType.selfMod(2) == 14 
    LPri.evalRL([1,3,9,10,12]) = 0;     % 14 has par2fit = [2,4,5:8,11,13]; 
  elseif  runType.selfMod(2) == 15 
    LPri.evalRL([1,9,10,13]) = 0;       % 15 has par2fit = [2,3,4,5:8,11,12]; 
  end
  lnPri = sum(LPri.evalRL);
end

% REM: below, HT1ll1 treat-each-other and feelings parts of likelihood,
%      llt and llf
if (p.synth(1) || p.synth(2)) && details  % to generate synthetic
    % data, 'details' has to be true AND what kind of synthetic 
    % data is needed has to be specified.
    runType.tSynth = p.synth(1);   % contribution (t[reat each other]) part
    runType.fSynth = p.synth(2);   % f[eelings] i.e. Evaluations part
    runType.detailed = 1;
   [llt,llf,Outp] = HT1ll1(ps, d, p, runType,runType.selfMod); 
else
   % ... no synthetic data from HT1ll1
   runType.tSynth = 0; 
   runType.fSynth = 0;
   try
       runType.detailed = details;
   catch
       runType.detailed = 0;
   end
   if runType.detailed == 0
      [llt,llf,~] = HT1ll1(ps, d, p, runType,runType.selfMod);  
      Outp = [];
   else
      [llt,llf,Outp] = HT1ll1(ps, d, p, runType,runType.selfMod);  
   end
end

%% outputs: sum log posterior 
%  and other details if asked for ...
appr2use = approvs2use;
% if there are no Other->Self eval, llf(3) should be 0
% and approvs2use only have 2 elements. So cater for that for dot products below:
if length(appr2use) < 3; appr2use = [appr2use 0]; end
if details
  mslf.sllt  = llt;
  mslf.sllf  = llf;
  mslf.mpostf = -(dot(llf,appr2use) + lnPri);
  mslf.ps = ps;  % this should include the contents of restp in 
                 % their right form, e.g. evalRL, feelm and feel u, etc.
  if ~isempty(Outp)  % this means that we have generated synthetic data
      try mslf.dsynth    =  Outp.dsynth; catch ; end
      try mslf.llreal = Outp.llreal; catch ; end   % ll for the real data, and trial by trial measures:
       % [ll, v.trial, v.harshLev, v.Spol, v.nmoves(1,:,v.harshLev),  v.nmoves(2,:,v.harshLev), goalDivergence ]
      mslf.feelPhd   = Outp.feelPhd;   % header for the below
      mslf.feelPrat  = Outp.feelPrat;  % trial by trial evaluation / feeling data in array form.
      try mslf.feelKeyHd = Outp.feelKeyHd; catch ; end
      try mslf.feelKey   = Outp.feelKey;   catch ; end
      mslf.contrEvalHd  = Outp.contrEvalHd;  % header for work contributions provided.
      mslf.contrEval    = Outp.contrEval;    % Work contrib. data, in the same row format as in d.evo.
  end
else
  mslf = -(dot(llf,appr2use) + lnPri);
end

return;  % end of function HT1lp02f
%% Instructions for running example / demo / debug 
% cwd=cd; cd('C:\Users\mmpsy\Dropbox\FIL_aux\MSc_iBSc_PhD_student\Gosalia_(Meera_iBSc)\SelfOtherEvalModelsProject-shared\SelfEvalModelCoding\dataHT1\resultsDec')
% %               1       =>2         3         =>4     =>5   =>6=====7     8       9       10       11  ==  12     13
% %            lambda      eta       wEx         sig    wOS11 wOS12 wOS21  wOS22  lambda2   lps      E0s     E0o   EvBlockLR
%   evalRL =[logit(1e-9),logit(0.1),logit(1e-9),log(0.5),0.33, 0.67, 0.67,  0,  logit(1e-9),logit(1e-9) ] 
% load('fitHT1t_fmri.mat'); ptN=3; prtn=2;
% d1=D{ptN}.here.d; p1=D{ptN}.here.p;  ps1=ps{ptN,prtn}; ps1.evalRL=evalRL; p1.synth=[0,1];
% detail=1; restp1 = evalRL([2,4,5,6]); approvs2use = [1,1]; p1.selfMod = [1,2]; % for RLish / PEkernel eval model
% mslf = HT1lp02f( restp1, ps1, d1, p1, nan, approvs2use, detail)

