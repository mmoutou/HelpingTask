function [ll, v, D] = llfeelHT1b(ps, P, D, v)
%  [ll, v, D] = llfeelHT1b(ps, P, D, v)
%   Derived from llfeelIC2b. Objective to provide a log likelihood for 
%   index participant for key set of params for the feeling ratings 
%   expressed by the Self, both towards self and towards Other, for a single new trial.
%   Version IC2b is the first version to use a PE kernel, following Rutledge, Will, Moutoussis.
%   To produce simulated data, P has to have field P.synth = [doDecs, doEvals]
%   NO inferred feelings of Other towards self. 
%   See end of file for testing instructions.
%
%   - Retrieves laboriously v.moves from D returns ...
%   - ps must have field evalRL, a vector with transformed entries for, maximally :
%          1      2       3      4      5       6     7     8      9         10     11    12
%       lambda   eta     wEx    sig    wOS11 wOS12 wOS21  wOS22  lambda2     lps    E0s   E0o
%     I will try to coalesce or fix as many of these as poss, esp. lambda2, the 
%     impact of PEs on evaluations, and lambda, the impact of PEs on learning values,
%     and also try fixing wOS22 to zero. Could even try lambda=eta (the evaluation decay / LR)
%
%  llfeelIC2b doesn't needs ps to have feelings-to-scale mapping parameters, 
%  like llfeelIC2 needed so that rat2resp(v.feelPr, ps.feelm, ps.feelu, p.MCQEvaLen) could be used.
%   
 
%%                                Outputs:
%   ll the log likelihood of v.Spol for this trial.   
%   v.nfeel     : (generated if need be) new feeling ratings. 
%                 Has 1->self, 2->other, 3->other-to-self, corresp. to col 11 to 13 in D.evo
%   - D.RLEval: 
%       Sets up, if v.trial==1, and fills in, an array with trial by trial
%       results of RL-ish variables of evaluation, D.RLEval and D.RLEvalHd
%   - various other fields of v

%  Maximum extend of evaluation / feeling param vector:
%          1      2     3      4     5      6     7     8       9       10   11    12
%       lambda   eta   wEx    sig   wOS11 wOS12 wOS21  wOS22  lambda2  lps   E0s  E0o
lambda = invlogit(ps.evalRL(1));    lambda2= invlogit(ps.evalRL(9));
eta    = invlogit(ps.evalRL(2));
wEx    = invlogit(ps.evalRL(3));
sig    = exp(ps.evalRL(4));
wOS = [ps.evalRL(5:6); ...  % Weights of returns on self and other on Self-Eval 
       ps.evalRL(7:8)];    % ... and of self ret. and other ret. on Other-eval.
% lapse rate:
try lps = invlogit(ps.evalRL(10)); catch lps=0;  end     
% Trait self- and other- evaluation constants:
try
    E0s = ps.evalRL(11); E0o = ps.evalRL(12);
catch
    E0s = 0; E0o = 0;
end

% Useful constant matrices, as the update will be 
% E(:,tr) = Wx*E(:,tr-1) + Wr*r(:,tr) ;
Wx  = (1-eta)*[[1-wEx, wEx]; [wEx, 1-wEx]];  
Wr  =    eta * wOS;

%% Some preliminaries
blN = v.block; 

% Retrieve / oorrect v.settingLev (which was a stupidly ambitious thing
% to vary in the experiiment ...) and v.moves, which MM forgot to 
% directly record in D.
v.settingLev = D.evo{blN}(P.evoPols{blN}(v.trial), 3);  % It's at third col, and all
             % rows dedicated to this trial in D.evo.
ncr = P.nCr;  % number of exchanges, or 'rounds', in this trial. In HT1 it's just 1.
for k=1:ncr
    sWorkInd= D.evo{blN}(P.evoSt{blN}(v.trial),8 ); % 8  is self-contribution
    oWorkInd= D.evo{blN}(P.evoSt{blN}(v.trial),21); % 21 is other-contribution 
    % was: find((P.ret(:,:,1,v.settingLev) == D.evo{blN}(P.evoSt(v.trial)+k-1,9)) & ...
    %                       (P.ret(:,:,2,v.settingLev) == D.evo{blN}(P.evoSt(v.trial)+k-1,10)));
    v.moves(:,k,v.settingLev) =  [sWorkInd, oWorkInd]'; 
end

% Initialize ...             
ll = [0 0];              % log-likelihoods
Z = max(abs(P.ret(:)));  % Constant for normalizing returns.
% Initialize the feeling ratings for this trial off the data D. 
% If synthetic data is to be generated, this will be over-written later, 
% after the ll is caclulated for this v.nfeel.
v.nfeel = D.evo{blN}(P.evoPols{blN}(v.trial)+1, 11:13);  % REM evoPols{blN}(v.trial)+1  
         % is the '902' code line with S/O/oSeval at col 11-13, S/Ofee at 22-23, etc.
feelReady = ~( sum(v.nfeel) == 0 ); % This will be false if running a completely
         % simulated / generated dataset, in which case evals will not be ready.
        
% REM: P.synth has two elements, the first for 'treating each other', the second
% for 'feelings' / evaluatios. If either is not 0 ...
try P.synth; catch, P.synth = [0,0]; end
fSynth = P.synth(1) || P.synth(2);  % Whether synthetic feeling ratings will be generated.

%% Begin main processing. If we are at the first trial, set up 
%  storage array and variables, some of which are quite non-trivial.
if v.trial == 1
    if blN == 1
        D.RLEval = {};
        D.meanRLEval = {}; 
    end
    % Make space in the data structure for an array with 
    % reinforcement learning variables of Evaluation. Headings:
    % trial, SEval decision variable, ..., Self Return expectation ... 
    %                                            simulate Self-Eval, ... ll inputted SEval ... ll Simulated SEval ...
    %                1       2          3         4          5        6        7      8        9         10         11
    D.RLEvalHd = {'trial','SEvalDV','OEvalDV','SretExp','OretExp','simSEv','simOEv','llSEv','llOEval','llSimSEv','llSimOEv'};
    D.RLEval{blN} = nan(P.trN(blN)+1,length(D.RLEvalHd));  % P.trN(blN) is the max. number of trials 
                                                       % for this partner / other (Helper in HT1)
    D.RLEval{blN}(:,1) = 0:P.trN(blN);  % each row is at the end of that trial, 
                              % so there is a 'zero'th row for initial values.
    [~, OAct0] = max(ps.prevPri);   % index of which action the participant believes is 
                                    % most prevalent (popular) in the population. 
                                    
    SAct0 = ps.Spref(v.settingLev, OAct0)+1; % the preferred action **index** of the pt. 
             % for the action index above, that they anticipate for the other, from 
             % the first trial, depending on the trial type.
   %  p.ret is indexed (myWork, yourWork, me or you, settingLev) :
   %  In HT1, p.ret(:,:,1,1) are the returns of self (help-seeker), 
   %      and p.ret(:,:,2,1) for Partner/ Helper / Other
   D.RLEval{blN}(1,4:5) = P.ret(SAct0,OAct0,1:2,v.settingLev) / Z;
   
   % Set the initial value of the evaluation as if the 
   % return above was a known equilibrium, with no prediction error:
   % (sadly the D.RLEval is in row vector form, the RHS below is 
   % naturally a col vec :/ , hence the transpositions...)
   % (line below is matlab faster version of 
   %  D.RLEval{blN}(1,2:3) = ( inv([[1 0]; [0 1]] - Wx) * Wr  *(1-lambda2)* D.RLEval{blN}(1,4:5)' )' ;    )
   D.RLEval{blN}(1,2:3)    = (   (([[1 0]; [0 1]] - Wx) \ Wr )*([E0s; E0o] + (1-lambda2)* D.RLEval{blN}(1,4:5)') )' ; 
   
   D.meanRLEval = nan(1,length(D.RLEvalHd)); % will hold means for trials 1-end
end

%% Retrieve actions of self and other, and resulting returns:
% rem v.moves had as row 1 the moves of the pt, row 2 the moves of the prtn.
act = v.moves(:,1:ncr,v.settingLev);
Ret = [0,0]';  % To hold average returns in this trial.
for k=1:ncr
   % REM p.ret is the return matrix, with indices
   % selfContr, otherContr, whoseReturn (1=self, 2=other), settingLev (1=ptTempted, 2=otherTempted)
    Ret = Ret + squeeze(P.ret( act(1,k), act(2,k), :, v.settingLev)); % accumulate. Normalizing etc. below.
end
Ret = Ret / (Z *ncr ) ;

%% feelings / evaluations:
%  First find prediction error (note vectorial form [self other]'):
PE = Ret - D.RLEval{blN}(v.trial,4:5)'; % D.RLEval{blN}(v.trial,4:5) is 
     % 'SretExp','OretExp' AT THE END OF (AS UPDATED BY) v.trial-1
D.RLEval{blN}(v.trial+1,4:5) = D.RLEval{blN}(v.trial,4:5) + lambda * PE' ;
% Now update the evaluations:
D.RLEval{blN}(v.trial+1,2:3) =  ( [E0s; E0o] + ...
                                  Wx * D.RLEval{blN}(v.trial,2:3)' + ...
                                  Wr*((1-lambda2)*Ret + lambda2*PE))' ; 

% Use the evaluation variables to find probability action vectors
% for the self-report measures. Take account of noise floor.
v.SfeelRespPr = logitBins4norm(D.RLEval{blN}(v.trial+1,2),sig);
v.SfeelRespPr = (1-lps)*v.SfeelRespPr  + lps/length(v.SfeelRespPr); 

v.OfeelRespPr = logitBins4norm(D.RLEval{blN}(v.trial+1,3),sig);
v.OfeelRespPr = (1-lps)*v.OfeelRespPr  + lps/length(v.OfeelRespPr); 

% Always calc. the log-lik if there is inputted eval/feelings data. 
% ll, v.Sfeel, v.Ofeel will be over-written if synthetic data is simulated.
if feelReady
    v.Sfeel = v.SfeelRespPr(v.nfeel(1));
    ll(1) = log(v.Sfeel);
    D.RLEval{blN}(v.trial+1,8) = ll(1);   % col llSev

    v.Ofeel = v.OfeelRespPr(v.nfeel(2));
    ll(2) = log(v.Ofeel);    
    D.RLEval{blN}(v.trial+1,9) = ll(2);   % col llOev
else
    D.RLEval{blN}(v.trial+1,8:9) = nan(1,2);
end

%%  ~~ if new evaluation data is to be simulated ~~
if fSynth   % Then create new v.nfeel - otherwise they  have
    % been estimated above.
    v.nfeel = [pBinSample(v.SfeelRespPr,1),0,0];  % last two zeros for mem. alloc.
    % store in D.evo, as in llfeelIC2 ...
    D.evo{blN}(P.evoPols{blN}(v.trial)+1, 11) = v.nfeel(1); % See above for P.evoPols{blN}(v.trial)+1
    D.RLEval{blN}(v.trial+1,6) = v.nfeel(1);  % ... and in D.RLEval, col simSEv
    v.Sfeel = v.SfeelRespPr(v.nfeel(1)); % probability of emitting the chosen self-eval.
    ll(1) = log(v.Sfeel);
    D.RLEval{blN}(v.trial+1,10) = ll(1);   % col llSimSev
  
    v.nfeel(2) = pBinSample(v.OfeelRespPr,1);
    D.evo{blN}(P.evoPols{blN}(v.trial)+1, 12) = v.nfeel(2);  % See above, ... 
    D.RLEval{blN}(v.trial+1,7) = v.nfeel(2);   % ... col simOEv. 
    v.Ofeel = v.OfeelRespPr(v.nfeel(2)); % probability of emitting the chosen self-eval.
    ll(2) = log(v.Ofeel);
    D.RLEval{blN}(v.trial+1,11) = ll(2);   % col llSimOev
end

%% if at last trial, also record the means:
if v.trial == P.trN(blN)
    D.meanRLEval = mean(D.RLEval{blN}(2:(P.trN(blN)+1),:)); 
end
return;  % end of whole function.

%%    Instructions how to test or demo this fn, llfeelIC2b  : 
% %               1         2         3            4        5      6     7     8       9       10          11    12
% %            lambda      eta       wEx          sig     wOS11 wOS12 wOS21  wOS22  lambda2    lps         E0s   E0o
%   evalRL =[logit(0.15),logit(0.1),logit(0.001),log(0.5), 0.33, 0.67, 0.67,  0, logit(0.0001),logit(0.001),0.00, 0.00 ] 
% cd('C:\Users\mmpsy\Dropbox\FIL_aux\MSc_iBSc_PhD_student\Gosalia_(Meera_iBSc)\SelfOtherEvalModelsProject-shared\SelfEvalModelCoding\dataIC2\resultsDec')
% load('fitIC2t_fmri.mat') 
% ptN=1; prtn=2; d1=D{ptN}.here.d; p1=D{ptN}.here.p;  ps1=ps{ptN,prtn}; ps1.evalRL=evalRL; p1.synth=[0,1];
% v1.trial=1; [ll,v2,d2] = llfeelIC2b(ps1, p1, d1, v1); for k=2:23; v1.trial=k; [ll,v2,d2] = llfeelIC2b(ps1, p1, d2, v1); end; 
% This should be Seval for real, then simulated data (can check d1.key_hd and d2.RLEvalHd for good measure).
% ind = find(d1.key_data(:,7)' ~= 0); plot([d1.key_data(ind,7)'; (d2.RLEval(2:end,6))']')
