
%%  TESTING BLOCK FOR FUNCTION BELOW
% clear runtype selfMod pS outp expt expll; 
% pS = expt25{1}.pS;    % expt25 is synthetic with same params and 25 pts
% or can do load process and set pS = p;
% % pS.SPartnPr = (discBetaMU(0.9,1,4));  
%pS.T = 0.1; 
%pS.prevPri =  discBetaMU(0.7,0.5,4);  % was [0.0880 0.6946 0.2075 0.0099]
%pS.Spref = [[1 0 1 3]; [1 1 1 2]];    % originally [[1 0 1 3]; [1 1 1 2]];
%pS.SPartnPr = discBetaMU(0.6,0.2,4);  % originally [0.2485 0.3111 0.2728 0.1675]
%pS.feelm = 2; pS.feelu = 0.75;     % feelings scale use parameters
%runtype.detailed=1; 
%runtype.tSynth=0;  runtype.fSynth=0;   runtype.maxPsynth=0; 
%expll = []; expt = {}; Sll =s[];

% sim a whole set of pts if need be:
%simN = 5;   % if ~runtype.synth; simN=1;  end;
%for ptN = 1:simN
  % may want to have different random seed for successive synthetic pts,
  % or comment this out :
%  runtype.rseed = 210622;             p.rseed = runtype.rseed + 7*ptN;
%  n2 = rem(ptN,25); if n2==0; n2=25; end;
%  d = expt25{n2}.dsynth;
%  disp(['Now on ' num2str(ptN)]); clear v outp ; 

%% comment out all above when generative testing complete, and also all
%  below 'endpiece' block of code at end of script / function :

function [sllt, sllf, outp] = HT1ll03(pS, d, p, runtype, selfMod)
% function HT1ll03(pS, d, p, runtype, selfMod)
%   (developped: IC2ll1 --> HT1ll1 --> this; relies on discBetaMU, not noisyBino.
%   sum log likelihoods for HT1 (Taylor's task0 data produced by one participant, interacting
%   with one partner only, derived from IC2fmri.m
%   'Helping Task' where participants are asked how much to commit to
%   the exchange. See also intro to IC2*.m
%   sllt is the sum-log-lik over the treating-each-other data,
%   sllf is the sum log like over the feelings (evaluations) data, with components
%         for Self, Other and (optional) inferred-Self-by-Other.
%   feelings data (the self report approval / interpersonal trust ....
%
%   pS is structure with self-parameters theta so P(d|theta,M) is calculated, 
%   whereas argument p \in M contains specification of the model/task dynamics.
% 
%   runtype is a structure potentially containing fields:
%    debugging
%    tSynth         Whether to run in syntetic-data-producing mode *rather than* real-data-ll-estimation mode
%                   for treating-each-other choices.
%    fSynth         As above for evaluations (feelings) choices.
%    synthPmax      Whether, if in synth mode, to use the peak-probability policy rather than sample.
%                   Hopefully useful for debugging ...
%    detailed       If not, only provide the ll (and basic synth data in synth mode)
%    rseed          If 'detailed' is not [], then try to use this to initiate random number generator.
%
%   outp contains a number of fields for synthetic data, but will always contain the 
%      contributinons data from the input argument d (which I 
%      forgot to encode in d during the expt, so re-calculated :/ )
%
% Instructions for running example / demo / debug after 'return' at end of function.
% Adapted for Helping Task from July 2026 on.


%% Default arguments. First one clears the decks too! 
try % check debugging flag
    runtype.debugging;  % check before proceeding !
    if isempty(runtype.debugging); runtype.debugging=0; end
    p.debug = runtype.debugging;
catch  
    runtype.debugging = 0; p.debug = runtype.debugging;
end
try % check synthetic data related flags
    runtype.tSynth;  % check before proceeding ! This is investment decisionmaking synth. flag.
    runtype.fSynth;  % This is feelings/evaluations synth. data flag.
    if isempty(runtype.tSynth); runtype.tSynth=0; end
    if isempty(runtype.fSynth); runtype.fSynth=0; end
    p.synth = [runtype.tSynth runtype.fSynth];
    try  % if maxPsynth=1, use the maximum-probability policy only to
         % generate synthetic data, rather than sampling from the policy prob. vector
      p.maxPsynth=runtype.maxPsynth; 
    catch
      p.maxPsynth=1; 
    end
catch  
    runtype.tSynth = 0;    runtype.fSynth = 0;  p.synth = [0, 0];
end
try % check detailed output flag
    runtype.detailed;
    if isempty(runtype.detailed)
       runtype.detailed=0; 
       p.rseed = sum(100*clock);
    else
      try  % if detailed output, try to use provided 
           p.rseed = runtype.rseed;
      catch
           p.rseed = [];
      end
      if isempty(p.rseed); p.rseed = round(sum(100*clock)); end
    end
    p.detailed = runtype.detailed;
catch
    runtype.detailed=0;  
    p.detailed = runtype.detailed;  p.rseed = round(sum(100*clock));
end
try % selfMod here is to try different generative models for Self.
    % selfMod(1) is re. treat-each-other model, default in llSpolIC2
    % selfMod(2) is the Evaluation model, 1 being belief & regret based, 2 being RLish 
    p.selfMod = selfMod; 
catch % set to empty for extra ease of testing.
    p.selfMod=[];
end
% if pS is empty, use the 'real' (fitted) or 'default' ones from d.
% For HT1, defaults are set in csv2pdHT1b
if isempty(pS)
  pS.Spref =  d.Spref;        pS.prevPri = d.prevPri;       
  pS.SPartnPr=d.PartnPr;      pS.T = d.T;
  pS.blockLR =d.blockLR;
  pS.feelm = d.feelm;         pS.feelu= d.feelu;     pS.evalRL = d.evalRL;
end
 
%% ******************** PARAMETERS  & DATA STRUCTURES **********************
% Cater for prehistoric matlab :
p.llmatVer = version; % matlab version used
if isempty(str2num(p.llmatVer(1:4)))
    p.llmatVer = str2num(p.llmatVer(1:3));
else
    p.llmatVer = str2num(p.llmatVer(1:4));
end
if p.detailed || p.synth(1) || p.synth(2) % don't bother if we don't have to sample random numbers ...
  if ( p.llmatVer > 7.13)  % Dell laptops have 7.14
    rng('default'); rng(p.rseed); % reset random number generator & reinitialize it via clock.
  else % for older versions ...
    RandStream.setDefaultStream(RandStream('mt19937ar','seed',p.rseed));
  end
end


%% Initialize v, outp, dsynth ------------------------------------------------
sllt=0;  sllf=[0,0,0];  outp=[];  % initialize outputs of whole fn.
totBlN = size(p.helperType,1);

if runtype.detailed ~= 0  % matrix with feelings / approvals related measures
    % inferred on basis of parameters: data (MCQ answers, ds for self, do for other etc),
    % log-lik of ratings but also predicted approval scores
    % e.g. likelihoods of actions under the participant's ideal behaviour in the circumstances:
  outp.feelPhd = {'Seval','Oeval','oSeval',...   % 1-3: new/generated evaluation/feelings ratings
      'lls','llo','llos',...                     % 4-6
      'prAppS1','prAppS2','prAppS3','prAppS4','prAppS5','prAppS6',...   % 7 -12
      'prAppO1','prAppO2','prAppO3','prAppO4','prAppO5','prAppO6',...     % 13-18
      'prAppOS1','prAppOS2','prAppOS3','prAppOS4','prAppOS5','prAppOS6'};   % 19-24
  for rn = 1:p.nCr; outp.feelPhd{end+1} = ['Scontr' num2str(rn)]; end
  for rn = 1:p.nCr; outp.feelPhd{end+1} = ['Ocontr' num2str(rn)]; end
  % was:  'Scontr1','Scontr2','Scontr3', 'Ocontr1','Ocontr2','Ocontr3'};       % 25-30
  outp.feelPrat = {}; 
  outp.feelPrat{1} = zeros(p.trN(1),length(outp.feelPhd));
  outp.feelPrat{2} = zeros(p.trN(2),length(outp.feelPhd));
  outp.contrEval = {};
  % Prepare to make a record of the work contributions
  for block = 1:totBlN
      outp.contrEval{block} = d.evo{block}(:,[1,2,3,9:13]);  % These are the initial 3 cols with trial number and type etc,
       % and then, merely as template, the retS, retO and *Eval columns, the ret* columns
       % to be replaced by sContr and oContr
  end
  contrCol = 4:5;  % columns in outp.workD where sContr
  outp.contrEvalHd = {'trial','round','trType','sContr','oContr','sEval','oEval','osEval'};  
end

% Initialize v ---------------------------------------------------------------
v.assets = [0,0]; % Amounts accumulated by S[elf] and P[ar]T[ner] respectively.
v.block = 0;
% moved into blocks loop below: v.trial = 0; % current trial (aka game)
v.practTrDone = 0; % number of practice trials done 
     % - when this goes to pars.practiceTrN, real task begins ...
     % adjusted separately if genMod is provided, in prep4PractOrMain...
v.Spol = nan(1,p.Nl);   % Self policy (as opposed to self action)
v.Ppol = nan(1,p.Nl);   % Partner (other, here Helper) policy (as opposed to action).
v.settingLev = p.invalid;           % current setting / context level
% Here, a number of additional v.ariables are set or reset at the start of every new partner / block.
v.demo = 0;   % never in demo mode in this log-lik function


%% extract in convenient format:
% pols = d.evo{block}(d.evo{block}(:,2)==901,[4 20]); % pols(trN,1) is self Spol and pol(trN,2) is Opol played.
% REM p.evoSt;  % starting row indices of each trial in d.evo
if ~isfield(p,'evoPols')  % if not already calculated ...
    p.evoPols = cell(1,3);
    p.evoEvalFee = p.evoPols;
    for k=1:length(p.evoSt)
        p.evoPols{k} = p.evoSt{k} + p.nCr;   % the '901' code lines with 
                                 % Spol at col 4 and Opol at 20
        p.evoEvalFee{k} = p.evoPols{k}+1;     % the '902' code lines with 
                           % S/O/oSeval at col 11-13 and S/Ofee at 22-23
    end
end

% If synthetic data is to be produced, prepare dsynth ---------------------------------
dsynth = d; % d-synthetic starts as copy of d; and remains so if we are not synthesizing data.
dEvalCol = 11:13; % reminder of cols in d.evo which  should be {'Seval'} {'Oeval'} {'oSeval'} 

if p.synth(1) % this is about synthesizing choices about treating-each-other, 
    % and also using in which case the synthetic feelings are also calculated.
    % Now zero the to-be-generated bits of dsyn : 
       dsynth.evo{block}(:,[4:19,21:23]) = 0;
end
if p.synth(2)   % this is about synthesizing choices about feelings (evaluations).
    % If only this is required, keep all the existing data - llfeelHT1b will
    % need it.
    % (Change below for new model version)
    if p.selfMod(2)==1 
       % Now zero the to-be-generated bits of dsyn : 
       dsynth.evo{block}(:,dEvalCol) = 0;   
    end
    if  p.selfMod(2) >=2 && p.selfMod(2) > 15
       error(['not ready for p.selfMod(2) ==' num2str(p.selfMod),',p.synth=', num2str(p.synth)]);
    end    
end

%% ********************** LOOP THROUGH EXPERIMENT **************************************
% Loop over blocks
for block = 1:totBlN
    v.trial = 0 ;     % Will be augmented near top of while loop over trials below.
    v.block = block;  % Store copy to pass to functions.
    v.moves = nan(2,p.nCr,p.settingl);  % will store most recent moves. mid dim is not so i.
              % (also see v.nmoves, auxiliary for moves newly created, before they are stored in v.moves)
    v.prMoves = v.moves;                % moves of trial before last played
    v.rets  = v.moves;
    v.prRets = v.rets;
    
    %%  Loop over trials - - - - - - - - - - - - - - - - - - - - - - - - - 
    % REM p.trN is in order of the blocks to be performed sequentially.
    while v.trial < p.trN(block)+p.settingl  % Note we are doing one set of contexts past the end,
      % so as to 'look back' at the last instance of each context and calc. various 
      % measures - which, unfortunately, I build into the updating that takes place
      % when each new decision is taken ...
    
      % Orientation - which trial and what kind. In HT1, this is for 
      % future reference, as there is only one type of trial.
      v.trial = v.trial + 1; 
      if v.trial <= p.trN(block)
          v.settingLev = p.setting{block}(v.trial); 
          v.crN = p.nCr;   % historical remnant from versions where v.crN could be variable ...
      else
        v.settingLev = v.trial - p.trN(block);  % so 1,2,3, etc. 
        % leave v.crN at whatever value it last had ...
      end
        
      %% retrieve actions performed by other (and self if not in syntetic mode)
    
      % REM: columns of d.evo, as per d.evo_hd, are:
      % Columns 1 through 12
      %  {'trial'}  {'round'}  {'trType'}  {'Spol'}  {'tStim1'}  {'tDec1'}  {'RT'}  {'sWork'}    {'retS'}    {'retO'}    {'Seval'}    {'Oeval'}
      %
      % Columns 13 through 22
      %  {'oSeval'}    {'infOpola'}    {'infOpolb'}    {'tScontrOn'}    {'tScontrOff'}    {'tOcontrOn'}    {'tOcontrOff'}    {'Opol'}    {'oWork'}    {'Sfee'}
      %
      % Columns 23 through 32
      %  {'Ofee'}    {'tSevalOn'}    {'tSevalOff'}    {'tSeval'}    {'tOevalOn'}    {'tOevalOff'}    {'tOeval'}    {'<blank>'}    {'slStim1'}    {'slDec1'}
      %
      % Columns 33 through 40
      %  {'slScontrOn'}    {'slScontrOff'}    {'slOcontrOn'}    {'slOcontOff'}    {'slSevalStim'}    {'slSeval'}    {'slOevalStim'}    {'slOeval'}
    
      %  in this trial:
        %  v.moves contains previous moves and v.nmoves contains 'new moves'
        %  v.nmoves(1=self or 2=other, round in trial, settingLev))  
        %  v.nmoves had as row 1 the moves of the pt, row 2 the moves of the prtn.
        %  p.ret is indexed (myWork, yourWork, me or you, settingLev)
        % Self and other return in columns 9 and 10 of d.evo, so:
      
      for rn=1:v.crN       % to update moves for every play (exchange) in this trial
                           % This is trivially 1 in the first version of Helping Task.
        if v.trial <= p.trN(block) 
          sWork= d.evo{block}(p.evoSt{block}(v.trial)+rn-1,8 ); % See above for 8. 
          oWork= d.evo{block}(p.evoSt{block}(v.trial)+rn-1,21); % See above for 21. 
          % was: find((p.ret(:,:,1,v.settingLev) == d.evo{block}(p.evoSt(v.trial)+rn-1,9)) & (p.ret(:,:,2,v.settingLev) == d.evo{block}(p.evoSt(v.trial)+rn-1,10)));
        else % use last of this type actually done, again.
          sWork = nan; oWork = nan;   % if past the end of actual trials, 'flag' values.
        end
        if isempty(sWork); error(['Could not retrieve v.nmoves for trial ', num2str(v.trial)]); end
        if (runtype.detailed ~= 0  && v.trial <= p.trN(block) )
           evoRow =  (v.crN+2)*(v.trial-1)+rn; % use same row of outp.contrEval as in d.evo
           outp.contrEval{block}( evoRow, contrCol ) = [sWork, oWork]; 
           if rn ==1  % also store the evaluations from the given data
              % So this should be stored in e.g. the last row for trial 1.
              outp.contrEval{block}( evoRow+v.crN+1 , contrCol(end)+(1:3) ) = d.evo{block}(evoRow+v.crN+1, dEvalCol ); 
           end
        end
        if p.synth(1)  % if in treating-each-other-choices-synthetic mode, self acts will be inferred 
                       % below by llSpolIC2, so here set to nan :
           v.nmoves(:,rn,v.settingLev) = [ nan , nan]; 
           % Now simulate Helper nmove :
           [v,dsynth,p] = simHelperHT1(v,dsynth,p);
        else
           v.nmoves(:,rn,v.settingLev) = [sWork, oWork]; 
        end
      end
      
      %% find policy probability for moves in this exchange based on moves already done,
      %  prior beliefs etc. and if necessary make synthetic data,
      %  keeping all the actions of the Other constant. If in synthetic mode,
      %  HT1llSpol03 also updates the Self v.nmoves. Onwards for self move : 
      if ~p.synth(1)
        if v.trial < p.trN(block); v.Spol = d.evo{block}(p.evoPols{block}(v.trial), 4);   % v.Spol needs to be
                  % retreived from d, not dsynth, as Spol has been wiped in the latter.
        else v.Spol = nan;  end  % flag value  'past the end'
      end
      % line below should cater, in version HT1ll03, for p.synth having
      % two entries, 1st for everything, 2nd for feelings/evals only.
      [ll, v, dsynth] = HT1llSpol03(pS, p, dsynth, v);  % rem ll has ll for v.Spol
      sllt = sllt + sum(ll); % this works even for the dummy trials at the end.
      
      %% update v.moves
      if v.trial <= p.trN(block)
         v.moves(:,1:p.nCr,v.settingLev) = v.nmoves(:,1:p.nCr,v.settingLev);
      end
      %% do the eval / feeling thing, if feeling related params are present! 
      if v.trial <= p.trN(block)   % I *think* that's right ...
          %  llf is 2 or 3-component, for Self->Self, Self->Other, and poss. inf_Other->Self
          if (isfield(pS,'feelm')) || (isfield(pS,'evalRL'))
            % if ~isempty(pS.feelm) && ~isnan(pS.feelm)  % This check whould not
            % be needed if p.selfMod has been provided OK as below ...
            % (Change below for new model version)
              if p.selfMod(2) == 1
                 % ll for feeling ratings for basic joint pref. model only
                 [llf, v, dsynth] = llfeelHT1(pS, p, dsynth, v); 
              elseif p.selfMod(2) >= 2  && p.selfMod(2) <= 15
                 % RLish / PE kernel based approval / emotion. Upper limit raised
                 % from 14 to 15 for RLEval15: without this, selfMod(2) = 15 fell
                 % through to the error below, which the try/catch in HT1fFit*
                 % swallowed as 'fmincon MAP fitting failed' for every participant.
                 [llf, v, dsynth] = llfeelHT1b(pS, p, dsynth, v); % rem ll has ll for feeling ratings only
              else
                 error(['p.selfMod(2)==' num2str(p.selfMod(2)) ' not catered for.']);
              end           
              sllf = sllf + [llf 0];   % the last 0 is for backwards compatibility - it's for the sOEval ll
            % end
          else
            sllf = nan;  
          end
          % REM outp.feelPhd = {'ds','do','dos','lls','llo','llos', ...
          %                     'prAppS','prAppO','prAppOS', ... x 6 each actually
          %                     's1','s2','s3','s4','s5','o1','o2','o3','o4','o5','so1',...};
          if runtype.detailed ~= 0 && isfield(pS,'feelm')  % 
            % before filling in outp.feelPrat, check that oSfeel* ('how much does your partner
            % approve of you' is valid - in fMRI and HT1 versions, this was not asked for :
            try v.oSfeel;       catch v.oSfeel = nan*v.Sfeel;               end
            try v.oSfeelRespPr; catch v.oSfeelRespPr = nan*v.SfeelRespPr;   end
            outp.feelPrat{block}(v.trial,:) = [v.nfeel,...  % new / generated feelings ratings (2 of them)
              mean(log(v.Sfeel)),mean(log(v.Ofeel)),mean(log(v.oSfeel)),...    % log-liks 
              v.SfeelRespPr,v.OfeelRespPr,v.oSfeelRespPr, ...                  % pmfs for responses
              v.moves(1,1:p.nCr,v.settingLev), ...        % moves (effort / treating each ...
              v.moves(2,1:p.nCr,v.settingLev)];           % ... other actions.
          end
          
      end % do feeling thing if v.trial <= p.trN(block)
      
    end % loop over trials
end % loop over blocks

if runtype.detailed  % provide summary of eval and average obs. data
   % % To summarize into one line:
   % blockLen = 5; % for Seval,Oeval,oSeval,sContrAv,oContrAv
   % % Prepare the header:
   % keyExpHd = cell(1,p.trN*5+3);
   % keyExp = nan(1,length(keyExpHd));
   % keyExpHd(1:3) = {'llSevalPerTr','llOevalPerTr','lloSevalPerTr'};
   % keyExp(1:3) = mean(outp.feelPrat(:,4:6));  % mean log likelihood per trial
   % for trN = 1:p.trN
   %     trStr = num2str(trN);
   %     keyExpHd( (4+(trN-1)*blockLen) : (3+trN*blockLen) ) = ...
   %         {['sContrAv' trStr],['oContrAv' trStr],...
   %          ['Seval' trStr],['Oeval' trStr],['oSeval' trStr] };
   %     rows = (trN-1)*(p.nCr + 2) + (1:p.nCr);   % rows of d.evo and outp.contrEval
   %     keyExp(4+(trN-1)*blockLen)    = mean(outp.contrEval(rows,4));
   %     keyExp(5+(trN-1)*blockLen)    = mean(outp.contrEval(rows,5));
   %     keyExp((6:8)+(trN-1)*blockLen)  = outp.contrEval(rows(end)+2,6:8);
   % end
   % outp.keyExp = keyExp;
   % outp.keyExpHd = keyExpHd; 
   outp.sllt = sllt;      outp.sllf=sllf; 
   outp.sll = sllt + sum(sllf);
end

if runtype.detailed && ~(p.synth(1) || p.synth(2))
  outp.llreal = v.llreal;
  outp.d = dsynth;       % if no synthetic choices at all, the interest here is not
         % in such synthetic choices, but in the evolution of beliefs about self and other
         % stored in  dsynth.OthInfdC etc ... 
end

if (p.synth(1) || p.synth(2))
   outp.p = p;            outp.dsynth = dsynth;    outp.pS = pS;  
   % % Further summarize into one line
   % blockLen = 5; % for Seval,Oeval,oSeval,sContrAv,oContrAv
   % % Prepare the header:
   % keyhd = cell(1,p.trN*5+3);
   % key = nan(1,length(keyhd));
   % keyhd(1:3) = {'llSevalPerTr','llOevalPerTr','lloSevalPerTr'};
   % key(1:3) = nan(1,3); % was: mean(outp.feelPrat(:,4:6));  % mean log likelihood per trial for data provided not synth.
   % for trN = 1:p.trN
   %     trStr = num2str(trN);
   %     keyhd( (4+(trN-1)*blockLen) : (3+trN*blockLen) ) = ...
   %         {['sContrAv' trStr],['oContrAv' trStr],...
   %          ['Seval' trStr],['Oeval' trStr],['oSeval' trStr] };
   %     key(4+(trN-1)*blockLen)    = mean(outp.feelPrat(trN,25:27));
   %     key(5+(trN-1)*blockLen)    = mean(outp.feelPrat(trN,28:30));
   %     key((6:8)+(trN-1)*blockLen)  = outp.feelPrat(trN,1:3);
   % end
   % outp.feelKey = key;
   % outp.feelKeyHd = keyhd; 
   % was, before dsynth was used: outp.synth = v.synth;   
   outp.sllt = sllt;      outp.sllf=sllf; 
   outp.sll = sllt + sum(sllf);
end

%% ***************************** End of main function **************************
return; 
%%  % Instructions for running example / demo / debug 
% In Dell G7: datDir = 'C:\Users\mmpsy\Documents\sci\student_archive\Gosalia_(Meera_iBSc)\SelfOtherEvalModelsProject-shared\SelfEvalModelCoding\dataIC2\resultsDec\';
% Or, in Cedric: datDir = '/home/mmoutou/Data/experiments/NeuroCompEvalOfOthers/Gosalia_Meera_related/SelfOtherEvalModelsProject-shared/SelfEvalModelCoding/dataIC2/resultsDec/' ;
% %               1         2         3           4       5      6     7     8       9
% %            lambda      eta       wEx         sig    wOS11 wOS12 wOS21  wOS22  lambda2
%   evalRL =[logit(0.15),logit(0.1),logit(0.2),log(0.5), 0.33, 0.67, 0.67,  0,   logit(0.30) ] 
%
% load('HT1DemoArguments.mat'); ptN=1; prtn=2;
%
% (Or, alternatively:
%  runType.tSynth = 0; runType.fSynth=1;  runType.debugging = 1; runType.detailed = 1; runType.synthPmax=0; runType.rseed = 123;
% selfMod = [1,2]; % [1 for standard mentalizing DM, 2] for RLish / PEkernel Eval model
%
% [sllt1, sllf1, outp1] = HT1ll2(ps1, d1, p1, runType, selfMod)

%% endpiece if simulating whole expt:
%if p.synth
%  expt{ptN} = outp;
%  expll(:,:,ptN) = outp.synth;
%  Sll(ptN)= sum(expll(:,1,ptN)); 
%  if ptN==simN
%    disp('Sll:  '); disp(num2str(Sll)); hist(Sll); 
%    disp(['Avg. predictability:  ' num2str(exp(mean(Sll)/p.trN))]);
%  end
%else
%  Sll(ptN) = sll;
%  expSll(ptN) = expt25{n2}.sll;
%  if ptN==simN
%    disp(['Sll:  ' num2str(Sll)]);
%    disp(['Mean fitted sll: ' num2str(mean(Sll))]); 
%    disp(['Mean original sll: ' num2str(mean(expSll))]); 
%  end
%end

% end

