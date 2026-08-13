function [v,D,P] = simHelperHT1(v,D,P)
%SIMHELPERHT1 generate Helper actions for simulating Taylor Kwan's verion of Helping Task.
%   Only new action for trial v.trial are.
%   Helper action are based on a simple probability table
%   The kind of Helper that is simulated should be encoded in both the Data 
%   (immediately) and the bird's eye view of the task, in P, e.g.
%    p.helperType:   2    10     <---- row 1 := block 1, 2=Uncaring, 10 trials.
%                    1    13     <---- row 2 := block 2, 1=Caring, 13 trials. 
%   p.helperPol(helpSeekerPrevWork,helperWorkNow,Uncaring=1,first half of task =1, 2nd=2)

% Retrieve (median of) previous play of help-seeker

% 
blN = v.block;
helperType =  D.evo{blN}(v.trial*P.nCr, 30);  % last row of block for current trial.

%% Select the policy, from which action will be chosen:
if v.trial == 1
   v.Ppol = P.helperInit(helperType,:);  % Partner's / other's policy
else  % i.e. not the first trial
   % In HT1, the Other's policy changes with the trials, the simplest form to attempt misleading.
   % So first find if we are in the first half of the trials:
   sect1TrN = floor(P.helperType(v.block,2)/2);  % Number of trials in the first half
              % of this block / partner 
   sect = 1+ (v.trial > sect1TrN ); % Which part of the helper policy array to use.
   % Select the corresponding matrix:
   v.Ppol = P.helperPol(v.moves(1), ...
                                   :, ...
                                   helperType,...
                                   sect);
end

%% Select the action from the policy :
if P.maxPsynth==1   % Again, if we have been instructed to enact 
                   % only the most probably policy :
  v.nmoves(2,1:v.crN,v.settingLev) = find(v.Ppol == max(v.Ppol),1);
else
  v.nmoves(2) = pBinSample(v.Ppol,1);
end


return;