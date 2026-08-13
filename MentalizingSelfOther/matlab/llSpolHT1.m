function [ll, v, D] = llSpolHT1(ps, P, D, v)
% [ll, v, D] = llSpolHT1(ps, P, D, v)
%   Derived from ll SpolIC2, which was in turn derived from
%     PpolIC1 - note, not directly ...IC2 as this is meant to operate
%     on merged data that is all in one block, unlike the task itself which
%   is in 2 blocks etc. Objective is still to provide a log likelihood for 
%   index participant (not Other-avatar, like Ppol* ) for key set of params
%   pertaining to Self, for a single new trial. 
%
%   find the ll of Self's preference based on priors and evidence so far ...

%% 
%                                Outputs:
%   ll the log likelihood of v.Spol for this trial.   
%   v.SpolPr    : policy probabilities for Self for this trial
%   D.pri4OthC
%   D.CbOth
%   D.SelfPersona{blN}(:,:,v.trial)  records v.whoSis for all trials.
%   D.minOthC
%
%   v.othExPol  
%   v.mapOthC  
%   v.currmapOthC{cx} working copy to be used within this fn. only for Dkl with v.mapOthC just 
%                    before copying v.mapOthC into the correct context of this one.
%  
%         ps to be the parameters on which the ll depends, as opposed to 
%              task specifications that are drawn from P.
%    ps.Spref      % self preferences
%    ps.prevPri    % prevalence priors about the Other that pt holds.
%    ps.SPartnPr   % preferences over types of partners
%         task parameters that need to be provided
%    P.settingl
%    P.peaks4c
%    P.polU
%    P.Nl
%           variables that need to be provided but which will be updated
%    v must contain uptodate v.moves
%    v.moves    % remain as is, but change use  of 1st index from 2 (other) 
%               % in this function compared to PoplIC1 to 1 (self) 
%    v.detail   % If this is nonzero, produce full details including synthetic self-choices
%    v.rseed
%    v.whoSis   % initiated in this fn. only if no previous rounds in this context
%    v.whoSwas  % ditto

%% Check block number and if at the first block and trial, declare multi-block objects 
blN = v.block;   % convenience copy for block (partner)
%  Other general preliminaries
patN = size(P.peaks4c,1);          % number of patterns in P.peaks4c

if v.trial ==1  && blN ==1   % declare some cell arrays (kinda lists)
    D.SelfPersona = {}; 
    if P.detailed
       D.surp = {}; 
       D.surp_hd = {'tr','cr','PEmyact','lnPmyact','wPEself','Dklself',...  % 1 to 6, incl surpr about self
                  'whoSisP1','whoSisP2','whoSisP3','whoSisP4',...         % 7 to 10
                  'PEothact','lnPothact','DklmapOthC','DklothC',...      % 11 to 14, surpr about other.
                  'exOthWrk0','exOthWrk1','exOthWrk2','exOthWrk3',...      % 15 to 18, expectn. of other action
                  'DklContext'};                                          %  19: intial surprise when new context shown 
    end
    D.surpSt = {}; % for indices of first exchange / interaction in each trial
    D.pri4OthC = {}; 
    D.OthInfdC = {}; 
    D.minOthC = {}; 
    D.SelfC = {};
    D.rwiseCs = {}; 
    D.ideals = {}; 
end

%% ~~~~~~~~~  Section 1 - Self infers their own public persona ~~~~~~~~~~~~~~~
%  Initially it's their own prevalence prior, later 'who their actions make them' 
%  Does not try to invert Other-agent's mental  process as such, but instead Self
%  assumes that the prevalence beliefs about types incorporate all relevant
%  influence, including peoples social and monetary preferences when face
%  with the context (returns matrix/ces) at hand. v.whoSis has P.settingl rows
%  as Other-agent (and pt) are allowed to be different in different contexts.
%  cf.  1.4.2. in WhoToBeWithYou_2

if v.trial ==1    % for initialization of v.whoSis etc.
  v.whoSwas = repmat(ps.prevPri,P.settingl,1);  % ps.prevPri is prevalence 
  % priors about the Other that pt holds. d.prevPri etc. are 1 x P.Nl vec over types  (locn. of peak)
  v.whoSis = v.whoSwas;
  D.SelfPersona{blN} = nan(size(v.whoSis,1),size(v.whoSis,2), ...
          P.trN(blN)+P.settingl);  % storage array for whoSis at each trial. 
             % Note that it has more rows than there are trials, as the last P.settingl ones say 
             % whoSis at the very end, at the start of what the next decision would have been.
  % if P.detailed
  D.surp{blN} = zeros(P.trN(blN)*P.nCr,19); % +(P.Nl*P.Nl)); % storage of surprise, PE and 
           % update measures for PE, Dkl etc. regressors. Poss. last P.Nl*P.Nl block for currmapOthC.
  D.surpSt{blN}  = [0 cumsum(P.nCr*ones(1,P.trN(blN)-1))] +1;  % indices of first crossing for each trial
  v.currmapOthC = {};   % This will store the previous v.mapOthC for each context so we can do Dkls
  % end
end

if ~isnan(v.moves(1,1,v.settingLev))  % if there is a preceding exchange to take 
  % into account, update on the basis of contributions made.
  % rem v.moves had as row 1 the moves of the pt, row 2 the moves of the prtn.
  % Will use standard P.poldens set in paramIC1
  v.whoSwas = v.whoSis; % keep un-updated copy to infer other's C
  % v.moves has actions from the last time each context was encountered,
  sll = nan(1,P.Nl);   % to store sum log likelihood
  if ~P.detailed   % we updated every time decision & feelings are needed, to fit answers given.
    for me=1:P.Nl
      sll(me)  = sum(log(P.poldens(me,v.moves(1,:,v.settingLev))));
    end
    v.whoSis(v.settingLev,:) = exp(sll) .* v.whoSis(v.settingLev,:);  % un-normalized
    % ... and normalised belief *about the peak* of who I am in this context:
    v.whoSis(v.settingLev,:) = v.whoSis(v.settingLev,:) / sum( v.whoSis(v.settingLev,:) );
  else % we want interaction-by-interaction updates for imaging regressors etc.
    swhonew = v.whoSis(v.settingLev,:);   % auxiliary working copies 'who S is'
    % All will be derived from last trial at the present settingLev, so find where to store:
    if v.trial <= P.trN(blN)  % the dummy trials at the end treated
       % separately as P.settings isn't valid past P.trN(blN)
       settingLevInds = P.setting{blN}(1:(v.trial-1))==v.settingLev; 
    else
       settingLevInds = P.setting{blN}(1:P.trN(blN))==v.settingLev;
    end
    lastSuchTr    = find(settingLevInds,  1,'last');
    lastNotSuchTr = find(1-settingLevInds,1,'last');
    for cr = 1:P.nCr  % it's the trial BEFORE this decision that the evidence comes from!
      swhoold = swhonew;  % reference copy of  self image to calc. surpise soon            
      % here done crossing by crossing, so 'other way around', vectorially over 
      % my possible actions. First, likelihood of different actions:
      lmyact = P.poldens(:,v.moves(1,cr,v.settingLev))' ; 
      % then, posterior 
      swhonew =  lmyact .* swhonew;   % un-normalized
      swhonew = swhonew / sum(swhonew); % normalized
      % record indices and *new* 'whoSis' :
      D.surp{blN}( D.surpSt{blN}(lastSuchTr)+cr-1,[1:2, 7:10])  = [lastSuchTr, cr, swhonew];
      % Dkl about self-type peak and (weighed) PE about self-type peak:
      D.surp{blN}( D.surpSt{blN}(lastSuchTr)+cr-1,6) = rowpDkl(swhonew,swhoold );
      % 1:4 below for consistency with v.moves coding
      D.surp{blN}( D.surpSt{blN}(lastSuchTr)+cr-1,5) = dot(1:4,swhonew) - dot(1:4,swhoold);
      % lnP and PE about self action - just take into account peak 
      % belief about self type and emission noise:
      pmyact = swhoold * P.poldens;  % probabilites of my different actions just due to emission noise
      % if v.trial >= 6;  % debug lines
      % disp('now at tr >= 6');  end;
      D.surp{blN}( D.surpSt{blN}(lastSuchTr)+cr-1,4) = log(pmyact(v.moves(1,cr,v.settingLev)));
      D.surp{blN}( D.surpSt{blN}(lastSuchTr)+cr-1,3) = v.moves(1,cr,v.settingLev) - dot(1:4,pmyact);
    end
    v.whoSis(v.settingLev,:) = swhonew;
  end

end
D.SelfPersona{blN}(:,:,v.trial)  = v.whoSis; % always record whoSis per trial (as opposed
                                             %  to each interaction, if several)

%% ~~~~~~~~~~~~~   Section 2 - beliefs about Other ~~~~~~~~~~~~~~~~~~
%% Section 2a -  If on first trial, initialize beliefs about C of Other: 
%  v. weakly egocentric within a large repertoire of possibilities
if v.trial ==1
  % Construct prior bel. about peaks in conditional and marginal parts of 
  % C of Other, D.pri4OthC{blN}, with pages for all 'setting' contexts, and with -
  % - an almost flat conditional part of C, with slight 'egocentric' peaks. 
  %   Little evidence should be able to overwhelm it!
  % - marginal along the 'if Other ...', 'horizontal' edge of pref. that
  %   Index pt is believed to want Other-agent to have, here ps.SPartnPr . 
  D.pri4OthC{blN} = zeros(P.Nl,P.Nl+1,P.settingl);
  for cx=1:P.settingl
     for i=1:P.Nl    % note no mult. by the marginal pref.
        D.pri4OthC{blN}(:,i,cx) = noisyBino((ps.Spref(cx,i)+1)/(P.Nl+1),3*P.Nl,P.Nl);
     end
     D.pri4OthC{blN}(:,P.Nl+1,cx) = ps.SPartnPr;
  end
  
  % Also construct the matrix in which beliefs about the Other's C matrix will be 
  % stored. D.OthInfdC{blN} to contain the posteriors at each trial, v.oInfdC just the last.
  % So D.OthInfdC{blN} will have the correct initial beliefs for the beginning of the
  % current trial, as it's posterior to the moves currently considered, which are
  % the ones up to now, not the ones in the current trial!
  % Note that in PpolIC1 much was recorded in d structure which is now stored
  % in D instead.
  D.OthInfdC{blN} = nan(patN,P.settingl,P.trN(blN));
  D.minOthC{blN}  = 0.0001/patN;         % additive 'noise' to avoid zero priors.
  v.oInfdC = nan(patN,P.settingl);

end    % preparations required in first trial in all blocks
%
% Very first time around, also store explicitly the repertoire of Cpt matrices that
% that Self will consider, C_beliefs_about_Other....
if v.trial ==1 && blN ==1
  % ... NB this is a repertoire so we don't need different blocks for different contexts /
  % setting levels. Note also that this is a distribution over actualized
  % 'who-we-ares', not just about our respective peaks.
  D.CbOth = nan(P.Nl,P.Nl,patN);
  for ci=1:patN
      % First construct the marginal out of the last col. of peaks4c. uOC is
      % a by-hand value for the default / initial sharpness of the distribution of 
      % the unknown Other's preferences.
      marPref = noisyBino(P.peaks4c(ci,P.Nl+1)/(P.Nl+1),P.uOC,P.Nl);
      for i=1:P.Nl   % this now counts rows of D.CbOth, so that the 
         % putative C matrices for the pt. have abcissa (cols) Other and 
         % ordinate (rows) pt :
         D.CbOth(i,:,ci) = noisyBino(P.peaks4c(ci,i)/(P.Nl+1),P.uOC,P.Nl) * marPref(i); 
      end
  end
  D.CbOth2 = reshape(D.CbOth,P.Nl*P.Nl,patN); % convenint 2-d version of same.
  
end    % preparations required in first trial. 


%% Section 2b - infer upon C for other & find MAP C for other  ---------------------------
%  Self (pt) will then infer how Other is likely to be with them, based on MAP update of Other-C ~~
%
%  I.e. find MAP C for other amongst D.CbOth and then consider how such a person
%  will treat the Self that Other infered Self to be, v.whoSis. 
%  Self knows starting preval. beliefs about 'how people are' and trial evidence of 
%  Other behaviour. Self simulates Other's (active) inference on the basis of
%  preference matri(ces) Cprtn. 
%
%  In  1.4.3. Starting from egocentric conditional pref. and try to update
%             Beliefs About the peaks of the joint pref of the Other, as seen 
%             by the Self - 'If so far Other believed X about Self and 
%             responded y, what Y is likely to underpin Self's (joint) tastes?
%             -- see also WhoToBeWithYou_2 

%     - - - - - - - - find MAP C among D.CbOth - - - - - - -
% (Self beliefs for Other C matrix). Update beliefs about all but store the
%  MAP estimate in v.mapOthC.
   
% prelims:
Spers =  v.whoSwas(v.settingLev,:); % convenience copy of Self *previous* persona
Spolb =  sum(repmat(Spers(:),1,P.Nl) .* P.poldens); % choice beliefs about ...
           % ... the Self, according to (marginalised over) their persona.
% For most trials, the prior is the most recent posterior; except the first of
% each kind, where we need to get prior belief over D.CbOth first for all contexts,
% not just the one that the partners will face. This is derived
% from the D.pri4OthC{blN}(:,i,h) calculated at v.trial==1 above.
if isnan(v.moves(1,1,v.settingLev)) %i.e. initialize if no previous data to infer upon.
   % We will need all 'context' columns of v.oInfdC, the prior/most recent posterior,
   % so that we can also calc. the surprise about facing a context as opposed to
   % another at the beginning of a trial, if there are more than one contexts (not in HT1). 
   if sum(isnan(v.moves(1,1,:))) == P.settingl  % i.e. nan moves for all contexts / settings levels
        if blN == 1
          for cx=1:P.settingl
             for ci=1:patN
                prd = 1.0;  % auxiliary to accumulate product.
                for j=1:(P.Nl+1)
                   prd=prd*D.pri4OthC{blN}(P.peaks4c(ci,j),j,cx); 
                end
                v.oInfdC(ci,cx) = prd;
                if ci==patN   % at loop end, normalize prior belief vectors and 
                              % add 'noise floor' :
                   v.oInfdC(:,cx) = v.oInfdC(:,cx)/sum(v.oInfdC(:,cx));
                   v.oInfdC(:,cx) = v.oInfdC(:,cx) + D.minOthC{blN};
                   v.oInfdC(:,cx) = v.oInfdC(:,cx)/sum(v.oInfdC(:,cx));
                end
             end  
             D.OthInfdC{blN}(:,cx,1) = v.oInfdC(:,cx); % save for posterity first time around ;-)
          end
        else % If there has been previous blocks, apply simple learning rate over 
             % partners to form the belief distribution over all the possible Other's C maps:
            D.OthInfdC{blN}(:,:,1) = (1-ps.blockLR)*D.OthInfdC{blN-1}(:,:,1) + ...
                                        ps.blockLR *D.OthInfdC{blN-1}(:,:,end);
            v.oInfdC(:,:) = D.OthInfdC{blN}(:,:,1);  
        end
   else     % save for posterity v.oInfdC already recorded, which we will change below :
      D.OthInfdC{blN}(:,v.settingLev,v.trial) = v.oInfdC(:,v.settingLev); 
   end
   % Now find max for actual context we are facing
   cx=v.settingLev; % for convenience
   ci =  v.oInfdC(:,cx)==max(v.oInfdC(:,cx));
   v.mapOthC = D.CbOth(:,:,ci);
   % Occasionally the above will not bring up a single max,
   % in which case just average over the equal-height peaks:
   if sum(ci)>1.5
       v.mapOthC = mean(v.mapOthC,3);
   end
   % Store as currmapOthC to do Dkls etc later:
   v.currmapOthC{cx} = v.mapOthC;
   
   if P.detailed  % calc. susprise that v.oInfdC(:,cx) materialised, assuming
                  % equal chances of all the contexts appearing.
       priInfdC = sum(v.oInfdC')/P.settingl;      %#ok<UDIM> % Note deliberate row form.
       D.surp{blN}(D.surpSt{blN}(v.trial),19) = rowpDkl(v.oInfdC(:,cx)',priInfdC);
   end

else  % calc. v.mapOthC based on previous round moves. To do
      % this we first need to find what the Other's policy othPk
   % IF their type was ci and the Self policy was Spolb.
   % This will allow us to calc the likelihood p(D|othPk) =
   % p(D|ci,Spolb).
   cx = v.settingLev;   % for convenience - not used as loop index here.
   othPk = nan(patN,1);   % To hold the other's action that optimizes divergence
                          % from their goals for each possible C enumerated by patN
   othls = nan(P.Nl,1);   % Other's likelihoods over data purely due to the 
                          % emission noise ... 
   if ~P.detailed % then infer on all the Other's observed actions together,
                  % not incrementally for each one observed.
     for pk=1:P.Nl      % pre-calc these as only few. Sadly
                       % no benefit to do in log space here.
        othls(pk) = prod(P.poldens(pk,v.moves(2,:,v.settingLev)));
     end
     % normalize to prevent underflows, though likelihoods don't add
     % up to 1 - any multiplicative factor would do here I think ...
     othls = othls / sum(othls);
     for ci=1:patN
        goal  = D.CbOth(:,:,ci); 
        % Assume that Other would choose correctly, i.e. would choose
        % policy 'who to be' so as to minimise the divergence bet. 
        % 'who WE would be' and goals:
        othPk(ci) = 1;        % initialize w. first policy level
        goalDiv = inf;
        for pk = 1:P.Nl   % a tiny search !
           wedbe = P.poldens(pk,:)' * Spolb; 
           gd = sum(sum( wedbe .* (log(wedbe) - log(goal))));
           if gd < goalDiv; goalDiv=gd; othPk(ci)=pk; end
        end

        % badly un-normalized ;-) posterior
        v.oInfdC(ci,cx)=v.oInfdC(ci,cx) * othls(othPk(ci));
        % if we've got all the un-normalized ones, normalize:
        if ci==patN
           v.oInfdC(:,cx)=v.oInfdC(:,cx)/sum(v.oInfdC(:,cx));
           % add 'noise floor', then normalize again:
           v.oInfdC(:,cx) = v.oInfdC(:,cx) + D.minOthC{blN};
           v.oInfdC(:,cx) = v.oInfdC(:,cx)/sum(v.oInfdC(:,cx));
        end
     end

     D.OthInfdC{blN}(:,v.settingLev,v.trial) = v.oInfdC(:,cx); % save for posterity ;-)
     % At last find max for actual context we are facing:
     ci =  v.oInfdC(:,cx)==max(v.oInfdC(:,cx));
     v.mapOthC = D.CbOth(:,:,ci);
     % Occasionally the above will not bring up a single max,
     % in which case just average over the equal-height peaks:
     if sum(ci)>1.5
         v.mapOthC = mean(v.mapOthC,3);
     end
  
     % no need to bother with  v.currmapOthC{cx} = v.mapOthC;  % refresh after each observation
     % as we do below, in the detailed case.
     
   else % P.detailed is true, do it observatation by observation
                    
     croInfdC = nan(patN,P.nCr); % crossings-other-Inferred-C 
       % hold the belief about other's C at the beginning of each interaction (legacy name = crosing). 
       % It is temporary to this trial, so does not need to remember all the contexts.
     for cr = 1:P.nCr  
       % copy v.oInfdC into temp storage to calc. the shift later
       croInfdC(:,cr) = v.oInfdC(:,cx);  % first time around this will be
       % from last trial where this was inferred. Later, from the last interaction (crossing).
       
       % Now pre-calc other's action likelihood (w.r.t. emission noise). Sadly
       % no benefit to do in log space here.
       othls = P.poldens(:,v.moves(2,cr,v.settingLev));
       % normalize to prevent underflows, though likelihoods don't add
       % up to 1 - any multiplicative factor would do here I think ...
       othls = othls / sum(othls);
       for ci=1:patN
          goal  = D.CbOth(:,:,ci); 
          % Assume that Other would choose correctly, i.e. would choose
          % policy 'who to be' so as to minimise the divergence bet.
          % 'who WE would be' and goals:
          othPk(ci) = 1;        % initialize w. first policy level
          goalDiv = inf;
          for pk = 1:P.Nl   % a tiny search !
             wedbe = P.poldens(pk,:)' * Spolb; 
             gd = sum(sum( wedbe .* (log(wedbe) - log(goal))));
             if gd < goalDiv; goalDiv=gd; othPk(ci)=pk; end;
          end

          % update badly un-normalized ;-) posterior
          v.oInfdC(ci,cx)=v.oInfdC(ci,cx) * othls(othPk(ci));
          % if we've got all the un-normalized ones, normalize:
          if ci==patN
             v.oInfdC(:,cx)=v.oInfdC(:,cx)/sum(v.oInfdC(:,cx));
             % add 'noise floor', then normalize again:
             v.oInfdC(:,cx) = v.oInfdC(:,cx) + D.minOthC{blN};
             v.oInfdC(:,cx) = v.oInfdC(:,cx)/sum(v.oInfdC(:,cx));
          end
       end     % end loop over different poss C patterns.

       % If blow by blow, update D.OthInfdC{blN}(:,v.settingLev,v.trial) only after the last interaction,
       % NOT HERE ! 
       % Find max for actual context we are facing:
       ci =  v.oInfdC(:,cx)==max(v.oInfdC(:,cx));
       v.mapOthC = D.CbOth(:,:,ci);
       % Occasionally the above will not bring up a single max,
       % in which case just average over the equal-height peaks:
       if sum(ci)>1.5
         v.mapOthC = mean(v.mapOthC,3);
       end
    
       %% Calc. 'surprise' / update about other-C :
       %  Dkl simplest based on MAP value only, in col. 13,
       %  this is unsigned-PE like in the sense that it is about the shift
       %  in modal belief, rather than how different the prior and posterior
       %  belief distros were, like oInfdC. See also D.surp_hd .
       D.surp{blN}( D.surpSt{blN}(lastSuchTr)+cr-1,13) = rowpDkl(v.mapOthC(:)',v.currmapOthC{cx}(:)');
       %  Dkl over whole distro over poss other-C, poss in col 14 - see D.surp_hd :
       D.surp{blN}( D.surpSt{blN}(lastSuchTr)+cr-1,14) = rowpDkl(v.oInfdC(:,cx),croInfdC(:,cr));  % use croInfdC, v.oInfdC
       % (may store the entire currmapOthC WHICH FORMED THE PRIOR EXPECTATION about mapOthC, 
       % not the updated version, to use to calc. crossing by crossing expected actions of Other: )
       % D.surp{blN}( D.surpSt{blN}(lastSuchTr)+cr-1,19:(18+P.Nl*P.Nl)) = v.currmapOthC{cx}(:); % retrieve
              % with reshape(D.surp( D.surpSt{blN}(lastSuchTr)+cr-1,19:(18+P.Nl*P.Nl)),P.Nl,P.Nl)
       
       %  Dkl over actions of other, based on expectations about other's 
       %  action before update. Also 'PEothact','slnPothact' simply based on 
       %  observed action and previously stored 
     
       v.currmapOthC{cx} = v.mapOthC;  % refresh mapOthC *after* each observation.
                                       % but croInfdC at the beginning of each loop over crossings.
       
     end    % loop over all the crossings in the prev. trial of this type
     
     % Now store the inital beliefs for this trial, before thinking about policies ...
     D.OthInfdC{blN}(:,v.settingLev,v.trial) = v.oInfdC(:,cx); % save for posterity ;-)
           
   end  % if statement re. whether to infer on all contributions together or blow-by-blow
end % if there is a previous exchange to update one's beliefs with.
  
% Now we can fill in expectations (post. of previous trials) from the 
% lastSuchTr+1 to current (v.trial) inclusive, and can calc. suprises for 
% lastSuchTr+1 to lastNotSuchTr :
if exist('lastSuchTr')
  if v.trial > lastSuchTr+1   % as lastSuchTr is done, and so is v.trial
   D.OthInfdC{blN}(:,v.settingLev,(lastSuchTr+1):(v.trial-1)) = repmat(D.OthInfdC{blN}(:,v.settingLev,v.trial),1,v.trial-lastSuchTr-1);
  end
  % now calc. initial suprises for trials of actual interaction
  if (lastNotSuchTr >= lastSuchTr+1) 
    lastTr = min([lastSuchTr+1, P.trN(blN)]);
    for trn = lastTr:lastNotSuchTr
       priOthInfdC = mean(D.OthInfdC{blN}(:,:,trn)'); %#ok<*UDIM>
       hLev = P.setting{blN}(trn);
       D.surp{blN}(D.surpSt{blN}(trn),19) = rowpDkl( D.OthInfdC{blN}(:,hLev,trn)', priOthInfdC );
    end
  end
end
%% Section 2c Calculate the expectation for how the Other will play (or
%  would play for the very last, dummy trial)
%  ... on the basis of which  Self will choose own policy to satisfy own C matrix later.
%  cf. 1.4.3. in WhoToBeWithYou_2 ... determines policy that minimizes
%  divergence from their own goals.

% . . . . First Calculate the expectation for how the Other will play . . . .
% This accummulates the actions that would transpire for each possible C
% matrix of Other, weighed by the belief that this C obtains.

% May very well do the following by vectorizing the 'wedbe' and 'goal'
% matrices and then performing both the Dkl and the finding of max.
% via matlab/octave sum and max operations, that naturally work on
% matrix columns. However here do it all explicitly for clarity.
% Onwards to accumulate policy expected to be performed by the pt!
v.othExPol = zeros(1,P.Nl);          % policy expected by the pt overall.
Spers =  v.whoSis(v.settingLev,:);     % Self's latest persona!
Spolb =  sum(repmat(Spers(:),1,P.Nl) .* P.poldens); % Self-policy-beliefs ie choice beliefs about Self
wedbe = nan(P.Nl,P.Nl,P.Nl);
for pk = 1:P.Nl; wedbe(:,:,pk)= P.poldens(pk,:)' * Spolb; end
  
% For 'blow by blow' I need to do the for loop below for each
% incremental value of in croInfdC ...
if ~P.detailed ||  isnan(v.moves(1,1,v.settingLev))
  for ci = 1:patN  % very similar loop to above ...
    goal  = D.CbOth(:,:,ci);
    % find Other-Agent's best 'who to be' policy so as to minimise the
    % divergence bet. 'who WE would be' and goals:
    othPk(ci) = 1;        % initialize w. first policy level
    goalDiv = inf;
    for pk = 1:P.Nl   % a tiny search !
        we = wedbe(:,:,pk);
        gd = sum(sum( we .* (log(we) - log(goal))));
        if gd < goalDiv; goalDiv=gd; othPk(ci)=pk; end
    end
    % accumulate marginal :
    v.othExPol = v.othExPol + P.poldens(othPk(ci),:) * v.oInfdC(ci,v.settingLev);
  end
else  % Interaction - by - interaction
  for cr = 1:P.nCr
    v.othExPol = zeros(1,P.Nl);   % reset for each crossing - experience should
    % be encoded in croInfdC(:,cr) as cr grows.
    for ci = 1:patN  % very similar loop to above ...
      goal  = D.CbOth(:,:,ci);
      % find Other-Agent's best 'who to be' policy so as to minimise the
      % divergence bet. 'who WE would be' and goals:
      othPk(ci) = 1;        % initialize w. first policy level
      goalDiv = inf;
      for pk = 1:P.Nl   % a tiny search !
          we = wedbe(:,:,pk);
          gd = sum(sum( we .* (log(we) - log(goal))));
          if gd < goalDiv; goalDiv=gd; othPk(ci)=pk; end
      end
      % Accummulate marginal summing over ci:
        v.othExPol = v.othExPol + P.poldens(othPk(ci),:) * croInfdC(ci,cr);
    end  % loop over poss. other C
      
    % Now store expected action of Other and also the surpise that
    % accrued once the relevant actions were observed:
    D.surp{blN}( D.surpSt{blN}(lastSuchTr)+cr-1,15:18) = v.othExPol;
    % REM D.surp 'PEothact','lnPothact','DklmapOthC','DklothC',...      % 11 to 14, surpr about other.
    % so, first the PE, then the lnPothact :
    D.surp{blN}( D.surpSt{blN}(lastSuchTr)+cr-1,11) = dot(1:4,v.othExPol) - v.moves(2,cr,v.settingLev);
    D.surp{blN}( D.surpSt{blN}(lastSuchTr)+cr-1,12) = log(v.othExPol( v.moves(2,cr,v.settingLev)));
      
  end % loop over crossings
    
end  % if detailed (crossing by crossing) or not

%% from here one, we only process trials with actual decisions 
%  not the last P.settingl ones where only 'looking back' surprises etc are needed.
if v.trial <= P.trN(blN)

  %% Section 3 - choose Self policy to min div from Self goal . . . .
  % As the Self's C matrix has NOT been formed, unlike d.OC(:,:,1 to P.settingl)
  % (which was in getIC1partner), if in the first round, construct it.
  % Also find the 'ideal choice' that Other can make to allow the pt. to
  % have themselves an 'ideal choice' available
  % so as to approch their Self C-matrix as closely as poss.
  % NB Other is cols, Self is rows.
  if v.trial ==1
    % Construct Cself and store it in D.SelfC for each of the two contexts
    v.uSC = P.polU;  % uncertainty modulation to be used in making Self C
    % matrix. Keep it fairly sharp, with the objective of making Self able to
    % be responsive but not spasmodic, with tuning of choice stochasticity
    % left to Temp/precision parameter ... P.polU should be OK, numerically ...
    D.SelfC{blN} = zeros([P.Nl,P.Nl,P.settingLevN]);
    for i=1:P.Nl  % NB   prior pref   * conditional   :
        for k=1:P.settingLevN
            D.SelfC{blN}(:,i,k) = ps.SPartnPr(i) * noisyBino((ps.Spref(k,i)+1)/(P.Nl+1),v.uSC,P.Nl);
        end 
    end
    
    %        ** find 'ideal choices' for Self and Other for each context **
    %
    %      * Construct the row-wise version of the SelfC matrix and store it *
    R = P.Nl; C=P.Nl; N=P.Nl; % hopefully to make role clearer ...
    D.rwiseCs{blN} = zeros([R*C,N*N,P.settingLevN]);  % row-wise lot of SelfCs
    for cx = 1:P.settingLevN
      % a vectorized version of repetitions of SelfC:
      bigCself = repmat(D.SelfC{blN}(:,:,cx),[R C]);  bigCself = bigCself(:);
      colbigCs = zeros(size(bigCself));   % make space ...
      colbigCs(P.lindPC) = bigCself;      % put the elements in the right place ...
      D.rwiseCs{blN}(:,:,cx)  = reshape(colbigCs,[N*N, R*C])'; % Here the last bit in
      % square brackets works even if it's the wrong way round - I *think* its the right way round ...
    end
    %   Now can look for 'ideal choices' :
    D.ideals{blN} = zeros([2, P.settingLevN]);       % rows for self=1, other=2
    for cx = 1:P.settingLevN
      rdkl = rowpDkl(P.rowPolComb, D.rwiseCs{blN}(:,:,cx)) ;
      minlind =  find(rdkl == min(rdkl),1,'first');
      [D.ideals{blN}(1,cx), D.ideals{blN}(2,cx)] = ind2sub([R,C],minlind);
    end
    
  end
  % Now can use this SelfC matrix, but gd now will record Dkl for each action :
  goal = D.SelfC{blN}(:,:,v.settingLev);
  % Speak = 1;
  % goalDiv = inf;
  gd  = nan(1,P.Nl);
  for i = 1:P.Nl
    wedbe = P.poldens(i,:)' * v.othExPol ;
    gd(i) = sum(sum( wedbe .* (log(wedbe) - log(goal))));
    % if gd(i) < goalDiv; goalDiv=gd(i); Speak=i; end;
  end
  
  
  %%                        Key outputs
  
  v.SpolPr = psoft(-gd,1/ps.T);   % Provisional simple softmax not policy precision yet ...
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if P.synth(1)  % Synthetic data: v.nmoves for self will be updated within
      % here, within this fn, but for Other it should be done (just) after this fn, in the parent fn.
    if P.maxPsynth==0  % here we haven't been told to use only the most likely policy
      v.Spol = pBinSample(v.SpolPr,1);   % just one policy for all rounds in trial ...
      % v.Spol should go from 1 on, not from 0 on ...
      v.nmoves(1,1:v.crN,v.settingLev) = pBinSample(P.poldens(v.Spol,:),v.crN); % ... gives
      % rise probabilistically to the work actions in the different trials .
    elseif P.maxPsynth==1  % if we've been told to use only the most likely policy :
      v.Spol = find(v.SpolPr == max(v.SpolPr));
      v.nmoves(1,1:v.crN,v.settingLev) = ones(1,v.crN)* v.Spol;
    else
      error(['Not ready for P.maxPsynth=', num2str(P.maxPsynth)]);
    end
    % Store selected policy:
    D.evo{blN}(P.evoPols{blN}(v.trial),4) = v.Spol;
    % Below, and in HT versions, no need to store the subset of data previously called key_data. 
    % Everything to be retrieved from d.evo .
    % play out crossings of this trial:
    [v, D] = playSimHT1game( P, D, v );
    v.synth(v.trial,:)  = [nan, v.trial,  v.settingLev, v.Spol, ...
      v.nmoves(1,1:v.crN,v.settingLev), ...
      v.nmoves(2,1:v.crN,v.settingLev), gd ];
    % Record the actions in D.evo -- col 8 is for sWork, 21 is for oWork : 
    D.evo{blN}(P.evoSt{blN}(v.trial):(P.evoSt{blN}(v.trial)+v.crN-1),8)  = v.nmoves(1,1:v.crN,v.settingLev);
    D.evo{blN}(P.evoSt{blN}(v.trial):(P.evoSt{blN}(v.trial)+v.crN-1),21) = v.nmoves(2,1:v.crN,v.settingLev);
  else     % if no synthetic data, look up v.Spol in D.evo:
    v.Spol = D.evo{blN}(P.evoPols{blN}(v.trial),4);
  end
  
end  % processing for full trials with decisions, not the dummies at the end.
%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if v.trial <= P.trN(blN)
   ll = log(v.SpolPr( v.Spol));  % v.nmoves irrelevant, as we know v.Spol !
   if P.synth(1); v.synth(v.trial,1)=ll; end % fill in place-held log-lik values
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   if P.detailed  && ~P.synth(1) % fuller record incl. and details of 
     % real data ll, actions etc.  
     v.llreal(v.trial,:) = [ll, v.trial, v.settingLev, v.Spol, ...
                         v.nmoves(1,:,v.settingLev),...
                         v.nmoves(2,:,v.settingLev), gd ];
   end
else  % for dummy trials at the end
   ll = 0; % we are sure that there will actually be no v.Spol ;-)
end
%%
return;

