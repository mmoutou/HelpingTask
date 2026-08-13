function [v, dsynth] = playSimHT1game( p, dsynth, v )
%  [v, dsynth] = playSimIC1game( p, dsynth, v )
%  Play out set of SIMULATED interactions. Needs self and other actions
%  to be already provided in v.nmoves, unlike playOutIC1game on which it's based.
%  Return the returns depending on the type of captain. Also calculate the contribution
%  to fees based on the returns here, so needs p.trWinCoeff

blN = v.block; 
cx = v.settingLev;

for crn = 1:v.crN
   v.rets(1,crn,cx) = p.ret(v.nmoves(1,crn,cx), v.nmoves(2,crn,cx), 1, cx);
   v.rets(2,crn,cx) = p.ret(v.nmoves(1,crn,cx), v.nmoves(2,crn,cx), 2, cx);
end  % loop over crossings

% Winnings:
v.Sfee = sum(v.rets(1,:,v.settingLev))*p.trWinCoeff;
v.Ofee = sum(v.rets(2,:,v.settingLev))*p.trWinCoeff;
 
%%  store indiv. returns 
%
%try   % if v.practIC has been defined, leave it alone; otherwise make it
%      % and set it to zero.
%  v.practIC;
%catch
%  v.practIC = 0; 
%end
%if v.practIC == 0   % i.e. otherwise, if this IS a practice run, don't store stuff.
  dsynth.evo{blN}(p.evoSt{blN}(v.trial):(p.evoSt{blN}(v.trial)+p.nCr-1),9:10) = v.rets(:,:,cx)';
  % In HT1, extra copy below not needed. Had:
  % dsynth.key_data(p.keySt(v.trial):(p.keySt(v.trial)+p.rN(v.trial)-1),5:6) = v.rets(:,:,cx)';
  % ... and fees:
  dsynth.evo{blN}(p.QRow{blN}(v.trial),[22, 23]) = [v.Sfee, v.Ofee];
  % again below commented out as per above.
  % dsynth.key_data(p.keySt(v.trial), [10, 11])= [v.Sfee, v.Ofee];
%end

return; % end of function playOutIC1game


