function msl = HT1lp03( restp, pS, d, p, psPr, details)
%HT1LP03 wrapper of HT1ll03 to fit the firs't', or neuroecon, part, i.e. 
%   prevPri, SPartnPri, T, and blockLR params with fmincon etc.
%   Derived from IC2lp2 ... fitting without emotional evaluations (see HT*lp*f.m for fitting emot. eval.)
%   returns mainly the minus-sum-log-posterior (prior being psPr over all params)
%   but also the sum log likelihood and a copy of the parameters as actually
%   tried (as opposed to the separate input baseline params pS and 
%   the components to of interest here, prevp, prevu, Spartp, Spartu, T, blockLR.

try
  details;
catch
  details = 0;
end
%% consturct new participant-parameter structure, keeping the 
%  crucial ps.Spref well alone :
ps = pS;
ps.prevp = restp(1);   % REM this is prevalence prior about others' type (biskBeta M = pSucc)
ps.prevu = restp(2);   %  ... spread (uncertainty) param. for the above.
ps.SPartp = restp(3);  % Spartp is Self preference peak about Other's type, Active Inf C map style.
ps.SPartu = restp(4);  %  ... and its spread / uncertainty
ps.T = restp(5);
ps.blockLR=restp(6);
ps.prevPri  = discBetaMU(ps.prevp,ps.prevu,p.Nl);    % prevalence prior distro.
ps.SPartnPr = discBetaMU(ps.SPartp,ps.SPartu,p.Nl);  % Self pref about others distro.

%% Only now, and if appropriate, calc log prior of params:
if isempty(psPr)
  lnPri = 0;
else
  %   In HT1pslPr03, the priors over the u parameters are over abs(u), as we want to avoid
  %   values very near 0 but to consider both positive and negative ones otherwise.
  lnPri = HT1pslPr03(ps, psPr, p);
end

%% outputs: sum log posterior 
%  and other details if asked for ...
if details
  msl.lik  = HT1ll03(ps, d, p);
  msl.mpost = -(msl.lik + lnPri);
  msl.ps = ps;
else
  msl = -(HT1ll03(ps, d, p) + lnPri);
end

return;  % end of function HT1lp03

