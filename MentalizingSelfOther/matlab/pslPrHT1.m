function [slPri, lpri] = pslPrHT1(ps, psPr, p)
%PSPRIHT1 return minus sum log priors and/or structure with component prior probs of ps under the prior psPr
%  Developped from pslPrIC2 . 
%  In this version, the priors over the u parameters are over abs(u), as we want to avoid values 
%      very near 0 but to consider both positive and negative ones otherwise.
%   Assume independent priors over the components of ps, as per e.g.:
% psPr :
%  psPr.Spref0 = ones(p.Nl,p.Nl,2)/p.Nl; % Spref rows, Owrk cols, setting (context) pages, flat. 
%  psPr.prevp0 = [1.01, 1.01];  % A and B for betapdf for pSucc of noisyBino describing prevPri
%  psPr.prevu0 = [2.0, 2.0];    % A and B for gampdf for U of noisyBino describing prevPri
%  psPr.SPartp0 = [1.05, 1.05]; % This and next for SpartnPr
%  psPr.SPartu0 = [2.0, 2.0];
%  psPr.T0 = [1.5, 0.5];        % A and B for gampdf on T
%  psPr.blockLR0 = [1.01, 1.01];  % A and B for betapdf for blockLR
%
% REM
%    ps.Spref      % self preferences
%    ps.prevPri    % prevalence priors about the Other that pt holds.
%    ps.SPartnPr   % preferences over types of partners

% use linear indexing and strides to access psPr.Spref. This is only
% of much use if psPr.Spref is 3D, which in HT1 it isn't - legacy and 'one day' code here!
str1 = (0:(p.settingLevN*p.Nl - 1))*p.Nl;   % stride derived offset for linear indices for ps.Spref
lispref = ps.Spref' + 1; 
lispref = lispref(:)' + str1;  % turn ps.Spref into linear indices 
% look up and put into original shape : 
lpri.spref = log(reshape(psPr.Spref0(lispref)',[p.Nl,p.settingLevN])'); 

% for prevPri:
lpri.prevp = -betalike(psPr.prevp0, ps.prevp );
lpri.prevu = - gamlike(psPr.prevu0, ...
                       abs(ps.prevu) );
% for SpartnPr :
lpri.spartnp = -betalike(psPr.SPartp0, ps.SPartp );
lpri.spartnu = - gamlike(psPr.SPartu0, ...
                         abs(ps.SPartu) );  % REM nlogL = gamlike(params,data) : negative of log-likelihood of params.

% for T:
lpri.t =  - gamlike( psPr.T0, ps.T);

% for blockLR: 
lpri.blockLR = -betalike(psPr.blockLR0, ps.blockLR );

slPri = lpri.t+lpri.spartnu+lpri.spartnp+sum(lpri.spref(:))+lpri.blockLR;

return;

