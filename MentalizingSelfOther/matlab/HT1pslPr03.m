function [slPri, lpri] = HT1pslPr03(ps, psPr, p)
%PSLPRHT1B return minus sum log priors and/or structure with component prior probs of ps under the prior psPr
%  Developped from pslPrIC2 . 
%  In this version, the priors over the u parameters are over abs(u), as we want to avoid values 
%      very near 0 but to consider both positive and negative ones otherwise.
%   Assumes independent priors over the components of ps, whose sufficient stats are
%      denoted by Spref0, prevp0 etc., as per e.g.:
%  psPr.Spref0 = ones(p.Nl,p.Nl,2)/p.Nl; % Spref rows, Owrk cols, setting (context) pages, flat. 
%  psPr.prevp0 = [1.01, 1.01];  % A and B for betapdf for pSucc of noisyBino describing prevPri
%  psPr.prevu0 = [1.5, 1];    % mu and sigma for 1/U of noisyBino describing prevPri
%  psPr.SPartp0 = [1.05, 1.05]; % This and next for SpartnPr
%  psPr.SPartu0 = [1.5,];   % mu and sigma for 1/U for SpartnPr
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

% OZ Aug 26: lpri.prevp and lpri.prevu are computed above but were not in this sum,
% so prevp and prevu were the only two parameters fitted without a prior while their
% counterparts SPartp and SPartu had theirs. The consequence is visible in the fits:
% prevu drifting towards zero raises the binomial to a huge power and collapses the
% prevalence belief onto a single type, e.g. prevu = -0.043 giving
% prevPri = [1, 1e-17, 2e-23, 4e-18]. The gamma prior on abs(prevu) penalises exactly
% that corner, 3.19 nats at |prevu| = 0.043 against 1.29 at |prevu| = 1.97, and it
% diverges at 0, so the degenerate solution becomes unreachable rather than merely
% unlikely. Both terms added below.
slPri = lpri.t+lpri.prevp+lpri.prevu+lpri.spartnu+lpri.spartnp+sum(lpri.spref(:))+lpri.blockLR;

return;

