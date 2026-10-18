function [slPri, lpri] = HT1pslPr03f(ps, psPr)
%HT1PSPRI03F return minus sum log priors and/or structure with the feelings component prior probs of ps 
%   only, under the prior psPr. Assume independent priors over these components of ps, as per e.g.:
% psPr :
%  psPr.feelm0 = [1.025, 10.0];    % A and B for gampdf for m of rat2resp - flat max near 1
%  psPr.feelu0 = [1.025, 10.0];    % A and B for gampdf for u of rat2resp - flat max near 1

if isfield(ps,'feelm') % if we are fitting model with feelm0, feelu0
    % for feelm :
    lpri.feelm = - gamlike(psPr.feelm0, ps.feelm );
    % for feelu :
    lpri.feelu = - gamlike(psPr.feelu0, ps.feelu );
    % sum :
    slPri = lpri.feelm + lpri.feelu;
elseif isfield(ps,'evalRL')  % if we are fitting model with evalRL
% Maximal repertoire - usually many of these will be coalesced/dropped/fixed:
%          1       2     3     4      5       6       7      8        9       10   11    12    13
%prHd={'lambda','eta','wEx','sig','wOS11','wOS12','wOS21','wOS22','lambda2','lps','E0s','E0o','EvBlockLR'};
    natP = [invlogit(ps.evalRL(1:3)),exp(ps.evalRL(4)),ps.evalRL(5:8),invlogit(ps.evalRL(9:10)),...
            ps.evalRL(11:12)];  % DO NOT WRITE :end) HERE - code defensively in case ps.evalRL is 
                                % the wrong length etc.
    % evalRL(13), EvBlockLR, is a learning rate on the unit interval like the others,
    % so give it a prior whenever the parameter and a prior column for it both exist.
    % priParEvalRL returns only 12 columns (its prHd lists 12 names), so use
    % priParEvalRLHT1, which returns 13, if you want parameter 13 regularised.
    if length(ps.evalRL) >= 13 && size(psPr.evalRL0,2) >= 13
        natP = [natP, invlogit(ps.evalRL(13))];
    end
    nP = length(natP);
    lpri.evalRL = log( dbetasc(natP,psPr.evalRL0(1,1:nP),psPr.evalRL0(2,1:nP),...
                                    psPr.evalRL0(3,1:nP),psPr.evalRL0(4,1:nP)));
    slPri = sum(lpri.evalRL) ;                           
else
    error('This form of ps is missing all fields catered for so far.')
end
    
return;

