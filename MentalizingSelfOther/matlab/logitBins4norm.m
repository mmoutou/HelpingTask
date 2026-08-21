function [pmf, normBounds] = logitBins4norm(mu,sig,unitIntBounds)
%LOGITBINS4NORM gives the mass density function created by binning a Gaussian. 
%
% [pmf, normBounds] = logitBins4norm(mu,sig,unitIntBounds)
% 
% Each bin comes from 'slicing' a Gaussian with mean mu and SD sig, sliced
%    into bins. However, we map the domain of the Gaussian onto the unit interval,
%    *** and we define the boundaries of our bins in the unit interval via unitBounds ***
%    If unitBounds is not provided, we take it to be [0,1/6,...5/6,1] by default.
%  pmf is a resulting Probability Mass Function. 
%  normBounds are the boundaries from unitIntBounds mapped to the real line  
%    via a logistic sigmoid. 

%% First, provide as default our 6 bins as per our equations:
try
    unitIntBounds;
catch
    unitIntBounds = [0, (1:5)/6, 1];
end

%% Defensive coding to catch various possible input errors:
if ~isvector(unitIntBounds)
    error('please provide unitIntBounds as a vector of at least 2 elements');
end
if sig < 0 
    error('sig must be >= 0 - it''ll be the SD as per normcdf(x,mu,sig)');
end
epsi = 1e-15;  % an 'almost zero' number, to mop up rounding errors
if unitIntBounds(1) < 0
    if unitIntBounds(1) < -epsi
        error('unitIntBounds(1) must not be negative');
    else  % ... there was probably a rounding error!
        unitIntBounds(1) = 0; 
    end
end
if unitIntBounds(end) > 1
    if unitIntBounds(end) > 1 + epsi
        error('unitIntBounds(end) must not be > 1');
    else % ... probable rounding error, so:
        unitIntBounds(end) = 1;
    end
end

%% Main calculations

% total number of bins is one less than the number of bounds:
totBoundN = length(unitIntBounds);
totBinN = totBoundN - 1;

% normBounds is formed by mapping unitIntBounds via the logit function:
normBounds = logit(unitIntBounds);

% We calculate cummulative densities in the following vector. 
% It might be quicker to set the ends to 0 and 1 by hand if the 
%    unitIntBounds goes from 0 to 1, but I leave this to you ;)
cdf = normcdf(normBounds,mu,sig);
pmf = cdf(2:end) - cdf(1:(end-1));

% emit a warning if the pmf doesn't sum up to 1 :
if abs( sum(pmf) - 1) > 1e-6
    warning('pmf didn''t add up to 1 - we normalize and proceed anyway');
end
pmf = pmf/sum(pmf);   
    
return;

