

function p = peaks4c (p)
% function p = peaks4c (p)
%  populating 'continuous' combinations of preference peaks in
%  the conditional part of a C matrix pref(me|you) and a grand
%  marginal pref(you-I-want-to-deal-with). The last is in the
%  last column of p.peaks4c that is returned. Only 'continuous'
%  cases of the conditional part are allowed, i.e. no slopes
%  more than +/- 1. It contains indices (1 to Nl), not values.
%
% Needs to receive p.Nl within p.
%
% Maybe can consider uncertainties re. taste in you too but 
% take e.g. noisyBino(yPref/(N+1),N/3,N) as reasonable default ...

try
  p.Nl;
  if p.Nl > 7
    error('p.Nl too big for exhaustive enumeration'); 
  end
catch
  p.Nl = 4;
end

Nrow=p.Nl^p.Nl;
a = 1:p.Nl;
pkar = nan(Nrow,p.Nl);
for cl = 1:p.Nl;
    nrow = p.Nl^(p.Nl-cl+1);     % no. of rows per block in this column
    blk = repmat(a,Nrow/(p.Nl^cl),1);   blk=blk(:);
    pkar(:,cl)=repmat(blk,Nrow/length(blk),1);
end
% now remove unwanted rows that contain too big a jump
for ro = 1:Nrow
    for cl = 2:p.Nl
        if abs(pkar(ro,cl-1) - pkar(ro,cl)) > 1.5
           pkar(ro,:) = nan;
           break;
        end
    end
end

pkar=pkar( ~isnan(pkar(:,1)),:);

%% finally replicate for each level of peak taste in you 
p.peaks4c = nan(size(pkar,1)*p.Nl,1+p.Nl);
p.peaks4c(:,end) = reshape(repmat(a,size(pkar,1),1),p.Nl*size(pkar,1),1);
p.peaks4c(:,1:p.Nl) = repmat(pkar,p.Nl,1);

return;
