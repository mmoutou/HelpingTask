function [pols, cents, U] = actUncert( nl, st )
% [pols, cents, U] = actUncert( nl, st ) provides an nl x nl
%   matrix such that the kth row has peak 
%   policy prob at k and standard dev approximately st, where
%   the policy centers cents are scaled, cents=((1:nl)-0.5)/nl
%   U are the uncertainties used by noisyBino to produce pols
%   Default arguments (nl,st)=(4,0.15)

try, st; catch st=0.15; end;
try, nl; catch nl=4;   end;
if nl<3 || nl > 10
  error('Ensure 3>=nl>=10');
end
if st < 0.01 || st > 0.25
  error('Ensure 0.01<=st<=0.25');
end
N = 100;
du = 1/N;
x = ((1:nl)-0.5)/nl ;     cents = x; 
U = nan(nl, 1);
for l=1:nl
  si=nan(20*N,1);
  for i=1:(20*N)
    u = du*i;
    p=noisyBino(l/(nl+1),u,nl); 
    mu=x*p'; 
    si(i)=sqrt(x.*x*p'-mu^2); 
    if i>1 && si(i)>=st; 
      % disp(['si(i)=' num2str(si(i))]); 
      % disp(['u=',num2str(u)]);
      break; 
    end;
  end
  % linearly interpolate
  cu = - (si(i)-st)*du/(si(i)-si(i-1));
  u = u + cu;
  U(l) = u;
  pols(l,:) = noisyBino(l/(nl+1),u,nl);

end   % end loop over different modal policies

return;

