function [priPar, prHd] = priParEvalRL(col2keep, col2ch,colVal)
%PRIPAR8BL provide array of priors on parametes in native space.
%   This should work with dens = dbetasc(x,a,b,lo,hi) as used in e.g. sel8bl03
%   Usually all arguments can be kept blank, to use defaults. 
%   Defaults are weak (regularizing) except for lapse rate, where there is
%       an informative prior, so that it is quite rare for it to be >0.5
%   col2keep columns can be selected from the full array, if less than all are needed.
%   col2ch is the columns to be replaced from their default values
%   colVal is an 4 x col2ch array with the new values. 
%
%   priPar is the output array, prHd its column labels/headings. 

% Maximal repertoire - usually many of these will be coalesced/dropped/fixed:
%          1       2     3     4      5       6       7      8        9       10    11    12
prHd = {'lambda','eta','wEx','sig','wOS11','wOS12','wOS21','wOS22','lambda2','lps','E0s','E0o'};
        
try
    col2keep;
catch
    col2keep = 1:length(prHd);
end
try 
    col2ch;    % Columns to change from default to values given by colVal.
catch
    col2ch = [];
    colVal = [];
end


%         1       2     3     4       5      6      7        8        9      10      11    12    13 
%prHd={'lambda','eta','wEx','sig','wOS11','wOS12','wOS21','wOS22','lambda2','lps', 'E0s', 'E0o','EvBlockLR' };
%                        < max at 1> <------ SD ~10----------->                                          <max at 4> 
priP= [[1.01,  1.01,  1.01, 1.01,    10,    10,     10,     10,    1.01,    1.01,   10,    10,   1.01];... %  10,    1.2  ]; ...  % a
       [1.01,  1.01,  1.01,  2,      10,    10,     10,     10,    1.01,    3.01,   10,    10,   1.01];... %  10,    5.8  ]; ...  % b
       [ 0,     0,     0,    0,     -46,   -46,    -46,    -46,      0,      0,     -46   -46,    0  ];... % -46,     0   ]; ...  % lo
       [ 1,     1,     1,   100,     46,    46,     46,     46,      1,      1,      46,   46     1] ];   %  46     100  ]];     % hi

if ~isempty(col2ch)
    priP(:,col2ch) = colVal;
end

priPar = priP(:,col2keep);
prHd = prHd(col2keep);

return;  % end function priParEvalRL


