function psPr0 = HT1psPr0Set( psPr0Code)
% PSPR0=  function to provide priors for t-fitting (neuroecon)
%         psPr0Code = 0 (default) for weak priors (uninformative for most)
%         add to the list below as needed.

try
    psPr0Code;
catch
    psPr0Code = 0;
end

Nl = 4; setLevN = 1;  % In HT1, have 4 levels of choice for effort and only 1 'context' (per Other)
psPr0.Spref0 = nan(Nl,Nl,setLevN); % Spref rows, Owrk cols, context pages, flat. 

if psPr0Code == 0   % As uninformative as practical
    
    for o=1:Nl
        for c=1:setLevN
            m = o/(1+Nl); u = 1.5*m*(1-m); 
            psPr0.Spref0(:,o,c) = discBetaMU(m,u,4); 
        end
    end
    psPr0.prevp0 = [1.05, 1.05];  % A and B for betalike for pSucc of noisyBino describing prevPri
    psPr0.prevu0 = [1.01, 1];   % [1.01, 1] rises very rapidly from 0 and then falls slowly.
    
    psPr0.SPartp0 =  [1.05, 1.05]; % [45, 10];  for highly constrained diagnostic. Uninf. is: psPr0.SPartp0 = [1.05, 1.05]; 
    psPr0.SPartu0 = [1.01 1]; % [2, 0.04] for diagnostic;  % check with x = 0:0.01:10; plot(x,gampdf(x,1.5,0.1)) % for less highly constrained.
    
    psPr0.T0 = [1.5, 0.5];        % A and B for gamlike on T
    psPr0.blockLR0 = [1.1, 1.4];  % A and B for betalike for betaLR. Modestly discrourages very high apparent LRs.

elseif psPr0Code == 1  % Setting SPartp0 SPartu0 to medians of Aug03 informative fit

    for o=1:Nl
        for c=1:setLevN
            m = o/(1+Nl); u = 1.5*m*(1-m); 
            psPr0.Spref0(:,o,c) = discBetaMU(m,u,4); 
        end
    end
    psPr0.prevp0 = [1.05, 1.05];  % A and B for betalike for pSucc of noisyBino describing prevPri
    psPr0.prevu0 = [1.01, 1];   % [1.01, 1] rises very rapidly from 0 and then falls slowly.
    
    SPartp0M = 0.7280;  conc=1000; 
    psPr0.SPartp0 =  [SPartp0M*conc (1-SPartp0M)*conc ];   
    % check with  x = 0:0.001:1; y=betapdf(x,psPr0.SPartp0(1), psPr0.SPartp0(2)); plot(x,y); sum(x .* y)*0.001

    SPartu0M = 0.44147;  conc=100; 
    psPr0.SPartu0 =  [SPartu0M*conc 1/conc ]; 
    % check with  x = 0:0.01:100; y=gampdf(x,psPr0.SPartu0(1),psPr0.SPartu0(2)); plot(x(1:1000),y(1:1000)); sum(x .* y)*0.01

    psPr0.T0 = [1.5, 0.5];        % A and B for gamlike on T
    psPr0.blockLR0 = [1.1, 1.4];  % A and B for betalike for betaLR. Modestly discrourages very high apparent LRs.


else
    error(['This prPr0Code, ' num2str(psPr0Code) ' is not catered for ']);
end

return;  %% ~~~~~~~~~~~~~~~~~~~ eof ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
