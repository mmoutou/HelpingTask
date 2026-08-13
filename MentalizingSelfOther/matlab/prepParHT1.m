function p = prepParHT1( p )
% prepare Parameters for Tyler's Helping Task - simple Caring / Uncaring avatar version.
% Add various fields to the parameter structure and 
% return them. paramHT1 also contains provision for 'partners' to come
% from explicitly different 'groups' and has a smaller number of trials.
% p must already have p.code, p.debug, p.HTdir, p.data_dir, p.helperType,
%   and optionally p.commDir, which is meant to be network dir for remote access.

cwd = cd;  
% Blocks, helper types, and trial numbers: these must
% be provided, only example here:
% p.helperType =  [[2, 10]; [1,  13]];  % Block 1 is of helper type 2 (=Caring) and 
%                                       % lasts for 10 trials, etc.
p.trN = p.helperType(:,2);

%%  fMRI scanning related parameters, for future reference
p.port = 0;              % 0 for emulscanning, 1 for scanner in serial port.
p.slicesPerVolume = NaN; % For fMRIT. Was 48;
p.dummySlices = NaN;     % For fMRIT. Was 6;
p.Slices2wait4 = p.dummySlices * p.slicesPerVolume + 1;  % to wait for equilibration
p.emulTotalSlices = NaN; % For fMRIT. Was 100000; % A large default, would cover 1h 40 min at 60 ms per slice
p.emulScanner = 'NONE';  % For fMRIT. Was 'allegra';
p.jitterSpan = NaN;      % For fMRIT. Was 4500; % O responses will vary within an interval of this many ms.

%% Filenames and where to find / write things:
% p.matlab = version;     % there should already be a numeric p.matVer
fs = filesep;
p.fileName = [p.code '_HT1_'  datestr(now,'yyyy-mm-dd_HH_MM')];
p.tUnitsPerSec = 1000; % Some packages measure time in ms, some in sec.
p.invalid = -6666;
p.sp0 = 100;

p.prog_dir = p.HTdir;
p.prtnPreFile = 'NONE'; % [p.code(1:(end-1)),num2str(3-str2num(p.code(end))) 'prefs.mat']; 
p.myPreFile   = [p.code 'prefs.mat'];
% directory where we can write files to comm-unicate from one pt. 
% in a pair to another (for versions where avatars really represent players):
% Was: p.commDir = '\\abba\fMRI_REDIT\RiverCrossing\tempData';
if ~isfield(p,'commDir') || ~exist(p.commDir,'dir')
  p.commDir = [getenv('HOMEDRIVE') getenv('HOMEPATH') fs 'Dropbox' fs 'task_code' fs 'IPD' fs 'comm'];
end

%% Display (was Cogent) and keyboard related:
% p.context = 1;    % 'context' here refers to HARDWARE CONTEXT, NOT THE TASK FACTOR, so 1 
%                   % is a default keyboard / layout, as e.g. for NSPN Dell laptops.
% p.resolution = [1920 1080];  % Typical for 2026 ... was [1024 768];
% % p.display_type = 1; %1 fullscreen, 0 windowed
% p.display_type = 1; % set to 0 for debug if needed below. 
% if p.debug ~=0; p.display_type = 0; end
% p.image_dir = [p.HTdir fs 'images' fs];
                

%% task related:  ---------------------------------------------------------------%
p.Nl = 4;       % levels of possible investment / effort in each round.
p.settingl = 1;   % Can have e.g. 2 task settings, or contexts, encoded by returns matrices.
                  % In the Helping Task, only one. Was p.harshl in Moutoussis, Gosalia et al '25 paper.
p.pol =0:(p.Nl-1);  % allowable policies, i.e. levels of contribution.
p.polsd = 0.15;     % alternative policies standard spread ... was 0.20 ...
p.poldens = actUncert(p.Nl,p.polsd); 
% prepare ready-made p indices and distros, to use for Dkl calcs. Row-wise version to use
% with rowpDkl .
x = p.poldens'; 
p.polcomb = x(:) * x(:)';   % big matrix with all probs. for all policy combinations.
R = p.Nl; C=p.Nl; N=p.Nl; % For clarity - but even clearer below where p.blS is calculated !!
p.blm = N*N*repcolmat(reshape(repcolmat((0:(R*C-1)),N)',[R*N, C]),N) ; % blocked index matrix
p.lindPC = p.blm + repmat(reshape(1:(N*N),[N N]),[R C]);
p.lindPC = p.lindPC(:);  % The k_th element of this is where in the vectorized array of matrices
       % to find the k_the element of the vectorized array of the row-wise variant of said array,
       % which is good for using rowpDkl
x = zeros(size(p.polcomb(:)));
x(p.lindPC) = p.polcomb(:);
p.rowPolComb = reshape(x, [R*C, N*N]);

p.polU = p.Nl/3.75;    % This is a spread coeff. to produce similar policies out of
     % specification of peaks. If p.pol are normalized to the unit interval, 
     % 4/3.5 corresponds v. roughly to sd=0.2, like above. It is NOT a preference coeff.
     % 4/3.75 corresp. v. roughly to sd ~ 0.15. 4/6.5 would be required for ~ 0.1
p.uOC = p.polU; % Instance of the above specifically to construct Other's C-map.
     % 
p = peaks4c(p); % construct matrix p.peaks4c i.e. a large repertoire of combinations 
                % of where the peaks of marginal preferences, 'what kind of person
                % I'd like to deal with' and conditional preferences, 'what kind of
                % person I'd like to be if they have peak preference X'.
% Prevalences of each policy in a model population ...
p.CgX =[ noisyBino(2/3,1,p.Nl);  
         noisyBino(2/3,1,p.Nl);  ] ;   %  Trivially the same - here it's 1D really

% Construct the returns array with indices:
% p.ret(selfContrib, otherContrib, whoseRet, settingLevel)
p.settingLevN = 1;    % this many exchange (previously harshLevN) setting levels.
% In the first (and in HT1, only) setting level, stakes are much higher for self than Helper:
% col=helper effort:  0      1/3      2/3     3/3 
p.ret(:,:,1,1) =  [  -10       0       8      14     ;     % row 1 helpSeeker effort = 0
                     -20      -5      10      16     ;     % row 2 helpSeeker effort = 1/3
		             -35      -20     12      18     ;     %  ... etc. 
		             -50      -35     14      20       ] ; 
% Stakes for Other, from the point of view of the Helper themselves,
% so still formatted as rows=self, cols=other.
% col=helpSeeker effort:  0      1/3      2/3     3/3 
p.ret(:,:,2,1) =  [      -4       0       4       10     ;    % row 1 helper effort = 0
                         -6      -2       6       12     ;    % row 1 helper effort = 1/3
		                 -8      -4       8       14     ;    %  ... etc.
		                 -12     -8      10       16       ] ;

% In the second context level, set everything to NA here
p.ret(:,:,1,2) = nan;  % 
p.ret(:,:,2,2) = nan;  % 

% Ready - made set of indices, created with 
%  for (i=1:50); p.OpolInd(:,i) = [pBinSample(noisyBino(0.8,2,5)) , pBinSample(noisyBino(0.75,1,5)) ]'; end
% If this is commented out, it will create random samples based on p.CgX
p.setting = {}; 
p.setting{1} = ones(1,50);      p.setting{2} = ones(1,50);    p.setting{3} = ones(1,50);   
   % {1},{2}... here is the number of partners/blocks in the experiment, 50 is 'loads of trials per block'
   % They all have the same setting=1, i.e. here return map. The setting is the 
   % obvious part of the environment to the participant, i.e. here the returns map,
   % not the nature of their partner / helper.
% Set trN : 
p.nCr = 1;  % Was 5 for IC1 behavioural work-up. Number of rounds played per trial,
% functionality removed, as never used: p.rNreal = p.nCr*ones(1,50);   % explicitly all exactly nCr.

% default matrix index for Dkl about feelings about Self, 
% hence most variables in this end in [P robability] S elf
R = p.Nl; C=1; Nr = p.Nl; Nc=p.nCr; % hopefully to make role clearer ...
p.blS = Nr*Nc*repcolmat(reshape(repcolmat((0:(R*C-1)),Nr)',[R*Nr, C]),Nc) ; % blocked index matrix
p.lindPS = p.blS + repmat(reshape(1:(Nr*Nc),[Nr Nc]),[R C]);
p.lindPS = p.lindPS(:);  % The k_th element of this is where in the vectorized array of matrices
       % to find the k_the element of the vectorized array of the row-wise variant of said array,
       % which is good for using rowpDkl
x = p.poldens';
p.polPS = x(:)*ones([1,p.nCr])/p.nCr;   % big matrix with all probs. for all policy combinations.
x = zeros(size(p.polPS(:)));
x(p.lindPS) = p.polPS(:);
p.rowPolS   = reshape(x, [R*C, Nr*Nc]);

% Practice trials not processed yet:
p.rNpract = NaN; % 

% Line below hacked for HT1 -- not sure whether it covers all bases.
p.retArr = zeros(2,2,max(p.trN));  % Returns matrix will be initialized as 
p.retArr(1,1,:) = -60;        % typical IPD but has space to change 
p.retArr(1,2,:) = -15;        % if need be for each different trial.
p.retArr(2,1,:) = -90;        % Rem 1 is for D, 2 is for C.
p.retArr(2,2,:) = -45;        % So format is:     DD     DC
                              %                   CD     CC
% params to translate to actual payments :
p.roundEnd = - p.retArr(2,1,1);  % Endownment per round = r_dd_1
p.trWinCoeff = 0.002;   % by trial and error ... for main winnings
p.OinfCoeff = 0.5; % So if e.g. scored 0.75 per guess of other policy (if asked), 
                   % get 0.75*2_guesses*0.5 = 0.75
p.flat = 0.0; % flat fee component                 


% A ready-made set of policies for 2 others that can be used as 
% partners or pairs to observe:
p.partPol = zeros(4,max(p.trN));
for i=1:max(p.trN)
   p.partPol(1,i)= p.pol(1, pBinSample(p.CgX(1,:)) ); 
   p.partPol(3,i)= p.pol(1, pBinSample(p.CgX(1,:)) ); 
   p.partPol(2,i)= p.pol(1, pBinSample(p.CgX(2,:)) ); 
   p.partPol(4,i)= p.pol(1, pBinSample(p.CgX(2,:)) );  
end

% Helper policies. First, initial policies:
%               0  1/3  2/3  3/3      % Contributions of Helper
p.helperInit = [0, 0.8, 0.2,  0  ;    % Uncaring Helper init. policy
                0,  0,  0.6,  0.4 ];  % Caring Helper init. policy probabilities
% Next, policies for most trials. First, 
% Uncaring, first half of trials:
p.helperPol(:,:,1,1) = [  0.8,   0.2,   0,   0  ;     % if helpSeeker contributed 0   in previous trial
                          0.1,   0.8,  0.1,  0  ;     % if helpSeeker contributed 1/3 "    "  ...
                          0,     0.1,  0.8, 0.1 ;
                          0,     0.1,  0.8 ,0.1   ];
% Uncaring agent, second half of trials:
p.helperPol(:,:,1,2) = [  0.8,   0.2,   0,   0 ;
                          0.1,   0.8,  0.1,  0 ;
                          0.1,   0.8,  0.1,  0 ; 
                          0.1,   0.8,  0.1,  0    ]  ;
% Caring agent, first half:
p.helperPol(:,:,2,1) =   [ 0,    0.5,  0.5,  0;
                           0,    0.2,  0.8,  0;
                           0,    0,    0.5, 0.5;
                           0,    0,    0.1, 0.9   ];
% Caring agent, second half:
p.helperPol(:,:,2,2) =    p.helperPol(:,:,2,1);


%% End piece - tidy up and write copies ----------------------                
try
    save([p.data_dir p.code '_paramHT1'],'p'); 
catch
    warning(['Could not save p to ' p.data_dir p.code '_paramHT1.mat'])
end
try 
    save([p.commDir fs [p.code '_paramHT1'] '.mat'],'p'); 
catch
    warning('Could not save params. to p.commDir'); 
end

return; % end of function ----------------------------------------------



