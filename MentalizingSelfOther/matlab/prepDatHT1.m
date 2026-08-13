function [v, d, p] = prepDatHT1( p, d)
% Set up the data structures for storing stuff
% for task Help-seeking Task mk 1 and fill in pt code and date.
% Also initialise key fields of the variable structure v

try
    if ~d.forReal 
        error('d.forReal must not be false')
    end
catch
    d = struct();
    d.forReal = 1; 
end

d.code = p.code; 
d.PID =  p.PID;
d.date = clock;
d.evo_hd = {'trial',   'round',   'trType',  'Spol',  'tStim1',...            % will be col 1 to 5 of d.evo
            'tDec1',   'RT', 'sWork', 'retS',    'retO', ...                  % 6 to 10
            'Seval',   'Oeval',   'oSeval',  'infOpola', 'infOpolb',...       % 11 to 15
            'tScontrOn','tScontrOff','tOcontrOn','tOcontrOff','Opol',...      % 16 to 20
            'oWork', 'Sfee', 'Ofee', 'tSevalOn', 'tSevalOff',...              % 21 to 25; 24-29 eval timings
            'RTSeval',  'tOevalOn', 'tOevalOff', 'RTOeval',  'partnerType', ... % 26 to 30; 
            'slStim1', 'slDec1', 'slScontrOn', 'slScontrOff', 'slOcontrOn' ,...  % 31 to 35; 31-40 event slices.
            'slOcontOff', 'slSevalStim', 'slSeval', 'slOevalStim', 'slOeval',... % 36 to 40                        
             'tSEvalOptions', 'tOEvalOptions'};                                  % 41 - 
d.evo = {};
d.feelm = []; d.evalRL=[]; % Declarations for person evaluation / interpersonal feeling parameters.

p.evoSt = {};

for helperN = 1:p.helperNum
    d.evo{helperN} = zeros(1,length(d.evo_hd)); % time evolution of task, (trial, round, data). recorded are:
      % index:   1        2       3       4        5        6      7      8             9
      % var      trial    round   trType  Spol     tStim1   tDec1  RT    tWhichContext  retS
      %  ...     10       11      12      13       14       15        16       17       18      19
      %  ...     retO     Seval   Oeval   oSeval   infOCCa  infOCCb  infODDa   infODDb  phys1  phys2
      %  ...     20       21      22      23       24       25
      %  ...     OpolCC   OpolDD  Sfee    Ofee     tSevalOn tSevalOff 
      %  ...     26       27        28       29      30
      %  ...     tSeval   tOevalOn tOevalOff tOeval partnerType
      %  ...     31       32      33       34          35
      %  ...     slStim1 slDec1 slScontrOn slScontrOff slOcontrOn      
      %  ...     36         37          38      39          40
      %  ...     slOcontOff slSevalStim slSeval slOevalStim slOeval
      % Created with one 'row' and 'col' only, will have 'rows' added as trials are completed.
      % tStim1,2 will record when first and second of either questions or moves
      % are displayed. 
      p.evoSt{helperN} = []; % record starting row for this trial in evo matrix
end

%% Now fill in with the very basics:
for helperN = 1:p.helperNum
    trialN = p.trN(helperN); 
    evoSt=1;  % Starting row for each array that we are filling in.

    for tr = 1:trialN
        p.evoSt{helperN}(tr) = evoSt;
        
        d.evo{helperN}(evoSt:(evoSt+p.nCr+1),1)  = tr;
        d.evo{helperN}(evoSt:(evoSt+p.nCr+1),3)  = p.settingl;  % trivial trial type
        % Here, the 'partnerType' records whether the Helper was uncaring=1 or caring=2 :
        d.evo{helperN}(evoSt:(evoSt+p.nCr+1),30) = p.helperType(helperN,1);
        d.evo{helperN}(evoSt:(evoSt+p.nCr-1),2)  = 1:p.nCr; % trivial - which interaction within-trial.
        d.evo{helperN}(evoSt+p.nCr,2) = 901;   % Code for self-policy data etc, taken pre-trial
        d.evo{helperN}(evoSt+p.nCr+1,2) = 902; % person-eval and inference Q's, taken post-trial
    
        % startpoint of next trial:
        evoSt = evoSt + p.nCr + 2;
    
    end % end loop over trials in each helper block.

end % end loop over different helpers (big blocks)

% Also record the rows in d.evo where the questions about the
% policy and about the (evaluations+inferences) are recorded:
p.polRow = {};
p.QRow = {}; 
for helperN = 1:p.helperNum
    [p.polRow{helperN}, ~] = find(d.evo{helperN}==901);
    [p.QRow{helperN}, ~] = find(d.evo{helperN}==902);
end

% d.expStartTime = time/p.tUnitsPerSec;  % not here but when we get
% initial input from user - e.g. confirmation of keys used just when ready4IC is called.

% initialise a number of 'var' variables and obtain codes for which keys to use.
% rand('seed',sum(100*clock)); not here !

%%  Initialize key parts of the variables structure:

v.assets = [0,0]; % Amounts accumulated by S and PT respectively.
v.trial = 0; % current trial (aka game, which in HT1 only has 1 round or interaction)
v.practTrDone = 0; % number of practice trials done 
     % - when this goes to pars.practiceTrN, real task begins ...
     % adjusted separately if genMod is provided, in prep4PractOrMain...
v.Spol = nan(1,p.Nl);        % Self (in HT1, help-seeker) policy (as opposed to action).
v.Ppol = nan(1,p.Nl);        % Partner (other, here Helper) policy (as opposed to action).
v.settingLev = p.invalid;         % current setting / context level
v.moves = nan(2,p.nCr,p.settingl);  % most recent moves. mid dim is not so i.
v.prMoves = v.moves;            % moves of trial before last played
v.rets  = v.moves;
v.prRets = v.rets;
v.wait4pt = nan;   % nan = actively wait for participant to OK having read instructions.

% This fn. is called before the expt. starts. d.forReal is already non-zero in HT1,
% i.e. no practice / prelim. trials processed here, so set 
% the v.practTrDone to invalid flag:
if d.forReal
    v.practTrDone=p.invalid; 
else
    error('No practice trials should be processed here')'
end

return;

