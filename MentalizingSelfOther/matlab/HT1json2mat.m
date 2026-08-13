function T = HT1json2mat(filepath, pid, task_design, renumber_blocks, outCsv)
% Parse one raw Helping Task Pavlovia/jsPsych export into a clean MATLAB table.
% MATLAB translation of Tyler's R parser, keeping column naming close to his output.
% If outCsv is empty, no csv file output is written.
% If renumber_blocks is empty or missing, then block numbering starts from the real,
%    not the practice, data.
% If pid is not given, an attempt is made to set it from the (hopefully) timestamp, last
%    section of the filepath string.

if nargin < 4 || isempty(renumber_blocks)
    renumber_blocks = true;
end
if nargin < 5
    outCsv = '';
end
if nargin < 3 || isempty(task_design)
    task_design = 'unknownDesign';
end
if nargin < 2 || isempty(pid)
   try  % Extract the month onwards timestamp as a convenient (almost random) ID number:
     pat = '\d{4}-(\d{1,2})-(\d{1,2})_(\d{1,2})h(\d{1,2})\.(\d{1,2})\.(\d{1,3})\.csv$';
     tok = regexp(filepath, pat, 'tokens', 'once');
     pid = strjoin(string(tok), '');
   catch
       pid = filepath;
       warning('pid setting did not work well - set to the whole filepath');
   end
end

lines = readlines(filepath, 'EmptyLineRule', 'skip');
lines = string(lines);

prolificId = extract_prolific_id(lines);

headerIdx = find(contains(lines, 'trial_history_columns'), 1, 'first');
if isempty(headerIdx)
    error("Could not find a 'trial_history_columns' line in %s.", filepath);
end
headerLine = lines(headerIdx);
headerTok = regexp(headerLine, 'trial_history_columns"":\s*""([^\"]*)""', 'tokens', 'once');
if isempty(headerTok)
    error('Found a trial_history_columns line, but could not extract the header string.');
end
colNames = split(string(headerTok{1}), ',');
colNames = strtrim(colNames(:));
nCols = numel(colNames);

candidateMask = ~cellfun('isempty', regexp(cellstr(lines), '^\s*""[0-9]+;', 'once'));
candidateLines = lines(candidateMask);
semiCounts = count(candidateLines, ';');
trialLines = candidateLines(semiCounts == (nCols - 1));
if isempty(trialLines)
    error('Found the header but no trial rows with a matching number of fields.');
end

parsed = cell(numel(trialLines), 1);
okLength = false(numel(trialLines), 1);
for i = 1:numel(trialLines)
    inner = regexprep(trialLines(i), '^\s*""(.*)""\s*,?\s*$', '$1');
    parts = split(string(inner), ';');
    parsed{i} = parts(:)';
    okLength(i) = numel(parts) == nCols;
end
parsed = parsed(okLength);
if isempty(parsed)
    error('Trial rows were found, but none had the expected number of fields.');
end

rawMat = strings(numel(parsed), nCols);
for i = 1:numel(parsed)
    rawMat(i, :) = parsed{i};
end

safeNames = matlab.lang.makeValidName(cellstr(colNames), 'ReplacementStyle', 'delete');
d = array2table(rawMat, 'VariableNames', safeNames);
nameMap = containers.Map(cellstr(colNames), safeNames);

numericCols = { ...
    'overall_idx', 'block_num', 'block_idx', ...
    'effort_selection_rxn_time', 'effort_points', 'partner_points', ...
    'rate_self_rxn_time', 'rate_self_value', ...
    'rate_partner_rxn_time', 'rate_partner_value'};
for i = 1:numel(numericCols)
    cc = numericCols{i};
    if isKey(nameMap, cc)
        vn = nameMap(cc);
        d.(vn) = str2double(d.(vn));
    end
end

blockNum = getVar(d, nameMap, 'block_num');
overallIdx = getVar(d, nameMap, 'overall_idx');
if renumber_blocks
    [~, ~, blockOut] = unique(blockNum, 'stable');
else
    blockOut = blockNum;
end

helperPolicy = strings(height(d),1);
if isKey(nameMap, 'caring/uncaring')
    helperPolicy = map_policy_labels(d.(nameMap('caring/uncaring')));
end

effortIntended = map_effort_labels(getVarStr(d, nameMap, 'effort_value_selected'));
effortActual   = map_effort_labels(getVarStr(d, nameMap, 'effort_value_actual'));
partnerEffort  = map_effort_labels(get_partner_effort_raw(d, nameMap));

T = table( ...
    repmat(string(pid), height(d), 1), ...
    repmat(string(task_design), height(d), 1), ...
    repmat(prolificId, height(d), 1), ...
    overallIdx, ...
    blockOut, ...
    helperPolicy, ...
    effortIntended, ...
    effortActual, ...
    partnerEffort, ...
    getVar(d, nameMap, 'effort_selection_rxn_time'), ...
    getVar(d, nameMap, 'rate_self_rxn_time'), ...
    getVar(d, nameMap, 'rate_partner_rxn_time'), ...
    getVar(d, nameMap, 'rate_self_value'), ...
    getVar(d, nameMap, 'rate_partner_value'), ...
    getVar(d, nameMap, 'partner_points'), ...
    getVar(d, nameMap, 'effort_points'), ...
    'VariableNames', { ...
        'PID', 'Task_Design', 'ProlificId', 'Trial', 'Block', 'Helper_Policy', ...
        'Effort_Intended', 'Effort_Actual', 'Partner_Effort', ...
        'RT_Effort_Decision_ms', 'RT_Satisfaction_Rating_ms', 'RT_Trust_Rating_ms', ...
        'Satisfaction_Rating', 'Trust_Rating', 'Points_For_Helper', 'Points_For_Help_Seeker'});

T.Properties.VariableNames = { ...
    'PID', 'Task Design', 'ProlificId', 'Trial', 'Block', 'Helper Policy', ...
    'Effort Intended', 'Effort Actual', 'Partner Effort', ...
    'RT Effort Decision (ms)', 'RT Satisfaction Rating (ms)', 'RT Trust Rating (ms)', ...
    'Satisfaction Rating', 'Trust Rating', 'Points For Helper', 'Points For Help-Seeker'};

if ~isempty(outCsv)
    writetable(T, outCsv, 'QuoteStrings', true);
end
end

function prolificId = extract_prolific_id(lines)
prolificId = missing;

% pat = '"participant"\s*:\s*"([^"]+)"';
pat = '""participant""\s*:\s*""([^"]+)""';
for i = 1:numel(lines)
    tok = regexp(lines(i), pat, 'tokens', 'once');
    if ~isempty(tok)
        prolificId = string(tok{1});
        return;
    end
end
end


function x = get_partner_effort_raw(T, nameMap)
preferred = {'partner_value_actual','partner_effort_value_actual','partner_effort_value','effort_value_partner','partner_effort_actual'};
for i = 1:numel(preferred)
    rawName = preferred{i};
    if isKey(nameMap, rawName)
        x = string(T.(nameMap(rawName)));
        return;
    end
end

altNames = {'show_partner_effort_results','showpartnereffortresults'};
for i = 1:numel(altNames)
    rawName = altNames{i};
    if isKey(nameMap, rawName)
        raw = string(T.(nameMap(rawName)));
        x = extract_effort_token(raw);
        return;
    end
end

x = strings(height(T), 1);
x(:) = missing;
end

function out = extract_effort_token(raw)
s = lower(erase(string(raw), '"'));
out = strings(size(s));
out(:) = missing;

for i = 1:numel(s)
    row = s(i);
    tok = missing;
    if contains(row, 'maximum')
        tok = "maximum";
    elseif contains(row, 'zero')
        tok = "zero";
    else
        m = regexp(char(row), '(?<!\d)(23|13)(?!\d)', 'match', 'once');
        if ~isempty(m)
            tok = string(m);
        elseif contains(row, '2/3')
            tok = "2/3";
        elseif contains(row, '1/3')
            tok = "1/3";
        end
    end
    out(i) = tok;
end
end

function x = getVar(T, nameMap, rawName)
if isKey(nameMap, rawName)
    x = T.(nameMap(rawName));
else
    x = nan(height(T), 1);
end
end

function x = getVarStr(T, nameMap, rawName)
if isKey(nameMap, rawName)
    x = string(T.(nameMap(rawName)));
else
    x = strings(height(T), 1);
    x(:) = missing;
end
end

function out = map_effort_labels(x)
s = strtrim(string(x));
s = erase(s, '"');
sLower = lower(strtrim(s));
out = nan(size(s));
out(sLower == "zero") = 1;
out(sLower == "13") = 2;
out(sLower == "23") = 3;
out(sLower == "maximum") = 4;
out(sLower == "minimum") = 1;
out(sLower == "1/3") = 2;
out(sLower == "2/3") = 3;
out(sLower == "1") = 1;
out(sLower == "2") = 2;
out(sLower == "3") = 3;
out(sLower == "4") = 4;
end

function out = map_policy_labels(x)
s = string(x);
out = s;
out(s == "C") = "CARING";
out(s == "U") = "UNCARING";
end
