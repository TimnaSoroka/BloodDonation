clear
% in hol_data there are all raw files excpect: 13,34,57 that did not donated

TBD=readcell('/Users/timnas/Documents/projects/BloodDonation/BloodDonationTable102.xlsx');
sex_idx=find(strcmpi(TBD(:,1),'sex'));

    numericVec = cellfun(@(x) double(x), TBD(1,:));

% Set directory containing CSV files
backupDir = '/Users/timnas/Documents/projects/BloodDonation/excluded_or_editted_data/';
csvDir= '/Users/timnas/Documents/projects/BloodDonation/hol_data/';
files = dir(fullfile(csvDir, '*.csv'));

% Step 1: QA- check that data were not duplicated
for i = 1:length(files)
    filePath = fullfile(csvDir, files(i).name);
    
    % Read first 10 lines to check for duplication
             opts = detectImportOptions(filePath, 'TextType', 'string');
                opts = setvaropts(opts, 'comment', 'Type', 'string'); % Ensure "comments" is read as text
                opts.VariableNamingRule='preserve';
                T = readtable(filePath, opts);

        previewData=T(1:10,:);
    if size(previewData, 1) < 2
        fprintf('File %s too short to check.\n', files(i).name);
        continue;
    end
    
    % Compare first row with next rows
    repeatedRows = 0;
    firstRow = previewData(1, :);
    for r = 2:size(previewData,1)
        if isequal(firstRow(1,2:4), previewData(r, 2:4))
            repeatedRows = repeatedRows + 1;
        else
            break;
        end
    end
    
    if repeatedRows > 0
        fprintf('File %s: %d duplicated rows detected at the top. Fixing...\n', files(i).name, repeatedRows);
        
        % Load full data and remove duplicated rows
        fullData = readtable(filePath,opts);
        cleanedData = fullData(1:repeatedRows+1:end, :);
        
        % Move old file to backup
        movefile(filePath, fullfile(backupDir, files(i).name));
        cleanedPath = fullfile(csvDir, [files(i).name(1:end-4) '_cleaned.csv']);
        writetable(cleanedData, cleanedPath);
    end
end


% Step 2: Group files by participant code (characters 5-7)
participantMap = containers.Map();

for i = 1:length(files)
    name = files(i).name;
    if length(name) >= 7
        code = name(5:7);  % Extract participant code
        filePath = fullfile(csvDir, name);
        if isKey(participantMap, code)
            currentList = participantMap(code);
            currentList{end+1} = filePath;
            participantMap(code) = currentList;  % Assign back
        else
            participantMap(code) = {filePath};
        end
    end
end

% Step 3: Process each participant
participantKeys = keys(participantMap);

for i = 1:length(participantKeys)
    code = participantKeys{i};
    fileList = participantMap(code);

    % Concatenate tables from all files of the same participant
    fullData = table();
    for j = 1:length(fileList)
        try
            opts = detectImportOptions(fileList{j}, 'TextType', 'string');
                opts = setvaropts(opts, 'comment', 'Type', 'string'); % Ensure "comments" is read as text
                opts.VariableNamingRule='preserve';
                T = readtable(fileList{j}, opts);
            fullData = [fullData; T];
        catch ME
            fprintf('Error reading file %s: %s\n', fileList{j}, ME.message);
        end
    end

    % Display basic info
    fprintf('\nParticipant: %s\n', code);
        subjData(i).code=code;
        subjData(i).Data=fullData;
    subjData(i).length=height(fullData);
    [subjData(i).in,subjData(i).out,validated_vals{i}]= extract_timings(fullData,code);    % Check for 'comment' column
    
    sbj_idx=find(numericVec==str2double(code));

    % add sex
    S=TBD(sex_idx,sbj_idx);
    if strcmpi(S{:},'זכר')
    subjData(i).sex=0;
    elseif strcmpi(S{:},'נקבה')
        subjData(i).sex=1;
    else
        fprintf('No gender match found for subject ID: %s \n', code)
    end

    % Add donation amount
    rowName=strcmpi('Pure_Donation_Amount',(TBD(:,1)));
tmp=TBD(rowName,sbj_idx);
 subjData(i).Donation_Amount=tmp{:};
    % Add weight
    rowName=strcmpi('Weight',(TBD(:,1)));
    tmp=TBD(rowName,sbj_idx);
    if isa(tmp{:},'double')
 subjData(i).Weight=tmp{:};
    else
        disp(tmp{:})
         subjData(i).Weight=nan;
    end
     % Add pulse1
    rowName=strcmpi('Pulse1st',(TBD(:,1)));
    tmp=TBD(rowName,sbj_idx);
 subjData(i).Pulse1st=tmp{:};

      % Add pulse2
    rowName=strcmpi('pulse2nd',(TBD(:,1)));
    tmp=TBD(rowName,sbj_idx);
 subjData(i).Pulse1nd=tmp{:};

       % Add Spo1
    rowName=strcmpi('Spo1st',(TBD(:,1)));
    tmp=TBD(rowName,sbj_idx);
 subjData(i).Spo1st=tmp{:};

        % Add Spo2
    rowName=strcmpi('Spo2nd',(TBD(:,1)));
    tmp=TBD(rowName,sbj_idx);
 subjData(i).Spo2nd=tmp{:};

    clearvars T fullData
end



stringMatrix = reshape([validated_vals{:}],2, []);
save('/Users/timnas/Documents/projects/BloodDonation/Holter_timings.mat','subjData');

function [in,out,validated_vals]= extract_timings(fullData,code)
if any(strcmpi(fullData.Properties.VariableNames, 'comment'))
    comments = fullData.comment;
elseif any(strcmpi(fullData.Properties.VariableNames, 'x_comment_'))
    comments = fullData.x_comment_;
else
    disp('No "comment" column found.');
end
in = find(contains(lower(comments), 'needle in', 'IgnoreCase', true));
out = find(contains(lower(comments), 'needle out', 'IgnoreCase', true));
if isempty(in) ||     isempty(out)
    [comm_sorted,comm_idx]=sort(comments);
    switch code
        case '001'
            in=comm_idx(8);
            out=comm_idx(2);
        case '002'
            in=comm_idx(6);
            out=comm_idx(5);
       case '076'
[comm_sorted,comm_idx]=sort(comments,'descend');
            in=comm_idx(8);
            out=comm_idx(6);
      case '078'
[comm_sorted,comm_idx]=sort(comments,'descend');
            in=comm_idx(6);
            out=comm_idx(5);
    case '080'
            [comm_sorted,comm_idx]=sort(comments,'descend');
            in=comm_idx(6);
            out=comm_idx(5);
                    case '030'
            in=comm_idx(11);
            out=comm_idx(12);
        case '031'
        in=comm_idx(9);
        out=comm_idx(10);
        otherwise
            comm_sorted(1:15)
    end
elseif size(in,1)>1 || size(out,1)>1
    switch code
        case '003'
            out=out(2);
                  case '077'
            in=in(2);
            case '029'
            out=out(2);
            case '055'
            out=out(2);
            case '065'
            out=out(2);
              case '014'
            in=in(2);
              case '056'
            in=in(2);
             case '042'
            in=in(1)-20*25;
              case '025'
            in=in(2);
            out=out(2);
                       case '095'
            out=in(2);
                        in=in(1);
                         case '103'
            out=out(1)-(20*25);
        otherwise
            disp(comments(in))
                        disp(comments(out))

    end

end
validated_vals={comments(in),comments(out)};
end
