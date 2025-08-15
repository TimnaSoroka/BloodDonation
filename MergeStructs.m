function mergedStruct = MergeStructs(varargin)
%% Validate input
input = varargin;
inputLength = numel(input);
if inputLength == 1
    error('MergeStructs should have at least 2 input items');
end

vectorLegnthForEachInput = cellfun(@numel, input);
uniqueVectorLegnthForEachInput = unique(vectorLegnthForEachInput);
if numel(uniqueVectorLegnthForEachInput) > 1
    error('MergeStructs should have the same vector-length in each input item');
end

%% Merge inputs
mergedStruct = input{1};
for inputIndex = 2:inputLength
    currentInputItemsToMerge = input{inputIndex};
    currentFieldsToMerge = fieldnames(currentInputItemsToMerge);
    for itemIndex = 1:uniqueVectorLegnthForEachInput
        for i = 1:length(currentFieldsToMerge)
            if ~isfield(mergedStruct(itemIndex), currentFieldsToMerge{i}) || isempty(mergedStruct(itemIndex).(currentFieldsToMerge{i}))
                mergedStruct(itemIndex).(currentFieldsToMerge{i}) = currentInputItemsToMerge(itemIndex).(currentFieldsToMerge{i});
            elseif ~isequal(mergedStruct(itemIndex).(currentFieldsToMerge{i}), currentInputItemsToMerge(itemIndex).(currentFieldsToMerge{i}))
                warning('Non-equal values for field %s for item #%d in input variable #%d', currentFieldsToMerge{i}, itemIndex, inputIndex);
            end
        end
    end
end