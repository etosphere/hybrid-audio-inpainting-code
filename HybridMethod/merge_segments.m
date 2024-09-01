function mergedSig = merge_segments(segments, toMergeFields)
%MERGE_SEGMENTS Merge the segments with given fields
getSegmentFieldLength = @(seg) (length(seg.(toMergeFields{1})));
segLen = cellfun(getSegmentFieldLength, segments);
mergedSig = struct;

for n = 1:length(toMergeFields)
  fieldName = toMergeFields{n};
  mergedSig.(fieldName) = zeros(sum(segLen), 1);
  endIdx = 0;
  for m = 1:length(segments)
    startIdx = endIdx + 1;
    endIdx = startIdx + segLen(m) - 1;
    mergedSig.(fieldName)(startIdx:endIdx) = segments{m}.(fieldName);
  end
end
end

