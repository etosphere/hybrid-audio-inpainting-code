function gapSegments = create_gaps(segments, gapLen, reliableLen)
%CREATE_GAPS Create gaps for each segments except the first and last one.
% The first and last segments still have the mask and gapSig, but will be
% reliable for all samples.
segLen = length(segments{1}.originSig);
assert(segLen > 2 * reliableLen + gapLen, ...
  'Signal length for segments is too short or gap length is too long.')

gapSegments = segments;

for n = 1:length(segments)
  gapSegments{n}.mask = true(size(gapSegments{n}.originSig));
  gapSegments{n}.gapSig = gapSegments{n}.originSig;
  if ~(n == 1 || n == length(segments))
    gapStartIdx = reliableLen + round(rand() * (segLen - reliableLen * 2 - gapLen + 1));
    gapEndIdx = gapStartIdx + gapLen - 1;
    gapSegments{n}.mask(gapStartIdx:gapEndIdx) = false;
    gapSegments{n}.gapSig(~gapSegments{n}.mask) = 0;
  end
end
