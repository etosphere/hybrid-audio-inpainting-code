function segments = create_segments(sig, segNum)
sigLen = length(sig);
segments = cell(segNum, 1);
segLen = round(sigLen / segNum);
lastSegLen = sigLen - segLen * (segNum - 1);
for n = 1:segNum
  if n < segNum
    idx = (1:segLen) + segLen * (n - 1);
  else
    idx = (1:lastSegLen) + segLen * (segNum - 1);
  end
  segments{n} = struct('originSig', sig(idx));
end

end
