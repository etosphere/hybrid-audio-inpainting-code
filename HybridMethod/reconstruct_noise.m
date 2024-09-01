function recSig = reconstruct_noise(sig, w, mask, type)
%RECONSTRUCT_NOISE Summary of this function goes here
%   Detailed explanation goes here
%
% mask: the reliable mask, set 1 for all reliable samples
% type: 1 for cross-fade, 2 for LAR interpolation

%% Initialize parameters
gapStartIdx = find(mask == 0, 1, 'first');
gapEndIdx = find(mask == 0, 1, 'last');
gapLength = gapEndIdx - gapStartIdx + 1;

%% auto-regression model approach
polesNum = 256;
crossFadeLen = 32;
leftReliableSig = sig(1:gapStartIdx-1);
rightReliableSig = sig(gapEndIdx+1:end);
[leftReliableLpcCoef, leftReliableLpcGain, leftReflectionCoef] = arburg(leftReliableSig, polesNum);
[rightReliableLpcCoef, rightReliableLpcGain, rightReflectionCoef] = arburg(rightReliableSig, polesNum);
noiseLen = gapLength + polesNum + crossFadeLen * 2;

if type == 1
  noiseSig1 = randn(noiseLen, 1);
  noiseSig1 = (noiseSig1 - mean(noiseSig1)) ./ std(noiseSig1);
  noiseSig2 = randn(noiseLen, 1);
  noiseSig2 = (noiseSig2 - mean(noiseSig2)) ./ std(noiseSig2);
  leftNoise = filter(sqrt(leftReliableLpcGain), leftReliableLpcCoef, noiseSig1);
  leftNoise = leftNoise(numel(leftReliableLpcCoef):end);
  leftNoise = leftNoise .* (std(leftReliableSig) / std(leftNoise));
  rightNoise = filter(sqrt(rightReliableLpcGain), rightReliableLpcCoef, noiseSig2);
  rightNoise = rightNoise(numel(rightReliableLpcCoef):end);
  rightNoise = rightNoise .* (std(rightReliableSig) / std(rightNoise));
  recNoiseSig = sqrt(rampup(length(rightNoise))) .* rightNoise + sqrt(rampdown(length(leftNoise))) .* leftNoise;
  recNoiseSig = rampsignal(recNoiseSig, crossFadeLen);
elseif type == 2
  leftReliableLAR = rc2lar(leftReflectionCoef);
  rightReliableLAR = rc2lar(rightReflectionCoef);

  % interpolate the LSF and LAR coefficients
  winLen = polesNum * 2;
  noiseLen = ceil(noiseLen / winLen) * winLen;
  noiseSig = randn(noiseLen, 1);
  recNoiseSig = zeros(size(noiseSig));
  interpLen = noiseLen / winLen; % the first buffer-filling frame and the last gap frame will use same coefficients as left and right reliable parts

  interpLAR = interp1([1, interpLen-1], [leftReliableLAR, rightReliableLAR].', [1:interpLen-1, interpLen-1]).';
  interpLpcCoef = zeros(size(interpLAR, 2), size(interpLAR, 1) + 1);
  for n = 1:size(interpLAR, 2)
    interpLpcCoef(n, :) = rc2poly(lar2rc(interpLAR(:, n)));
  end

  interpGain = interp1([1, interpLen], [leftReliableLpcGain, rightReliableLpcGain], 1:interpLen);
  state = zeros(polesNum, 1);
  for n = 1:interpLen
    idx = (1:winLen) + (n - 1) * winLen;
    [recNoiseSig(idx), state] = filter(sqrt(interpGain(n)), interpLpcCoef(n, :), noiseSig(idx), state);
  end
  recNoiseSig = recNoiseSig((1:gapLength+2*crossFadeLen) + polesNum);
  recNoiseSig = rampsignal(recNoiseSig, crossFadeLen);
end

recSig = sig;
recSig(~mask) = 0;
recSig(gapStartIdx-crossFadeLen:gapStartIdx-1) = recSig(gapStartIdx-crossFadeLen:gapStartIdx-1) .* sqrt(rampdown(crossFadeLen));
recSig(gapEndIdx+1:gapEndIdx+crossFadeLen) = recSig(gapEndIdx+1:gapEndIdx+crossFadeLen) .* sqrt(rampup(crossFadeLen));
recSig(gapStartIdx-crossFadeLen:gapEndIdx+crossFadeLen) = recSig(gapStartIdx-crossFadeLen:gapEndIdx+crossFadeLen) + recNoiseSig;
end