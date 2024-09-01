function optimalLambda = get_optimal_lambda(insig, mask, sigType)
sigRMSEnergy = sqrt(norm(insig(mask))^2 / sum(mask));

if sigType == 1 || strcmp(sigType, 'long')
  Fflatness = frametight(frame('dgtreal', gabwin({'tight', 'hann'}, 32, 128, 128), 32, 128, 'freqinv'));

  sigCoefTf = framecoef2tf(Fflatness, frana(Fflatness, insig));
  sigCoefTf(:, 1:2) = 0;
  sigCoefTf(:, end-1:end) = 0;
  [gradX, gradY] = gradient(abs(sigCoefTf));
  gradDistance = hypot(gradX, gradY);
  flatness.origin = generalized_spectral_flatness(abs(sigCoefTf));
  flatness.gradDistance = generalized_spectral_flatness(gradDistance);

  firstNanIdx = find(isnan(flatness.origin(3:end-2)), 1, 'first') + 2;
  lastNanIdx = find(isnan(flatness.origin(3:end-2)), 1, 'last') + 2;
  flatness.origin(firstNanIdx-3:lastNanIdx+3) = NaN;
  flatness.gradDistance(firstNanIdx-5:lastNanIdx+5) = NaN;

  meanFlatnessGradDist = mean(flatness.gradDistance, "omitmissing");
  
  optimalLambda = 10 ^ (1.739988 * log10(meanFlatnessGradDist * sigRMSEnergy) + 0.225884);
elseif sigType == 2 || strcmp(sigType, 'short')
  Ftrans = frametight(frame('dgtreal', gabwin({'tight', 'hann'}, 16*2, 64*2, 64*2), 16*2, 64*2));
  resCoefTf = framecoef2tf(Ftrans, frana(Ftrans, insig));
  resCoefTf(:, 1:2) = 0;
  resCoefTf(:, end-1:end) = 0;
  flatness.origin = generalized_spectral_flatness(abs(resCoefTf));
  firstNanIdx = find(isnan(flatness.origin(3:end-2)), 1, 'first') + 2;
  lastNanIdx = find(isnan(flatness.origin(3:end-2)), 1, 'last') + 2;
  [~, gradY] = gradient(abs(resCoefTf));
  flatness.gradY = generalized_spectral_flatness(abs(gradY));
  flatness.gradY(firstNanIdx-2:lastNanIdx+2) = NaN;
  prc90Flatness = prctile(flatness.gradY(~isnan(flatness.gradY)), 90);

  optimalLambda = 10 ^ (1.209760 * log10(prc90Flatness * sigRMSEnergy) + 0.547436);
end
end
