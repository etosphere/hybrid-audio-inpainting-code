function flatness = generalized_spectral_flatness(coef)
%GENERALIZED_SPECTRAL_FLATNESS Calculate the "spectral flatness" for any
% given matrix.

freqBinNum = size(coef, 1);
coefAmp = abs(coef);
flatness = (exp(sum(log(coefAmp + 1e-15), 1) / freqBinNum)) ./ (sum(coefAmp, 1) / freqBinNum);
flatness(isinf(flatness)) = NaN;
end