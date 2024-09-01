function [predictY, varargout] = predict_curve(fitX, fitY, predictX, enableRobust, sigType, isMatched)
%PREDICT_CURVE Fit a trend + period curve and predict it. The trend curve is
%fitted using linear regression. The period curve is fitted using AR model
%(fillgaps in matlab)
% The shapes of fitX and fitY have to be the same.
% The fitX should be a subset of predictX!
%
% predictY = predict_curve(fitX, fitY, predictX, enableRobust, sigType)
% [predictY, trendY] = predict_curve(...) get the trend curve 
% [predictY, trendY, arY] = predict_curve(...) get both trend curve and periodic
% curve (auto-regressive model)

assert(all(size(fitX) == size(fitY)), "fitX and fitY have to be the same shape.");

arOrder = min(ceil(length(fitX) / 2), 64);
doInterp = false;
if strcmpi(sigType, "frequency")
  minRSquare = 0.92;
  maxRSquare = 0.98;
  fitY = freq_to_erb(fitY / 2 / pi * 44100);
elseif strcmpi(sigType, "amplitude")
  minRSquare = 0.86;
  maxRSquare = 0.95;
end

if enableRobust
  normalize = "on";
  robust = "LAR";
else
  normalize = "off";
  robust = "off";
end

% remove NaN or Inf
fitX = fitX(~isinf(fitY) & ~isnan(fitY));
fitY = fitY(~isinf(fitY) & ~isnan(fitY));

maxLen = []; % maximum length for estimating AR curve

if isscalar(fitX)
  trendFit = @(x) (ones(size(x)) .* fitY);
  trendGoF = struct('rsquare', 1);
elseif length(fitX) > 1
  [trendFit, trendGoF] = fit(fitX, fitY, 'poly1', Normalize=normalize, Robust=robust);
end

totalX = union(fitX, predictX);
fitXIdx = ismember(totalX, fitX);
predictXIdx = ismember(totalX, predictX);
if length(fitX) > 4 && trendGoF.rsquare > minRSquare
  trendCurve = feval(trendFit, totalX);
else
  [ipt, residual] = findchangepts(fitY(1:end), Statistic='linear', MaxNumChanges=2);
  if ~isempty(ipt) && residual < 0.05
    maxLen = max(diff([0; ipt; length(fitY)]));
  end
  trendCurve = zeros(length(totalX), 1);
end
if length(fitX) > 4 && trendGoF.rsquare > maxRSquare
  doInterp = isMatched;
  arCurve = zeros(size(trendCurve));
else
  detrendCurve = NaN(length(totalX), 1);
  detrendCurve(fitXIdx) = fitY - trendCurve(fitXIdx);
  if ~isempty(maxLen)
    arOrder = min(arOrder, maxLen);
  end
  arCurve = fillgaps(detrendCurve, maxLen, arOrder);
end

if doInterp
  predictY = interp1(fitX, fitY, totalX, 'pchip', 'extrap');
else
  predictY = arCurve + trendCurve;
end
predictY = predictY(predictXIdx);

if strcmpi(sigType, "frequency")
  predictY = audtofreq(predictY, 'erb') / 44100 * 2 * pi;
end

if nargout == 2
  varargout{1} = trendCurve(predictXIdx);
elseif nargout == 3
  varargout{1} = trendCurve(predictXIdx);
  varargout{2} = arCurve(predictXIdx);
end
