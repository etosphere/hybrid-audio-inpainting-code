function [predictedPartials, partialSig] = partial_prediction(leftPartialIdx, rightPartialIdx, totalFrameNum, matchedPartials, firstGapFrameIdx, extraFrameNum, lastGapFrameIdx, fs, N, HopFactor, minSlopeRatio, maxBirthLen, maxDeltaAmp, minAmpThreshold, maxDeathLen, timeLoc, sigLen, paddingLen)
matchedPartialIdx = intersect(leftPartialIdx, rightPartialIdx);
unmatchedBirthPartialIdx = setdiff(rightPartialIdx, leftPartialIdx);
unmatchedDeathPartialIdx = setdiff(leftPartialIdx, rightPartialIdx);
mergedPartialIdx = [matchedPartialIdx; unmatchedBirthPartialIdx; unmatchedDeathPartialIdx];
mergedPartialType = [ones(size(matchedPartialIdx)); 2 * ones(size(unmatchedBirthPartialIdx)); 3 * ones(size(unmatchedDeathPartialIdx))]; % 1: matched, 2: unmatched birth, 3: unmatched death
gappedPartials = cellfun(@(~)double.empty(0, 3), cell(length(mergedPartialIdx), 1), 'UniformOutput', false);

for n = 1:totalFrameNum
  for k = 1:length(mergedPartialIdx)
    rowIdx = find(matchedPartials{n}(:, 4) == mergedPartialIdx(k), 1);
    if isempty(rowIdx)
      gappedPartials{k} = [gappedPartials{k}; NaN(1, 3)];
    else
      gappedPartials{k} = [gappedPartials{k}; matchedPartials{n}(rowIdx, 1:3)];
    end
  end
end

% add extra NaNs
for n = 1:length(gappedPartials)
  firstValidIdx = min(max(4, find(~isnan(gappedPartials{n}(:, 1)), 1, 'first')), firstGapFrameIdx);

  gappedPartials{n}(1:firstValidIdx - 1, 1:2) = NaN;
  gappedPartials{n}(end-3:end, 1:2) = NaN;
end

wrapToPi = @(phi)(phi - 2 * pi * floor((phi + pi)/(2 * pi)));

% predict the partials in the gap
for n = 1:length(gappedPartials)
  % predict frequency for both left and right parts
  predictX = (firstGapFrameIdx - extraFrameNum:lastGapFrameIdx + extraFrameNum)';
  notNanArray = ~isnan(gappedPartials{n}(:, 2));
  leftFitX = find(notNanArray(1:firstGapFrameIdx-1), 1, 'first'):find(notNanArray(1:firstGapFrameIdx-1), 1, 'last');

  % remove attack
  [changePoint, slope, intercept] = ischange(gappedPartials{n}(leftFitX, 2)', 'linear', MaxNumChanges = 1);
  changePointIdx = find(changePoint, 1);
  if ~isempty(changePointIdx) % might include an attack
    attackSlope = polyfit(1:changePointIdx, gappedPartials{n}(leftFitX(1:changePointIdx), 2), 1);
    attackSlope = attackSlope(1);
    decaySlope = polyfit(1:numel(leftFitX)-changePointIdx+1, gappedPartials{n}(leftFitX(changePointIdx:end), 2), 1);
    decaySlope = decaySlope(1);
    if attackSlope > -decaySlope * minSlopeRatio && decaySlope < 0
      while (gappedPartials{n}(leftFitX(changePointIdx), 2) < gappedPartials{n}(leftFitX(changePointIdx+1), 2))
        changePointIdx = changePointIdx + 1;
      end
      leftFitX = leftFitX(changePointIdx+1:end);
    end
  end
  rightFitX = (find(notNanArray(lastGapFrameIdx+1:end), 1, 'first'):find(notNanArray(lastGapFrameIdx+1:end), 1, 'last')) + lastGapFrameIdx;

  if (mergedPartialType(n) ~= 1)
    % predict frequency for both left and right parts
    if ~isempty(leftFitX)
      leftFreqPrediction = predict_curve(leftFitX', gappedPartials{n}(leftFitX, 1), ...
        predictX, false, "frequency", false);
    end
    if ~isempty(rightFitX)
      rightFreqPrediction = predict_curve(rightFitX', gappedPartials{n}(rightFitX, 1), ...
        predictX, false, "frequency", false);
    end
  end
  freqFitX = [leftFitX, rightFitX];

  % remove frames with unreliable amplitudes for prediction
  if ~isempty(leftFitX)
    if (leftFitX(1) <= firstGapFrameIdx - 1 - extraFrameNum)
      leftFitX = leftFitX(leftFitX <= firstGapFrameIdx-1-extraFrameNum);
    else
      leftFitX = leftFitX(1);
    end
    if (mergedPartialType(n) ~= 1)
      [leftAmpPrediction, leftAmpTrend, leftAmpAr] = predict_curve(leftFitX', gappedPartials{n}(leftFitX, 2), ...
        predictX, true, "amplitude", false);
    end
  end
  
  rightFitX = (find(notNanArray(lastGapFrameIdx+1:end), 1, 'first'):find(notNanArray(lastGapFrameIdx+1:end), 1, 'last')) + lastGapFrameIdx;
  if ~isempty(rightFitX)
    if (rightFitX(end) >= lastGapFrameIdx + 1 + extraFrameNum)
      rightFitX = rightFitX(rightFitX >= lastGapFrameIdx+1+extraFrameNum);
    else
      rightFitX = rightFitX(end);
    end
    if (mergedPartialType(n) ~= 1)
      [rightAmpPrediction, rightAmpTrend, rightAmpAr] = predict_curve(rightFitX', gappedPartials{n}(rightFitX, 2), ...
        predictX, true, "amplitude", false);
    end
  end
  ampFitX = [leftFitX, rightFitX];

  if (mergedPartialType(n) == 1) % matched partials
    freqPrediction = predict_curve(freqFitX', gappedPartials{n}(freqFitX, 1), ...
      predictX, false, "frequency", true);
    ampPrediction = predict_curve(ampFitX', gappedPartials{n}(ampFitX, 2), ...
      predictX, true, "amplitude", true);
    gappedPartials{n}(predictX, 1) = freqPrediction;
    gappedPartials{n}(predictX, 2) = max(ampPrediction, gappedPartials{n}(predictX, 2));

    % build the phase
    deltaPhase = ceil(round((N + mod(N, 2))/HopFactor)) .* movmean(gappedPartials{n}(:, 1), 2); % quadratic polynomial, use movmean to calculate (f[k+1] + f[k])/2
    wrappedDeltaPhase = wrapToPi(deltaPhase);

    nanArray = isnan(gappedPartials{n}(:, 3));
    firstNanIndex = find(nanArray(firstGapFrameIdx-extraFrameNum:lastGapFrameIdx+extraFrameNum), 1, 'first') + firstGapFrameIdx - extraFrameNum - 1;
    lastNanIndex = find(nanArray(firstGapFrameIdx-extraFrameNum:lastGapFrameIdx+extraFrameNum), 1, 'last') + firstGapFrameIdx - extraFrameNum - 1;
    
    if isempty(firstNanIndex) && isempty(lastNanIndex)
      fprintf('Skip build the phase for partial %d\n', n);
      continue;
    end
    while nanArray(firstNanIndex - 1)
      firstNanIndex = firstNanIndex - 1;
    end
    while nanArray(lastNanIndex + 1)
      lastNanIndex = lastNanIndex - 1;
    end
    for k = firstNanIndex:lastNanIndex
      gappedPartials{n}(k, 3) = wrapToPi(gappedPartials{n}(k - 1, 3)+deltaPhase(k));
    end
    phaseEstimation = wrapToPi(gappedPartials{n}(lastNanIndex, 3)+deltaPhase(lastNanIndex+1));
    phaseError = wrapToPi(gappedPartials{n}(lastNanIndex + 1, 3)-phaseEstimation);
    phaseSpreadError = phaseError .* linspace(0, 1, lastNanIndex-firstNanIndex+3)';
    gappedPartials{n}(firstNanIndex:lastNanIndex, 3) = wrapToPi(gappedPartials{n}(firstNanIndex:lastNanIndex, 3)+phaseSpreadError(2:end-1));
  elseif (mergedPartialType(n) == 2) % unmatched birth partials
    % predict frequency based on only right part
    gappedPartials{n}(predictX, 1) = rightFreqPrediction;

    if length(rightAmpTrend) >= 2
      trendSlope = rightAmpTrend(2) - rightAmpTrend(1);
      trendSlope = max(trendSlope, -maxDeltaAmp / 1.5);
    else
      trendSlope = 0;
    end

    if trendSlope < 0 % parabola shape to simulate a short attack
      parabolaCoef = [(lastGapFrameIdx + extraFrameNum)^2, lastGapFrameIdx + extraFrameNum, 1; ...
        2 * (lastGapFrameIdx + extraFrameNum), 1, 0; ...
        (lastGapFrameIdx + extraFrameNum - maxBirthLen)^2, lastGapFrameIdx + extraFrameNum - maxBirthLen, 1] \ ...
        [rightAmpPrediction(end); max(trendSlope, -maxDeltaAmp); minAmpThreshold];
      attackCurve = parabolaCoef(1) * predictX.^2 + parabolaCoef(2) * predictX + parabolaCoef(3);
      gappedPartials{n}(predictX, 2) = attackCurve;
    else % keep increasing from threshold in the gap
      if rightAmpPrediction(end-maxBirthLen+1) < minAmpThreshold
        gappedPartials{n}(predictX(end-maxBirthLen+1:end), 2) = rightAmpPrediction(end-maxBirthLen+1:end);
      else
        fadeCurve = linspace(1, 0, maxBirthLen+1)' .^ 2 .* (minAmpThreshold - rightAmpPrediction(end-maxBirthLen+1));
        gappedPartials{n}(predictX(end-maxBirthLen+1:end), 2) = rightAmpPrediction(end-maxBirthLen+1:end) + fadeCurve(1:end-1);
      end
    end

    % build the phase
    deltaPhase = ceil(round((N + mod(N, 2))/HopFactor)) .* movmean(gappedPartials{n}(:, 1), 2); % quadratic polynomial, use movmean to calculate (f[k+1] + f[k])/2
    nanArray = isnan(gappedPartials{n}(:, 3));
    firstNanIndex = find(nanArray(firstGapFrameIdx-extraFrameNum:lastGapFrameIdx+extraFrameNum), 1, 'first') + firstGapFrameIdx - extraFrameNum - 1;
    lastNanIndex = find(nanArray(firstGapFrameIdx-extraFrameNum:lastGapFrameIdx+extraFrameNum), 1, 'last') + firstGapFrameIdx - extraFrameNum - 1;
    if isempty(firstNanIndex) && isempty(lastNanIndex)
      continue;
    end
    while nanArray(lastNanIndex + 1)
      lastNanIndex = lastNanIndex - 1;
    end
    for k = lastNanIndex:-1:firstNanIndex
      gappedPartials{n}(k, 3) = wrapToPi(gappedPartials{n}(k + 1, 3)-deltaPhase(k));
    end
  elseif (mergedPartialType(n) == 3) % unmatched death partials
    % predict frequency based on only left part
    gappedPartials{n}(predictX, 1) = leftFreqPrediction;

    if length(leftAmpTrend) >= 2
      trendSlope = leftAmpTrend(2) - leftAmpTrend(1);
    else
      trendSlope = 0;
    end

    if trendSlope > 0 % parabola shape for the attack and decay
      parabolaCoef = [(firstGapFrameIdx - extraFrameNum)^2, firstGapFrameIdx - extraFrameNum, 1; ...
        2 * (firstGapFrameIdx - extraFrameNum), 1, 0; ...
        lastGapFrameIdx^2, lastGapFrameIdx, 1] \ ...
        [leftAmpPrediction(1); min(trendSlope, maxDeltaAmp); minAmpThreshold];
      decayCurve = parabolaCoef(1) * predictX.^2 + parabolaCoef(2) * predictX + parabolaCoef(3);
      gappedPartials{n}(predictX, 2) = decayCurve;
    else % keep decreasing to the threshold in the gap
      if leftAmpPrediction(maxDeathLen) < minAmpThreshold
        gappedPartials{n}(predictX(1:maxDeathLen), 2) = leftAmpPrediction(1:maxDeathLen);
      else
        fadeCurve = linspace(0, 1, maxDeathLen+1)' .^ 2 .* (minAmpThreshold - leftAmpPrediction(maxDeathLen));
        gappedPartials{n}(predictX(1:maxDeathLen), 2) = leftAmpPrediction(1:maxDeathLen) + fadeCurve(2:end);
      end
    end

    % build the phase
    deltaPhase = ceil(round((N + mod(N, 2))/HopFactor)) .* movmean(gappedPartials{n}(:, 1), 2); % quadratic polynomial, use movmean to calculate (f[k+1] + f[k])/2
    nanArray = isnan(gappedPartials{n}(:, 3));
    firstNanIndex = find(nanArray(firstGapFrameIdx-extraFrameNum:lastGapFrameIdx+extraFrameNum), 1, 'first') + firstGapFrameIdx - extraFrameNum - 1;
    lastNanIndex = find(nanArray(firstGapFrameIdx-extraFrameNum:lastGapFrameIdx+extraFrameNum), 1, 'last') + firstGapFrameIdx - extraFrameNum - 1;
    if isempty(firstNanIndex) && isempty(lastNanIndex)
      continue;
    end
    while nanArray(firstNanIndex - 1)
      firstNanIndex = firstNanIndex - 1;
    end
    for k = firstNanIndex:lastNanIndex
      gappedPartials{n}(k, 3) = wrapToPi(gappedPartials{n}(k - 1, 3)+deltaPhase(k));
    end
  end
end

for k = firstGapFrameIdx - extraFrameNum:lastGapFrameIdx + extraFrameNum
  tempPartials = zeros(length(gappedPartials), 4);
  for n = 1:length(gappedPartials)
    tempPartials(n, 1:3) = gappedPartials{n}(k, :);
    tempPartials(n, 4) = mergedPartialIdx(n);
  end
  matchedPartials{k} = tempPartials;
end

% remove rows with NaN
for n = 1:length(matchedPartials)
  toRemoveRows = sum(isnan(matchedPartials{n}), 2) ~= 0;
  matchedPartials{n}(toRemoveRows, :) = [];
end

predictedPartials = matchedPartials;

partialSig = jun_synthesize_partials(predictedPartials, timeLoc, sigLen+sum(paddingLen));
% Remove zero-padding
partialSig = partialSig(paddingLen(1)+1:end-paddingLen(2));
end