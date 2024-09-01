function [reconnectedPartials, leftPartialIdx, rightPartialIdx] = partial_matching(reconnectedPartials, firstGapFrameIdx, lastGapFrameIdx, extraFrameNum, totalFrameNum, minPartialLen, minSlopeRatio, maxDeltaFreq, maxDeltaAmp, matchPreference, fs)
leftPartialIdx = reconnectedPartials{firstGapFrameIdx-1}(:, 4);
rightPartialIdx = reconnectedPartials{lastGapFrameIdx+1}(:, 4);
leftPartialOffset = zeros(size(leftPartialIdx));
rightPartialOffset = zeros(size(rightPartialIdx));
for n = 1:extraFrameNum
  tempLeftIdx = setdiff(reconnectedPartials{firstGapFrameIdx-1-n}(:, 4), leftPartialIdx);
  tempRightIdx = setdiff(reconnectedPartials{lastGapFrameIdx+1+n}(:, 4), rightPartialIdx);
  leftPartialIdx = [leftPartialIdx; tempLeftIdx];
  rightPartialIdx = [rightPartialIdx; tempRightIdx];
  leftPartialOffset = [leftPartialOffset; n * ones(size(tempLeftIdx))];
  rightPartialOffset = [rightPartialOffset; n * ones(size(tempRightIdx))];
end

leftPartials = cellfun(@(~)double.empty(0, 3), cell(length(leftPartialIdx), 1), 'UniformOutput', false);
rightPartials = cellfun(@(~)double.empty(0, 3), cell(length(rightPartialIdx), 1), 'UniformOutput', false);

for idx = 1:firstGapFrameIdx - 1
  for p = 1:size(reconnectedPartials{idx}, 1)
    tempIdx = find(leftPartialIdx == reconnectedPartials{idx}(p, 4), 1);
    if ~isempty(tempIdx)
      leftPartials{tempIdx} = [leftPartials{tempIdx}; reconnectedPartials{idx}(p, 1:3)];
    end
  end
end
for n = find(leftPartialOffset ~= 0, 1):length(leftPartials)
  leftPartials{n} = [leftPartials{n}; NaN(leftPartialOffset(n), 3)];
end
for idx = lastGapFrameIdx + 1:totalFrameNum
  for p = 1:size(reconnectedPartials{idx}, 1)
    tempIdx = find(rightPartialIdx == reconnectedPartials{idx}(p, 4), 1);
    if ~isempty(tempIdx)
      rightPartials{tempIdx} = [rightPartials{tempIdx}; reconnectedPartials{idx}(p, 1:3)];
    end
  end
end
for n = find(rightPartialOffset ~= 0, 1):length(rightPartials)
  rightPartials{n} = [NaN(rightPartialOffset(n), 3); rightPartials{n}];
end

leftPartialLen = cellfun(@(x)(sum(~isnan(x(:, 1)))), leftPartials);
leftPartials = leftPartials(leftPartialLen >= minPartialLen);
rightPartialLen = cellfun(@(x)(sum(~isnan(x(:, 1)))), rightPartials);
rightPartials = rightPartials(rightPartialLen >= minPartialLen);
discardPartialIdx = [leftPartialIdx(leftPartialLen < minPartialLen); rightPartialIdx(rightPartialLen < minPartialLen)];
leftPartialIdx = leftPartialIdx(leftPartialLen >= minPartialLen);
rightPartialIdx = rightPartialIdx(rightPartialLen >= minPartialLen);
% remove useless partials
for n = [(-minPartialLen + 1:-1) + firstGapFrameIdx, (1:minPartialLen - 1) + lastGapFrameIdx]
  reconnectedPartials{n} = reconnectedPartials{n}(~ismember(reconnectedPartials{n}(:, 4), discardPartialIdx), :);
end

predictLen = lastGapFrameIdx - firstGapFrameIdx + 1;
predictType = ["frequency", "amplitude"];
leftPartialPrediction = zeros(length(leftPartials), predictLen, 2); % (:,:,1)
for n = 1:length(leftPartials)
  % linear regression for the trend, AR for the sinusoids
  [changePoint, slope, intercept] = ischange(leftPartials{n}(:, 2)', 'linear', MaxNumChanges = 1);
  changePointIdx = find(changePoint, 1);
  if ~isempty(changePointIdx) % might include an attack
    attackSlope = polyfit(1:changePointIdx, leftPartials{n}(1:changePointIdx, 2), 1);
    attackSlope = attackSlope(1);
    decaySlope = polyfit(1:size(leftPartials{n}, 1)-changePointIdx+1, leftPartials{n}(changePointIdx:end, 2), 1);
    decaySlope = decaySlope(1);
    if attackSlope > -decaySlope * minSlopeRatio && decaySlope < 0
      fitX = (-size(leftPartials{n}, 1) + changePointIdx - 1:-1)' + firstGapFrameIdx;
    else
      fitX = (-size(leftPartials{n}, 1):-1)' + firstGapFrameIdx;
    end
  else
    fitX = (-size(leftPartials{n}, 1):-1)' + firstGapFrameIdx;
  end

  if length(fitX) < 2
    continue;
  end

  for k = 1:2 % 1 for frequency, 2 for amplitude
    if k == 1
      tempIdx = (-length(fitX) + 1:0) + size(leftPartials{n}, 1);
    else
      ignoreIdx = min(length(fitX), extraFrameNum+1) - 1;
      tempIdx = (-length(fitX) + 1:-ignoreIdx) + size(leftPartials{n}, 1);
      fitX = fitX(1:end-ignoreIdx);
    end
    tempFill = predict_curve(fitX, leftPartials{n}(tempIdx, k), firstGapFrameIdx:lastGapFrameIdx, true, predictType(k), false);
    leftPartialPrediction(n, :, k) = tempFill;
  end
end
rightPartialPrediction = zeros(length(rightPartials), predictLen, 2);
for n = 1:length(rightPartials)
  % linear regression for the trend, AR for the sinusoids
  fitX = (1:size(rightPartials{n}, 1))' + lastGapFrameIdx;
  if length(fitX) < 2
    continue;
  end

  for k = 1:2 % 1 for frequency, 2 for amplitude
    if k == 1
      tempIdx = 1:size(rightPartials{n}, 1);
    else
      ignoreIdx = min(length(fitX), extraFrameNum+1) - 1;
      tempIdx = ignoreIdx + 1:size(rightPartials{n}, 1);
      fitX = fitX(ignoreIdx+1:end);
    end
    tempFill = predict_curve(fitX, rightPartials{n}(tempIdx, k), firstGapFrameIdx:lastGapFrameIdx, true, predictType(k), false);
    rightPartialPrediction(n, :, k) = tempFill;
  end
end

freqEucliDistMat = zeros(size(leftPartialPrediction, 1), size(rightPartialPrediction, 1));
ampEucliDistMat = zeros(size(leftPartialPrediction, 1), size(rightPartialPrediction, 1));
for row = 1:size(freqEucliDistMat, 1)
  for col = 1:size(freqEucliDistMat, 2)
    freqDiff = freq_to_erb(leftPartialPrediction(row, :, 1)) - freq_to_erb(rightPartialPrediction(col, :, 1));
    freqEucliDistMat(row, col) = norm(freqDiff) / sqrt(predictLen) / ...
      (1 + std(freq_to_erb(leftPartialPrediction(row, :, 1))) + std(freq_to_erb(rightPartialPrediction(col, :, 1))));
    ampDiff = leftPartialPrediction(row, :, 2) - rightPartialPrediction(col, :, 2);
    ampEucliDistMat(row, col) = norm(ampDiff) / sqrt(predictLen) / ...
      (1 + std(leftPartialPrediction(row, :, 2)) + std(rightPartialPrediction(col, :, 2)));
  end
end

varFreq = -maxDeltaFreq^2 * log((matchPreference - 1)/(matchPreference - 2));
varAmp = -maxDeltaAmp^2 * log((matchPreference - 1)/(matchPreference - 2));

usefulCostMat = 1 - exp(-freqEucliDistMat.^2/varFreq-ampEucliDistMat.^2/varAmp);
spuriousCostMat = 1 - (1 - matchPreference) * usefulCostMat;
multiCostMat = min(usefulCostMat, spuriousCostMat);
costTypeMat = usefulCostMat == multiCostMat; % true: a useful trajectory, false: a spurious assignment

assignment = munkres(multiCostMat);

tempPartialIdx = NaN(length(assignment), 1);
for idx = lastGapFrameIdx + 1:totalFrameNum
  partialIdx = reconnectedPartials{idx}(:, 4);
  for n = 1:length(assignment)
    if (assignment(n) ~= 0 && costTypeMat(n, assignment(n)))
      partialIdx(partialIdx == rightPartialIdx(assignment(n))) = leftPartialIdx(n);
      tempPartialIdx(n) = leftPartialIdx(n);
    end
  end
  reconnectedPartials{idx}(:, 4) = partialIdx;
end

rightPartialIdx = [];
for n = lastGapFrameIdx + 1:lastGapFrameIdx + 1 + extraFrameNum
  rightPartialIdx = union(rightPartialIdx, reconnectedPartials{n}(:, 4));
end
end