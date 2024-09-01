function [reconnectedPartials, partialsInOrder] = partial_reconnection(partials, partialsInOrder, minPartialConnectLen, totalFrameNum, maxTimeJump, maxOverlapNum, freqThreshold, ampThreshold, freqOrAmp, fs)
partialLens = cellfun(@(p) (size(p, 1)), partialsInOrder);
[partialLenSort, partialIdx] = sort(partialLens, 'descend');
partialIdx = partialIdx(partialLenSort >= minPartialConnectLen); % ignore short partials for prediction and connection

% ignore the phase in the prediction part
for n = 1:length(partialIdx)
  p = partialIdx(n);
  if isempty(partialsInOrder{p})
    continue;
  end

  while (true)
    tempPartialEst = zeros(totalFrameNum, 2); % freq, amp
    fitX = partialsInOrder{p}(:, 4);
    predictX = 1:totalFrameNum;
    tempPartialEst(:, 1) = predict_curve(fitX, partialsInOrder{p}(:, 1), predictX, false, "frequency", false);
    tempPartialEst(:, 2) = predict_curve(fitX, partialsInOrder{p}(:, 2), predictX, false, "amplitude", false);

    freqDist = Inf(size(partialsInOrder));
    ampDist = Inf(size(partialsInOrder));
    p1LeftIdx = partialsInOrder{p}(1, 4);
    p1RightIdx = partialsInOrder{p}(end, 4);
    for p2 = 1:length(partialsInOrder)
      if isempty(partialsInOrder{p2}) || p2 == p
        % ignore empty partials
        continue;
      end
      p2LeftIdx = partialsInOrder{p2}(1, 4);
      p2RightIdx = partialsInOrder{p2}(end, 4);
      minIdx = min(p1LeftIdx, p2LeftIdx);
      maxIdx = max(p1RightIdx, p2RightIdx);
      overlapNum = (p1RightIdx - p1LeftIdx + 1 + p2RightIdx - p2LeftIdx + 1) - (maxIdx - minIdx + 1);
      % if two partials are too far in time or in the same time range, ignore this pair
      if (~((p1LeftIdx >= p2LeftIdx - maxTimeJump - 1 && p1LeftIdx <= p2RightIdx + maxTimeJump + 1) || ...
          (p1RightIdx >= p2LeftIdx - maxTimeJump - 1 && p1RightIdx <= p2RightIdx + maxTimeJump + 1)) || ...
          (overlapNum >= min(maxOverlapNum+1, min(p1RightIdx-p1LeftIdx+1, p2RightIdx-p2LeftIdx+1))))
        continue;
      end
      % if the boundary frequencies of two partials are larger than 6%, ignore this pair
      if maxIdx == p1RightIdx
        p1RefFreq = partialsInOrder{p}(1 + floor(max(overlapNum, 0)/2), 1);
        p2RefFreq = partialsInOrder{p2}(end -floor(max(overlapNum, 0)/2), 1);
      elseif maxIdx == p2RightIdx
        p1RefFreq = partialsInOrder{p}(end -floor(max(overlapNum, 0)/2), 1);
        p2RefFreq = partialsInOrder{p2}(1 + floor(max(overlapNum, 0)/2), 1);
      end
      if abs(freq_to_erb(p1RefFreq/2/pi*fs)-freq_to_erb(p2RefFreq/2/pi*fs)) > 0.4
        continue;
      end
      if (overlapNum >= 1)
        [fitX, nonOverlapIdx] = setdiff(partialsInOrder{p}(:, 4), partialsInOrder{p2}(:, 4));
        tempPartialEst(:, 1) = predict_curve(fitX, partialsInOrder{p}(nonOverlapIdx, 1), predictX, false, "frequency", false);
        tempPartialEst(:, 2) = predict_curve(fitX, partialsInOrder{p}(nonOverlapIdx, 2), predictX, false, "amplitude", false);
      end
      freqEucliDist = norm((freq_to_erb(tempPartialEst(p2LeftIdx:p2RightIdx, 1)) - freq_to_erb(partialsInOrder{p2}(:, 1)))) / ...
        sqrt(p2RightIdx-p2LeftIdx+1);
      freqDist(p2) = freqEucliDist / (1 + std(freq_to_erb(tempPartialEst(p2LeftIdx:p2RightIdx, 1))) + std(freq_to_erb(partialsInOrder{p2}(:, 1))));
      ampEucliDist = norm(tempPartialEst(p2LeftIdx:p2RightIdx, 2)-partialsInOrder{p2}(:, 2)) / sqrt(p2RightIdx-p2LeftIdx+1);
      ampDist(p2) = ampEucliDist / (1 + std(tempPartialEst(p2LeftIdx:p2RightIdx, 2)) + std(partialsInOrder{p2}(:, 2)));
    end

    freqDist = freqDist ./ freqThreshold;
    freqDist(freqDist > 1) = Inf;
    ampDist = ampDist ./ ampThreshold;
    ampDist(ampDist > 1) = Inf;
    totalDist = freqDist * freqOrAmp + ampDist * (1 - freqOrAmp);
    p2 = find(totalDist == min(totalDist) & ~isinf(totalDist), 1);

    if ~isempty(p2)
      % merge two partials
      p2LeftIdx = partialsInOrder{p2}(1, 4);
      p2RightIdx = partialsInOrder{p2}(end, 4);
      minIdx = min(p1LeftIdx, p2LeftIdx);
      maxIdx = max(p1RightIdx, p2RightIdx);
      overlapNum = (p1RightIdx - p1LeftIdx + 1 + p2RightIdx - p2LeftIdx + 1) - (maxIdx - minIdx + 1);

      mergePartial = NaN(maxIdx-minIdx+1, 4);
      mergePartial(:, 4) = minIdx:maxIdx;
      if (overlapNum > 0)
        % has overlaps
        rampUpCurve = rampup(overlapNum+1);
        rampUpCurve = rampUpCurve(2:end);
        rampDownCurve = rampdown(overlapNum+1);
        rampDownCurve = rampDownCurve(2:end);
        if (p1LeftIdx == minIdx) % p1: ramp down, p2: ramp up
          mergePartial(1:p2LeftIdx-minIdx, 1:3) = partialsInOrder{p}(1:p2LeftIdx - minIdx, 1:3);
          mergePartial(p2LeftIdx-minIdx+1:p1RightIdx-minIdx+1, 1:2) = ...
            partialsInOrder{p}(end - overlapNum + 1:end, 1:2) .* rampDownCurve + ...
            partialsInOrder{p2}(1:overlapNum, 1:2) .* rampUpCurve;
          mergePartial(p1RightIdx-minIdx+2:end, 1:3) = partialsInOrder{p2}(overlapNum + 1:end, 1:3);
        else % p1: ramp up, p2: ramp down
          mergePartial(1:p1LeftIdx-minIdx, 1:3) = partialsInOrder{p2}(1:p1LeftIdx - minIdx, 1:3);
          mergePartial(p1LeftIdx-minIdx+1:p2RightIdx-minIdx+1, 1:2) = ...
            partialsInOrder{p2}(end - overlapNum + 1:end, 1:2) .* rampDownCurve + ...
            partialsInOrder{p}(1:overlapNum, 1:2) .* rampUpCurve;
          mergePartial(p2RightIdx-minIdx+2:end, 1:3) = partialsInOrder{p}(overlapNum + 1:end, 1:3);
        end
      else
        % no overlap
        mergePartial((p1LeftIdx:p1RightIdx)-minIdx+1, 1:3) = partialsInOrder{p}(:, 1:3);
        mergePartial((p2LeftIdx:p2RightIdx)-minIdx+1, 1:3) = partialsInOrder{p2}(:, 1:3);
        
      end
      partialsInOrder{p} = mergePartial;
      partialsInOrder{p2} = double.empty(0, 4);
    else
      break;
    end
  end
end

partialsInOrder = partialsInOrder(cellfun(@(x) ~isempty(x), partialsInOrder));

% get new partials
reconnectedPartials = cellfun(@(~)double.empty(0, 4), cell(size(partials)), 'UniformOutput', false);
for p = 1:length(partialsInOrder)
  for n = 1:size(partialsInOrder{p}, 1)
    bin = partialsInOrder{p}(n, 4);
    reconnectedPartials{bin} = [reconnectedPartials{bin}; [partialsInOrder{p}(n, 1:3), p]];
  end
end
end