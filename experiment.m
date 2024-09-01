%% Load data
clc;
clearvars;
rng(7, "twister"); % set the random number generator to the default seed

addpath('HybridMethod');
addpath('SPAIN');
addpath('Janssen');
addpath('SocialSparsity');
addpath('ReweightedInpainting');
addpath(genpath('FastPartialTracking'));
addpath(genpath('PEMO-Q'));

sigName = {
  "09_Viola", ...
  "23_Horn", ...
  "35_Glockenspiel", ...
  "44_Soprano", ...
  "56_Organ", ...
  "58_Guitar", ...
};

methodName = ["hybrid", "A-SPAIN", "w-CP", "re-CP", "f-Janssen"];

warning('off');

gapLenInMs = [10, 25, 50, 100, 150]'; % in ms

progressbar = waitbar(0, 'Starting...');

for gapId = 1:length(gapLenInMs)
  sigMetrics = cell(10, 1);

  for sigId = 1:length(sigName)
    fprintf("Gap length = %d ms, signal name = %s\n", gapLenInMs(gapId), sigName{sigId});

    [data, fs] = audioread(strcat('./signals/', sigName{sigId}, '.wav'));
    data = data ./ max(abs(data));
    
    gapLen = round(gapLenInMs(gapId) / 1000 * fs);
    if gapLenInMs(gapId) <= 200
      segmentNum = 10; % number of segments
      segments = create_segments(data, segmentNum);
      segments = create_gaps(segments, gapLen, 0.30*fs);
    else
      segmentNum = 6;
      segments = create_segments(data, segmentNum);
      segments = create_gaps(segments, gapLen, 0.53*fs);
    end
    
    time.hybrid = zeros(segmentNum-2, 1);
    snr.hybrid = zeros(segmentNum-2, 1);
    odg.hybrid = zeros(segmentNum-2, 1);
    time.aspain = zeros(segmentNum-2, 1);
    snr.aspain = zeros(segmentNum-2, 1);
    odg.aspain = zeros(segmentNum-2, 1);
    time.wcp = zeros(segmentNum-2, 1);
    snr.wcp = zeros(segmentNum-2, 1);
    odg.wcp = zeros(segmentNum-2, 1);
    time.recp = zeros(segmentNum-2, 1);
    snr.recp = zeros(segmentNum-2, 1);
    odg.recp = zeros(segmentNum-2, 1);
    time.fjanssen = zeros(segmentNum-2, 1);
    snr.fjanssen = zeros(segmentNum-2, 1);
    odg.fjanssen = zeros(segmentNum-2, 1);

    for segId = 2:segmentNum-1
      waitbar((segId-1 + (segmentNum-2)*(sigId-1) + (segmentNum-2)*length(sigName)*(gapId-1)) / ...
        ((segmentNum-2)*length(sigName)*length(gapLenInMs)), ...
        progressbar, ...
        sprintf('Gap %d, Sig %d, Seg %d', gapId, sigId, segId));
      reliableMask = segments{segId}.mask;
      gapStartIdx = find(~reliableMask, 1, 'first');
      gapEndIdx = find(~reliableMask, 1, 'last');
      clipSig = segments{segId}.gapSig;
      
      if ismember("hybrid", methodName)
        t = tic;
        %% Define frame parameters
        w = [1024 / 1, 64]; % window size
        M = w; % channel
        a = [w(1) / 8, w(2) / 4]; % time shift

        g = {};
        g{1} = gabwin({'tight', 'hann'}, a(1), w(1), w(1));
        g{2} = gabwin({'tight', 'hann'}, a(2), w(2), w(2));

        Ftonal = frametight(frame('dgtreal', g{1}, a(1), M(1), 'freqinv'));
        Ftrans = frametight(frame('dgtreal', g{2}, a(2), M(2)));

        offsetMode = 'half';
        switch offsetMode
          case 'full'
            c = ceil((gapStartIdx + gapEndIdx)/2);
            k = floor((c - 1)/a(1));
            d = 1 + k * a(1);
            offset = c - d;
          case 'half'
            c = ceil((gapStartIdx + gapEndIdx)/2);
            k = floor((c - 1)/a(1));
            d = 1 + k * a(1) + ceil(a(1)/2);
            offset = c - d;
          otherwise
            offset = 0;
        end
        [leftBound, rightBound, ~, ~, ~, ~, ~, ~, ~, ~, ~] = min_sig_supp_2( ...
          w(1), a(1), M(1), gapStartIdx, gapEndIdx, length(clipSig), 1, offset);
        if mod(rightBound-leftBound, 2) == 0
          rightBound = rightBound - 1;
        end
        boundExpansionSize = (ceil(gapLen/a(1)) + 50) * a(1);
        leftBound = max(leftBound-boundExpansionSize, 1);
        rightBound = min(rightBound+boundExpansionSize, length(clipSig));
        clipSig = clipSig(leftBound:rightBound);
        reliableMask = reliableMask(leftBound:rightBound);
        gapStartIdx = find(~reliableMask, 1, 'first');
        gapEndIdx = find(~reliableMask, 1, 'last');

        normalizedFactor = 1 / max(abs(clipSig));
        clipSig = clipSig .* normalizedFactor;

        %% multilayered expansions
        optimalLambda = get_optimal_lambda(clipSig, reliableMask, 1);
        if exist('optimalLambda', 'var')
          lambda = max(optimalLambda, 2e-5);
        end

        [tonalSig, tonalResidual] = social_analysis(Ftonal, clipSig, ...
          optimalLambda, reliableMask, struct('shrinkage', 'PEW'), ...
          'tolerence', 0.75);
        tonalCoef = frana(Ftonal, tonalSig);

        tonalSig = tonalSig(1:length(clipSig)); % align with the shape of the mask
        tonalResidual = tonalResidual(1:length(clipSig));
        tonalResidual(~reliableMask) = 0;

        %% Fast Partial tracking: initialization
        % Analysis Window Length
        if gapLenInMs(gapId) <= 25
          N = 1024 - 1;
        else
          N = 2048 - 1;
        end
        % Oversampling Factor
        OverSample = 1;
        % Hop Size Factor (HopSize = N/HopFactor)
        HopFactor = 8;
        % Magnitude Threshold for Peak-Picking (dB)
        Peak_dB = -50;
        % Polynomial Order Q of the short-term sinusoidal model
        % 1 = Frequency, Damping
        % 2 = Frequency Derivative, Damping Derivative
        % 3 = Frequency 2nd Derivative, Damping 2nd Derivative, and so on..
        Q = 3;
        % These parameters control the assignment cost function
        delta = .5; % 0 < delta < 1
        zeta_f = 0.2;
        zeta_a = 6; % dB
        zeta_fs = 0.005;

        [partials, timeLoc, paddingLen, sigLen, stftCoef, reliableFrameMask, extraFrameNum] = jun_track_partials( ...
          tonalSig, fs, N, HopFactor, OverSample, Peak_dB, Q, delta, zeta_f, zeta_a, zeta_fs, reliableMask);
        [partials, timeLoc, paddingLen, sigLen, ~, num_tracks] = jun_read_partials('demo_partials.bin');
        % jun_plot_partials(partials, timeLoc-paddingLen(1), fs, num_tracks);

        partialsInOrder = cellfun(@(~)double.empty(0, 4), cell(num_tracks, 1), 'UniformOutput', false); % freq, amp, phase, time bin
        for n = 1:length(partials)
          for p = 1:size(partials{n}, 1)
            pIdx = partials{n}(p, 4);
            partialsInOrder{pIdx} = [partialsInOrder{pIdx}; [partials{n}(p, 1:3), n]];
          end
        end

        %% re-connect fragmented partials
        minPartialConnectLen = 4;
        maxTimeJump = 4;
        maxOverlapNum = 3;
        freqThreshold = 3e-4;
        ampThreshold = 5;
        ampThreshold = ampThreshold / 20 * log(10);
        freqOrAmp = 0.5; % prediction tendency: a freq, 1-a amp

        totalFrameNum = length(partials);

        [reconnectedPartials, partialsInOrder] = partial_reconnection(...
          partials, partialsInOrder, minPartialConnectLen, totalFrameNum, ...
          maxTimeJump, maxOverlapNum, freqThreshold, ampThreshold, freqOrAmp, fs);
        % jun_plot_partials(reconnectedPartials, timeLoc-paddingLen(1), fs, num_tracks);

        %% Partial matching
        minSlopeRatio = 2.5; % detect if it's a damped partial with attack, the slope of "attack" part should be at least 2 times than the decay part
        minPartialLen = 4;
        maxDeltaFreq = 0.0003;
        maxDeltaAmp = 10;
        maxDeltaAmp = maxDeltaAmp / 20 * log(10);
        matchPreference = 0.5;

        firstGapFrameIdx = find(~reliableFrameMask, 1, 'first');
        lastGapFrameIdx = find(~reliableFrameMask, 1, 'last');

        if ~(isempty(firstGapFrameIdx) || isempty(lastGapFrameIdx))
          [matchedPartials, leftPartialIdx, rightPartialIdx] = partial_matching(...
            reconnectedPartials, firstGapFrameIdx, lastGapFrameIdx, extraFrameNum, ...
            totalFrameNum, minPartialLen, minSlopeRatio, maxDeltaFreq, maxDeltaAmp, ...
            matchPreference, fs);
        end

        %% partial prediction
        maxBirthLen = lastGapFrameIdx - firstGapFrameIdx + 1 + extraFrameNum * 2 - 1; % the maximum length of the extrapolation of partial birth
        maxDeathLen = lastGapFrameIdx - firstGapFrameIdx + 1 + extraFrameNum * 2 - 1; % the maximum length of the extrapolation of partial death
        minAmpThreshold = log(db2mag(-80)); % -80 dB
        maxDeltaAmp = -minAmpThreshold / (lastGapFrameIdx - firstGapFrameIdx + 1 + extraFrameNum * 2 - 1);

        if ~(isempty(firstGapFrameIdx) || isempty(lastGapFrameIdx))
          [predictedPartials, partialSig] = partial_prediction(...
            leftPartialIdx, rightPartialIdx, totalFrameNum, matchedPartials, ...
            firstGapFrameIdx, extraFrameNum, lastGapFrameIdx, fs, N, HopFactor, ...
            minSlopeRatio, maxBirthLen, maxDeltaAmp, minAmpThreshold, maxDeathLen, ...
            timeLoc, sigLen, paddingLen);
        else
          partialSig = jun_synthesize_partials(partials, ...
            timeLoc, sigLen + sum(paddingLen));
          partialSig = partialSig(paddingLen(1)+1:end-paddingLen(2));
        end

        %% Shorten tonal residual signal further for transient decomposition
        shortNeighborSize = max(ceil(gapLen/4), w(1)) * 2;
        shortenedLeftIdx = gapStartIdx - shortNeighborSize;
        shortenedRightIdx = gapEndIdx + shortNeighborSize;
        shortenedGapStartIdx = shortNeighborSize + 1;
        shortenedGapEndIdx = shortenedGapStartIdx + gapLen - 1;

        shortClipSig = clipSig(shortenedLeftIdx:shortenedRightIdx);
        shortTonalResidual = tonalResidual(shortenedLeftIdx:shortenedRightIdx);
        shortReliableMask = reliableMask(shortenedLeftIdx:shortenedRightIdx);

        %% Transient decomposition
        optimalLambda = get_optimal_lambda(shortTonalResidual, shortReliableMask, 2);

        [transSig, transResidual] = social_analysis(Ftrans, shortTonalResidual, optimalLambda, shortReliableMask, ...
          struct('shrinkage', 'PEW', 'weighting', false, 'reweighting', false, 'neighbor', ones(25, 1), 'center', [13, 1]), 'tolerence', 0.25, 'transient');
        transSig = transSig(1:length(shortClipSig)); % align with the shape of the mask
        transSig(isnan(transSig)) = 0;
        transResidual(isnan(transResidual)) = 0;
        transResidual(~shortReliableMask) = 0;

        %% noise inpainting
        % ignore "reliable" signals near the gap
        noiseRecMethod = 1;
        ignoreLen = ceil(w(1) / 4) * 2;
        % ignoreLen = 0;
        noiseWindowSize = 1024;
        noiseReliableMask = true(size(reliableMask));
        noiseReliableMask(gapStartIdx-ignoreLen:gapEndIdx+ignoreLen) = false;
        noiseReliableMask = noiseReliableMask(shortenedLeftIdx:shortenedRightIdx);

        shortNoiseSig = reconstruct_noise(transResidual, noiseWindowSize, noiseReliableMask, noiseRecMethod);

        %% re-construct the whole restored signal
        recClip = partialSig;
        recClip(shortenedLeftIdx:shortenedRightIdx) = recClip(shortenedLeftIdx:shortenedRightIdx) + shortNoiseSig;
        shortRecClip = recClip(shortenedLeftIdx:shortenedRightIdx);

        wPhase = 1024;
        aPhase = wPhase / 32;
        gPhase = gabwin({'tight', 'hann'}, aPhase, wPhase, wPhase);
        FPhaseLock = frametight(frame('dgtreal', gPhase, aPhase, wPhase, 'timeinv'));
        glaSigLen = framelength(FPhaseLock, length(shortClipSig));
        shortRecClipCoefTf = framecoef2tf(FPhaseLock, frana(FPhaseLock, shortRecClip));
        [shortRecClipCoefPhaseRec, shortRecClipPhaseRec, ~, ~] = gla(shortRecClipCoefTf, gPhase, aPhase, wPhase, 200, 'fgla', 'alpha', 0.99, 'timemod', ...
          @(sig) postpad([shortClipSig(1:shortenedGapStartIdx-1); sig(shortenedGapStartIdx:shortenedGapEndIdx); shortClipSig(shortenedGapEndIdx+1:end)], glaSigLen), 'Ls', length(shortClipSig));
        recClip(shortenedLeftIdx:shortenedRightIdx) = shortRecClipPhaseRec;

        % recClip = partialSig;
        % recClip(shortenedLeftIdx:shortenedRightIdx) = recClip(shortenedLeftIdx:shortenedRightIdx) + shortNoiseSig;

        % cross-fading
        crossFadeLen = round(gapLenInMs(gapId) / 1000 * fs * 0.08);
        % crossFadeLen = 128;
        if crossFadeLen > 0
          crossFadeCurveUp = rampup(crossFadeLen, 'hann');
          crossFadeCurveDown = rampdown(crossFadeLen, 'hann');
          leftFadeIdx = (-crossFadeLen + 1:0) + gapStartIdx - 1;
          rightFadeIdx = (1:crossFadeLen) + gapEndIdx;
          crossFadeOriginSig = clipSig;
          crossFadeRecSig = recClip;
          crossFadeOriginSig(leftFadeIdx) = crossFadeOriginSig(leftFadeIdx) .* crossFadeCurveDown;
          crossFadeOriginSig(~reliableMask) = 0;
          crossFadeOriginSig(rightFadeIdx) = crossFadeOriginSig(rightFadeIdx) .* crossFadeCurveUp;
          crossFadeRecSig(leftFadeIdx) = crossFadeRecSig(leftFadeIdx) .* crossFadeCurveUp;
          crossFadeRecSig(rightFadeIdx) = crossFadeRecSig(rightFadeIdx) .* crossFadeCurveDown;
          crossFadeRecSig(1:leftFadeIdx(1)-1) = 0;
          crossFadeRecSig(rightFadeIdx(end)+1:end) = 0;

          recClip = crossFadeRecSig + crossFadeOriginSig;
        else
          recClip(reliableMask) = clipSig(reliableMask);
        end

        segments{segId}.hybridRecSig = segments{segId}.originSig;
        segments{segId}.hybridRecSig(leftBound:rightBound) = recClip ./ normalizedFactor;

        time.hybrid(segId - 1) = toc(t);
      end

      if ismember("A-SPAIN", methodName)
        %% A-SPAIN
        t = tic;
        red = 4; % redundancy of DFT in each block
        param.w = 2800 + 1400 * (floor(gapLenInMs(gapId) / 28) - 1);
        param.a = param.w / 4; % window shift (300 in original version)
        param.wtype = 'hann';
        param.M = param.w;
        param.offset = 'half'; % 'half' / 'full' / 'none'
        param.mask = reliableMask;
        param.Ls = length(clipSig);

        % DFT parameters
        param.F = frame('dft');
        param.F.redundancy = red; % non-native, our own additional parameter
        param.F.frana = @(insig)dft([insig; zeros(length(insig)*(param.F.redundancy - 1), 1)]);
        param.F.frsyn = @(insig)postpad(idft(insig), length(insig)/param.F.redundancy);

        % paramsolver parameters
        paramsolver.s = 1; % increment of k
        paramsolver.r = 1; % every r-th iteration increment k by s
        paramsolver.epsilon = 0.1; % stopping criterion of termination function
        % paramsolver.maxit = ceil(floor(param.w*param.F.redundancy/2+1)*paramsolver.r/paramsolver.s); % maximum number of iterations
        paramsolver.maxit = 2000;
        paramsolver.store_snr = true;
        paramsolver.store_obj = false;
        paramsolver.f_update = 'h';

        s = gapStartIdx;
        f = gapEndIdx;
        h = f - s + 1;
        switch param.offset
          case 'full'
            c = ceil((s + f)/2);
            k = floor((c - 1)/param.a);
            d = 1 + k * param.a;
            offset = c - d;
          case 'half'
            c = ceil((s + f)/2);
            k = floor((c - 1)/param.a);
            d = 1 + k * param.a + ceil(param.a/2);
            offset = c - d;
          otherwise
            offset = 0;
        end
        [q, Q, ~, ~, ~, ~, ~, ~, ~, ~, param.Ls] = min_sig_supp_2(param.w, param.a, param.M, s, f, param.Ls, 1, offset);

        param.mask = reliableMask(q:Q); % for SPAIN algorithm, where we work with shortened signal
        recClip_SPAIN = clipSig; % ... for A-SPAIN
        data_gapped = clipSig(q:Q);

        param.algorithm = 'aspain';
        [recClip_SPAIN(q:Q), ~, ~] = spain_segmentation(data_gapped, param, paramsolver, clipSig(q:Q));

        segments{segId}.aspainRecSig = segments{segId}.originSig;
        segments{segId}.aspainRecSig(leftBound:rightBound) = recClip_SPAIN ./ normalizedFactor;
        time.aspain(segId - 1) = toc(t);
      end

      if ismember("w-CP", methodName)
        %% weighted Chambolle-Pock
        param.type = 'analysis';  % 'analysis' / 'synthesis'
        param.w = 2800 + 1400 * (floor(gapLenInMs(gapId) / 28) - 1);
        param.F = frametight(frame('dgt', {'hann', param.w}, param.w/4, param.w, 'timeinv'));
        param.offset = 'half';  % 'full' / 'half' / 'none'
        % 'norm' for DR, 'energy' for CP
        param.weighting = 'energy';  % 'none'/ 'supp' / 'abs' / 'norm' / 'energy'
        param.reweighting = false;

        t = tic;
        recClip_wCP = clipSig;
        recClip_wCP(q:Q) = reweighted(clipSig(q:Q), reliableMask(q:Q), param, [], [], []);

        segments{segId}.wcpRecSig = segments{segId}.originSig;
        segments{segId}.wcpRecSig(leftBound:rightBound) = recClip_wCP ./ normalizedFactor;
        time.wcp(segId - 1) = toc(t);
      end

      if ismember("re-CP", methodName)
        %% iteratively reweighted Chambolle-Pock
        param.type = 'analysis';
        param.w = 2800 + 1400 * (floor(gapLenInMs(gapId) / 28) - 1);
        param.F = frametight(frame('dgt', {'hann', param.w}, param.w/4, param.w, 'timeinv'));
        param.offset = 'half';
        param.weighting = 'none';
        param.reweighting = true;

        % parameters of the main cycle
        RWparamsolver.maxit = 10;
        RWparamsolver.epsilon = 1e-3;
        RWparamsolver.delta = 0.01;

        % algorithm
        t = tic;
        recClip_reCP = clipSig;
        recClip_reCP(q:Q) = reweighted(clipSig(q:Q), reliableMask(q:Q), param, [], [], RWparamsolver);

        segments{segId}.recpRecSig = segments{segId}.originSig;
        segments{segId}.recpRecSig(leftBound:rightBound) = recClip_reCP ./ normalizedFactor;
        time.recp(segId - 1) = toc(t);
      end

      if ismember("f-Janssen", methodName)
        %% frame-wise Janssen
        JanssenParam.w = 2800 + 1400 * (floor(gapLenInMs(gapId) / 28) - 1);
        JanssenParam.a = JanssenParam.w / 4;
        JanssenParam.wtype = 'hann';
        JanssenParam.Ls = length(clipSig(q:Q));
        JanssenParam.mask = reliableMask(q:Q);

        % solver settings
        janssenParamSolver.Nit = 20; % number of iterations
        % janssenParamSolver.p = 2*JanssenParam.a;

        % algorithm
        t = tic;
        recoveredFrameJanssen = clipSig;
        [~, recoveredFrameJanssen(q:Q)] = evalc('janssen(clipSig(q:Q), JanssenParam, janssenParamSolver);');

        segments{segId}.fjanssenRecSig = segments{segId}.originSig;
        segments{segId}.fjanssenRecSig(leftBound:rightBound) = recoveredFrameJanssen ./ normalizedFactor;
        time.fjanssen(segId - 1) = toc(t);
      end

      %% calculate SNR, TV-ISD, and ODG
      if ismember("hybrid", methodName)
        snr.hybrid(segId-1) = get_snr(segments{segId}.originSig(~segments{segId}.mask), segments{segId}.hybridRecSig(~segments{segId}.mask));
        [~, ~, ~, odg.hybrid(segId-1), ~] = evalc('audioqual(segments{segId}.originSig, segments{segId}.hybridRecSig, fs);');
      end
      if ismember("A-SPAIN", methodName)
        snr.aspain(segId-1) = get_snr(segments{segId}.originSig(~segments{segId}.mask), segments{segId}.aspainRecSig(~segments{segId}.mask));
        [~, ~, ~, odg.aspain(segId-1), ~] = evalc('audioqual(segments{segId}.originSig, segments{segId}.aspainRecSig, fs);');
      end
      if ismember("w-CP", methodName)
        snr.wcp(segId-1) = get_snr(segments{segId}.originSig(~segments{segId}.mask), segments{segId}.wcpRecSig(~segments{segId}.mask));
        [~, ~, ~, odg.wcp(segId-1), ~] = evalc('audioqual(segments{segId}.originSig, segments{segId}.wcpRecSig, fs);');
      end
      if ismember("re-CP", methodName)
        snr.recp(segId-1) = get_snr(segments{segId}.originSig(~segments{segId}.mask), segments{segId}.recpRecSig(~segments{segId}.mask));
        [~, ~, ~, odg.recp(segId-1), ~] = evalc('audioqual(segments{segId}.originSig, segments{segId}.recpRecSig, fs);');
      end
      if ismember("f-Janssen", methodName)
        snr.fjanssen(segId-1) = get_snr(segments{segId}.originSig(~segments{segId}.mask), segments{segId}.fjanssenRecSig(~segments{segId}.mask));
        [~, ~, ~, odg.fjanssen(segId-1), ~] = evalc('audioqual(segments{segId}.originSig, segments{segId}.fjanssenRecSig, fs);');
      end
    end
    
    mergedSigName = {"originSig", "gapSig"};
    rowNum = length(snr.hybrid);
    methodNum = length(methodName);
    methodColumn = [];
    snrRow = [];
    odgRow = [];
    timeRow = [];

    if ismember("hybrid", methodName)
      segments{1}.hybridRecSig = segments{1}.originSig;
      segments{end}.hybridRecSig = segments{end}.originSig;
      mergedSigName = [mergedSigName, "hybridRecSig"];
      methodColumn = [methodColumn; repmat("hybrid", rowNum, 1)];
      snrRow = [snrRow; snr.hybrid];
      odgRow = [odgRow; odg.hybrid];
      timeRow = [timeRow; time.hybrid];
    end
    if ismember("A-SPAIN", methodName)
      segments{1}.aspainRecSig = segments{1}.originSig;
      segments{end}.aspainRecSig = segments{end}.originSig;
      mergedSigName = [mergedSigName, "aspainRecSig"];
      methodColumn = [methodColumn; repmat("A-SPAIN", rowNum, 1)];
      snrRow = [snrRow; snr.aspain];
      odgRow = [odgRow; odg.aspain];
      timeRow = [timeRow; time.aspain];
    end
    if ismember("w-CP", methodName)
      segments{1}.wcpRecSig = segments{1}.originSig;
      segments{end}.wcpRecSig = segments{end}.originSig;
      mergedSigName = [mergedSigName, "wcpRecSig"];
      methodColumn = [methodColumn; repmat("w-CP", rowNum, 1)];
      snrRow = [snrRow; snr.wcp];
      odgRow = [odgRow; odg.wcp];
      timeRow = [timeRow; time.wcp];
    end
    if ismember("re-CP", methodName)
      segments{1}.recpRecSig = segments{1}.originSig;
      segments{end}.recpRecSig = segments{end}.originSig;
      mergedSigName = [mergedSigName, "recpRecSig"];
      methodColumn = [methodColumn; repmat("re-CP", rowNum, 1)];
      snrRow = [snrRow; snr.recp];
      odgRow = [odgRow; odg.recp];
      timeRow = [timeRow; time.recp];
    end
    if ismember("f-Janssen", methodName)
      segments{1}.fjanssenRecSig = segments{1}.originSig;
      segments{end}.fjanssenRecSig = segments{end}.originSig;
      mergedSigName = [mergedSigName, "fjanssenRecSig"];
      methodColumn = [methodColumn; repmat("f-Janssen", rowNum, 1)];
      snrRow = [snrRow; snr.fjanssen];
      odgRow = [odgRow; odg.fjanssen];
      timeRow = [timeRow; time.fjanssen];
    end
    
    mergedSig = merge_segments(segments, mergedSigName);

    % export metrics as a table
    metricTable = table(repmat(gapLenInMs(gapId), rowNum * methodNum, 1), ...
      repmat(sigName{sigId}, rowNum * methodNum, 1), ...
      repmat((2:segmentNum-1)', methodNum, 1), ...
      methodColumn, snrRow, odgRow, timeRow, ...
      'VariableNames', {'GapLength', 'SigName', 'SegmentNo', 'Method', 'SNR', 'ODG', 'Time'});
    writetable(metricTable, sprintf("./results/csv/%dms-%s.csv", gapLenInMs(gapId), sigName{sigId}));

    audiowrite(sprintf("./results/wav/%dms-%s-origin.wav", gapLenInMs(gapId), sigName{sigId}), mergedSig.originSig, fs);
    audiowrite(sprintf("./results/wav/%dms-%s-gap.wav", gapLenInMs(gapId), sigName{sigId}), mergedSig.gapSig, fs);
    if ismember("hybrid", methodName)
      audiowrite(sprintf("./results/wav/%dms-%s-hybrid.wav", gapLenInMs(gapId), sigName{sigId}), mergedSig.hybridRecSig, fs);
    end
    if ismember("A-SPAIN", methodName)
      audiowrite(sprintf("./results/wav/%dms-%s-aspain.wav", gapLenInMs(gapId), sigName{sigId}), mergedSig.aspainRecSig, fs);
    end
    if ismember("w-CP", methodName)
      audiowrite(sprintf("./results/wav/%dms-%s-wcp.wav", gapLenInMs(gapId), sigName{sigId}), mergedSig.wcpRecSig, fs);
    end
    if ismember("re-CP", methodName)
      audiowrite(sprintf("./results/wav/%dms-%s-recp.wav", gapLenInMs(gapId), sigName{sigId}), mergedSig.recpRecSig, fs);
    end
    if ismember("f-Janssen", methodName)
      audiowrite(sprintf("./results/wav/%dms-%s-fjanssen.wav", gapLenInMs(gapId), sigName{sigId}), mergedSig.fjanssenRecSig, fs);
    end
  end
end
