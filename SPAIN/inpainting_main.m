%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%      INPAINTING MAIN FILE      %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% comparison of SPAIN variants:
%     A-SPAIN
%     S-SPAIN OMP
%     S-SPAIN H

close all;
clear variables;
clc;

%% input file settings
fprintf('Loading audio: ')

% available samples:
%     acoustic_guitar
%     double_bass
%     flute
%     organ
%     piano
%     PianoVoice_clean
%     sample
%     string_quartet
%     violin
%     vivaldi_storm

sample = 'simple_test5';
% [data, fs] = audioread(strcat('samples/', sample, '.wav'));
saveFigure = false;
saveAudioFile = false;

[data, fs] = audioread(strcat('test/', sample, '.wav'));
% fs = 44100;
% sample = 'test-signal';
signalF0 = 440;
signalDeltaFreq = 880;
getTwoComponentSignal = @(f0, deltaF, x) ...
  ((sin(2 * pi * f0 .* x) + sin(2 * pi * (f0 + deltaF) .* x)) ./ 2)';
% data = getTwoComponentSignal(signalF0, signalDeltaFreq, (0:1/fs:5-1/fs));
param.Ls = length(data);

fprintf(strcat(sample, '.wav\n'));

%% gap parameters
gap_length = 1600; % gap length in samples
gap_start = 2; % gap start in s

%% settings
fprintf('Setting up the frame parameters\n')

% frame parameters
red = 4; % redundancy of DFT in each block
param.w = 2800; % window length (1200 in original version)
param.a = param.w / 4; % window shift (300 in original version)
param.wtype = 'hann';
param.M = param.w;
param.offset = 'full'; % 'half' / 'full' / 'none'

% some customize parameters
param.signalFreqs = [signalF0, signalF0 + signalDeltaFreq];
param.isDrawing = true;
param.crossFadeLen = 0;
% for time-domain compensation
TDCparamsolver.segs = 10;
TDCparamsolver.gaps = 4;
TDCparamsolver.shift = param.w / 2;
TDCparamsolver.lens = gap_length / 4;
enableTDC = true;

% The offset parameter ensures that the processed signal blocks are
% distributed symmetrically with respect to the center of the gap.
% It is done by moving the starting point of the transform in time domain.
% Offset 'full' corresponds to the situation where middle of the gap
% overlaps with middle index of a block. Offset 'half', on the other hand,
% has the middle of the gap on the axis of two adjacent blocks.

% DFT parameters
param.F = frame('dft');
param.F.redundancy = red; % non-native, our own additional parameter
param.F.frana = @(insig)dft([insig; zeros(length(insig)*(param.F.redundancy - 1), 1)]);
param.F.frsyn = @(insig)postpad(idft(insig), length(insig)/param.F.redundancy);

% paramsolver parameters
paramsolver.s = 1; % increment of k
paramsolver.r = 1; % every r-th iteration increment k by s
paramsolver.epsilon = 0.1; % stopping criterion of termination function
paramsolver.maxit = ceil(floor(param.w*param.F.redundancy/2+1)*paramsolver.r/paramsolver.s); % maximum number of iterations
paramsolver.store_snr = true;
paramsolver.store_obj = false;

%% degrading the signal
fprintf('Generating degraded signal\n\n')
fprintf('   Gap start:   %3.2f s\n', gap_start);
fprintf('   Gap length: %4.2f ms (%d samples)\n', 1000*gap_length/fs, gap_length);
s = round(fs*gap_start);
h = gap_length;
f = s + h - 1;
param.mask = true(param.Ls, 1);
param.mask(s:f) = false;
data_gapped = gapping(data, param.mask); % degrading the original signal

%% re-degradation
degradationType = 'gapping';
switch degradationType
  case 'gapping'
    % do nothing
  case 'clipping'
    data_gapped = data;
    clippingParam.preNormalization = 6;
    data_gapped(s:f) = degradationUnit_applyClipping(data(s:f), fs, 0, clippingParam);
  case 'noise'
    data_gapped = data;
    noiseParam.normalizeOutputAudio = 1;
    noiseParam.snrRatio = 12;
    noiseParam.noiseColor = 'white';
    data_gapped(s:f) = degradationUnit_addNoise(data(s:f), fs, 0, noiseParam);
  case 'compression'
    data_gapped = data;
    data_gapped(s:f) = data(s:f) .* 0.1;
  case 'mp3'
    data_gapped = data;
    mp3Param.LameOptions = '--preset cbr 64';
    data_gapped(s:f) = degradationUnit_applyMp3Compression(data(s:f), fs, 0, mp3Param);
  case 'saturation'
    data_gapped = data;
    saturationTime = 3;
    for n = 1:saturationTime
      data_gapped(s:f) = sin(pi * data_gapped(s:f));
    end
  case 'quantization'
    data_gapped = data;
    bitRate = 4;
    quantization_data = uencode(data(s:f), bitRate, 1, 'signed'); 
    data_gapped(s:f) = double(quantization_data) ./ (2^(bitRate-1)); 
  otherwise
    % do nothing
end

%% pre-select atoms according to the adjacent frames around the gap
leftRange = (s-param.w:s-1);
rightRange = (f+1:f+param.w);

param.atomCoefMask = preselectAtoms(data_gapped(leftRange), data_gapped(rightRange), param);

%% shortening the degraded signal
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

% expand range for TDC
if enableTDC == true
  q = q - param.w * (ceil(TDCparamsolver.gaps / 2) + 1);
  Q = Q + param.w * (ceil(TDCparamsolver.gaps / 2) + 1);
  param.Ls = Q - q + 1;
end

mask = param.mask; % for SNR evaluation
param.mask = param.mask(q:Q); % for SPAIN algorithm, where we work with shortened signal
% initialization of restored signal...
data_rec_1 = data_gapped; % ... for A-SPAIN
data_rec_2 = data_gapped; % ... for S-SPAIN OMP
data_rec_3 = data_gapped; % ... for S-SPAIN H
data_gapped = data_gapped(q:Q); % shortening gapped data for the algorithm in the same way as with the mask
fprintf('\n-----------------------------------------------\n');

%% signal descriptions
yyaxis left;
plot(data_gapped, Color='#9F9F9F');
hold on;
plot(amplitude_envelope(data_gapped), '-', LineWidth=1);
hold on;
yyaxis right;
plot(zerocrossrate_envelope(data_gapped), LineWidth=1);
xlim([0, length(data_gapped)]);
xlabel('Samples (n)');
legend('Degraded Data', 'Amplitude Envelope', 'Zero-cross Rate');
xline([s-q, f-q], '--', LineWidth=1, HandleVisibility="off");

%% A-SPAIN with pre-selection
param.algorithm = 'aspain';
fprintf('\nStarting the A-SPAIN algorithm with pre-selection\n\n');

data_rec_1_pre = data_rec_1;
t = tic;
[data_rec_1_pre(q:Q), snr_iter_1_pre, obj_iter_1_pre] = spain_segmentation(data_gapped, param, paramsolver, data(q:Q));
time = toc(t);

% time
fprintf('   Result obtained in %4.3f seconds.\n', time);
% SNR in the gap
snr_rec = snr_n(data(~mask), data_rec_1_pre(~mask));
fprintf('   SNR of the reconstructed signal is %4.3f dB.\n', snr_rec);
fprintf('\n------------------------------------------------\n');

%% Gradual A-SPAIN
fprintf('\nStarting the Gradual A-SPAIN algorithm\n\n');
data_rec_1_gradual = data_rec_1;

gapStartIdx = s;
gapEndIdx = f;

gappedData = data_gapped;
gradualParam = param;

t = tic;
while gapStartIdx < gapEndIdx
  [data_rec_1_gradual(q:Q), ~, ~] = spain_segmentation(gappedData, gradualParam, paramsolver, data(q:Q));
  
  gapStartIdx = gapStartIdx + gradualParam.a;
  gapEndIdx = gapEndIdx - gradualParam.a;
  gradualParam.mask = true(gradualParam.Ls, 1);
  gradualParam.mask(gapStartIdx-q:gapEndIdx-q) = false;
  gappedData = gradualParam.mask .* data_rec_1_gradual(q:Q);
end
time = toc(t);
% time
fprintf('   Result obtained in %4.3f seconds.\n', time);

snr_rec = snr_n(data(~mask), data_rec_1_gradual(~mask));
fprintf('   SNR of the reconstructed signal is %4.3f dB.\n', snr_rec);
fprintf('\n------------------------------------------------\n');

%% A-SPAIN (plain)
param.atomCoefMask = true(size(param.atomCoefMask));

param.algorithm = 'aspain';
fprintf('\nStarting the A-SPAIN algorithm\n\n');

t = tic;
[data_rec_1(q:Q), snr_iter_1, obj_iter_1] = spain_segmentation(data_gapped, param, paramsolver, data(q:Q));
time = toc(t);

% time
fprintf('   Result obtained in %4.3f seconds.\n', time);
% SNR in the gap
snr_rec = snr_n(data(~mask), data_rec_1(~mask));
fprintf('   SNR of the reconstructed signal is %4.3f dB.\n', snr_rec);
fprintf('\n------------------------------------------------\n');

%% plain A-SPAIN adding TDC
fprintf('\nStarting the A-SPAIN algorithm with TDC\n\n');

data_rec_1_TDC = data_rec_1;

t = tic;
data_rec_1_TDC(q:Q) = tdc(data_rec_1(q:Q), param.mask, param, paramsolver, TDCparamsolver);
time = toc(t);

% time
fprintf('   Result obtained in %4.3f seconds.\n', time);

snr_rec = snr_n(data(~mask), data_rec_1_TDC(~mask));
fprintf('   SNR of the reconstructed signal is %4.3f dB.\n', snr_rec);
fprintf('\n------------------------------------------------\n');

%% A-SPAIN with filterbank
param.algorithm = 'aspain';
fprintf('\nStarting the A-SPAIN algorithm with filterbank\n\n');

data_rec_1_filterbank = data_gapped;

[filters, rate, centerFreqs] = cqtfilters(fs, 400, fs / 2, 1, Q-q+1, 'uniform');
tightFilters = filterbankrealtight(filters, rate, filterbanklength(Q-q+1, rate));
filterCoef = filterbank(data_rec_1_filterbank, tightFilters, rate);
recFilterCoef = filterCoef;
originCoef = filterbank(data(q:Q), tightFilters, rate);

tempDataRec = zeros(size(data_gapped));

gapLeftIdx = ceil((s-q) / rate(1));
gapRightIdx = ceil((f-q) / rate(1));
filterbankParam = param;
filterbankParam.Ls = length(filterCoef{1});
filterbankParam.mask = true(size(filterCoef{1}));
filterbankParam.mask(gapLeftIdx:gapRightIdx) = false;

t = tic;
for band = 1:length(filterCoef)
  fprintf('Processing band %d...\n', centerFreqs(band));
  [tempDataRec, ~, ~] = spain_segmentation(real(filterCoef{band}), filterbankParam, paramsolver, real(originCoef{band}));
  tempDataRec = real(tempDataRec) + imag(filterCoef{band});
  recFilterCoef{band} = tempDataRec;
end
time = toc(t);

data_rec_1_filterbank = ifilterbank(recFilterCoef, tightFilters, rate, 'real');

% time
fprintf('   Result obtained in %4.3f seconds.\n', time);
% SNR in the gap
% snr_rec = snr_n(data(~mask), data_rec_1_filterbank(~mask(q:Q)));
% fprintf('   SNR of the reconstructed signal is %4.3f dB.\n', snr_rec);
fprintf('\n------------------------------------------------\n');

%% S-SPAIN OMP
param.algorithm = 'sspain';
paramsolver.f_update = 'OMP';
fprintf('\nStarting the S-SPAIN OMP algorithm\n(this may take up to a few minutes)\n\n');

t = tic;
% [data_rec_2(q:Q), snr_iter_2, obj_iter_2] = spain_segmentation(data_gapped, param, paramsolver, data(q:Q));
time = toc(t);

% time
fprintf('   Result obtained in %4.3f seconds.\n', time);
% SNR in the gap
snr_rec = snr_n(data(~mask), data_rec_2(~mask));
fprintf('   SNR of the reconstructed signal is %4.3f dB.\n', snr_rec);
fprintf('\n-----------------------------------------------\n');

%% S-SPAIN H
param.algorithm = 'sspain';
paramsolver.f_update = 'H';
fprintf('\nStarting the S-SPAIN H algorithm\n\n');

t = tic;
[data_rec_3(q:Q), snr_iter_3, obj_iter_3] = spain_segmentation(data_gapped, param, paramsolver, data(q:Q));
time = toc(t);

% time
fprintf('   Result obtained in %4.3f seconds.\n', time);
% SNR in the gap
snr_rec = snr_n(data(~mask), data_rec_3(~mask));
fprintf('   SNR of the reconstructed signal is %4.3f dB.\n', snr_rec);
fprintf('\n-----------------------------------------------\n');

%% Janssen
fprintf('Starting Janssen...\n');

% solver settings
janssenParamSolver.Nit = 50; % number of iterations
% Janssenparamsolver.p = 2*a;    

% algorithm
problemData.x = data_gapped;
problemData.IMiss = ~param.mask;

% shortening the data
problemData.x = problemData.x(s-q-param.w*2:f-q+param.w*2);
problemData.IMiss = problemData.IMiss(s-q-param.w*2:f-q+param.w*2);

t = tic;
recoveredJanssen = inpaintFrame_janssenInterpolation(problemData, janssenParamSolver);
% recoveredJanssen = janssen(data_gapped, param, janssenParamSolver);
time = toc(t);

% updating the segment solution
data_rec_2(s-param.w*2:f+param.w*2) = recoveredJanssen;

% time
fprintf('   Result obtained in %4.3f seconds.\n', time);
% SNR in the gap
snr_rec = snr_n(data(~mask), data_rec_2(~mask));
fprintf('   SNR of the reconstructed signal is %4.3f dB.\n', snr_rec);
fprintf('\n-----------------------------------------------\n');

%% plot of original signal and reconstruction
figure();
sliceRange = (s-200:f+200);
t = sliceRange / fs;
plot(t, data(sliceRange), Color='#CFCFCF', LineWidth=1);
hold on;
plot(t, data_rec_1(sliceRange), Color='#0072BD', LineWidth=1);
hold on;
plot(t, data_rec_1_pre(sliceRange), Color='#D95319', LineWidth=1);
% plot(t, data_rec_2(sliceRange), '--');
% hold on;
% plot(t, data_rec_3(sliceRange));
% hold on;
xline([s, f] ./ fs, Color="#77AC30", LineWidth=1);
title({'Original signal and reconstruction using SPAIN'; '(plotting the filled gap and surroundings)'});
% legend('original', 'A-SPAIN', 'S-SPAIN H');
% legend('original', 'A-SPAIN', 'S-SPAIN OMP', 'S-SPAIN H');
% legend('original', 'A-SPAIN', 'Janssen', 'S-SPAIN H');
legend('original', 'A-SPAIN', 'A-SPAIN with pre-selection');
xlim([t(1), t(end)]);
xlabel('time [s]');

if saveFigure == true
  savefig(strcat('results/time_', num2str(gap_length), '_', ...
    num2str(signalF0), '+', num2str(signalDeltaFreq), '_', ...
    num2str(snr_n(data(~mask), data_rec_1(~mask)), '%.2f'), '_', ...
    num2str(snr_n(data(~mask), data_rec_2(~mask)), '%.2f'), '_', ...
    num2str(snr_n(data(~mask), data_rec_3(~mask)), '%.2f'), '.fig'));
end

%% plot objective functions
% randomly choose one of the processed blocks
% [~, n] = size(obj_iter_1);
% indexes = 1:n;
% processed_blocks = indexes(~isnan(obj_iter_1(1, :)));
% block = processed_blocks(randi(length(processed_blocks)));

% plot
% figure;
% plot(obj_iter_1(:, block));
% hold on;
% plot(obj_iter_2(:, block));
% hold on;
% plot(obj_iter_3(:, block));

% legend('A-SPAIN', 'S-SPAIN OMP', 'S-SPAIN H');
% title(sprintf('Norm of the residuum in randomly chosen block (number %d)', block));
% xlabel('Iteration number');
% ylabel('Norm of the residuum');

% if saveFigure == true
%   savefig(strcat('results/residuum_', num2str(gap_length), '_', ...
%     num2str(signalF0), '+', num2str(signalDeltaFreq), '_', ...
%     num2str(snr_n(data(~mask), data_rec_1(~mask)), '%.2f'), '_', ...
%     num2str(snr_n(data(~mask), data_rec_2(~mask)), '%.2f'), '_', ...
%     num2str(snr_n(data(~mask), data_rec_3(~mask)), '%.2f'), '.fig'));
% end

%% plot spectrograms
% frame setting for the spectrogram
number = 100;  % dynamic range of the plots
margin = 1000;
F_plot = frametight(frame('dgtreal', {param.wtype, param.w}, param.a, param.M, 'timeinv'));

framePlotYLim = 4000;

% spectrograms as subplots
figure();
subplot(2, 2, 1);
colormap jet
plotframe(F_plot, frana(F_plot, data), fs, number);
ylim([20, framePlotYLim]);
xlim([s - margin, f + margin]/fs);
xline([s, f] ./ fs, Color="#77AC30", LineWidth=1);
title('Original signal');

subplot(2, 2, 2);
colormap jet
plotframe(F_plot, frana(F_plot, data_rec_1), fs, number);
ylim([20, framePlotYLim]);
xlim([s - margin, f + margin]/fs);
xline([s, f] ./ fs, Color="#77AC30", LineWidth=1);
title('A-SPAIN');

subplot(2, 2, 3);
colormap jet
plotframe(F_plot, frana(F_plot, data_rec_2), fs, number);
ylim([20, framePlotYLim]);
xlim([s - margin, f + margin]/fs);
xline([s, f] ./ fs, Color="#77AC30", LineWidth=1);
% title('S-SPAIN OMP');
title('Janssen');

subplot(2, 2, 4);
colormap jet
plotframe(F_plot, frana(F_plot, data_rec_1_pre), fs, number);
ylim([20, framePlotYLim]);
xlim([s - margin, f + margin]/fs);
xline([s, f] ./ fs, Color="#77AC30", LineWidth=1);
title('A-SPAIN with pre-selection')
% title('S-SPAIN H');

if saveFigure == true
  savefig(strcat('results/spectrogram_', num2str(gap_length), '_', ...
    num2str(signalF0), '+', num2str(signalDeltaFreq), '_', ...
    num2str(snr_n(data(~mask), data_rec_1(~mask)), '%.2f'), '_', ...
    num2str(snr_n(data(~mask), data_rec_2(~mask)), '%.2f'), '_', ...
    num2str(snr_n(data(~mask), data_rec_3(~mask)), '%.2f'), '.fig'));
end

%% save audio files
if saveAudioFile == true
  normalizeAmp = @(sig)(sig ./ max(abs(sig)));
  audiowrite(strcat('results/', sample, '_ASPAIN.wav'), normalizeAmp(data_rec_1(fs*1:fs*3)), fs);
  audiowrite(strcat('results/', sample, '_SSPAINOMP.wav'), normalizeAmp(data_rec_2(fs*1:fs*3)), fs);
  audiowrite(strcat('results/', sample, '_SSPAINH.wav'), normalizeAmp(data_rec_3(fs*1:fs*3)), fs);
end

%% plot the frame coefficients
% data_frame_1 = data_rec_1(q:Q-1);
% frame_coef_1 = frana(param.F, data_frame_1);
% plotframe(param.F, frame_coef_1, 'dynrange', 100);