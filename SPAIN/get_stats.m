files = dir(fullfile('result/', 'time_*.fig'));

filenames = struct2cell(files);
filenames = filenames(1,:);

data = zeros(length(files), 5);

for idx = (1:length(files))
  filename = files(idx).('name');
  filename = filename(6:end-4);
  splitedFilename = regexp(filename,'_','split');
  gapLength = str2double(cell2mat(splitedFilename(1)));
  deltaFreq = regexp(splitedFilename(2), '+', 'split');
  deltaFreq = str2double(deltaFreq{1}{2});
  aspainSNR = str2double(cell2mat(splitedFilename(3)));
  janssenSNR = str2double(cell2mat(splitedFilename(4)));
  sspainhSNR = str2double(cell2mat(splitedFilename(5)));
  data(idx,:) = [gapLength, deltaFreq, aspainSNR, janssenSNR, sspainhSNR];
end

subplot(3, 1, 1);
gap600Data = data(data(:, 1) == 600, :);
gap600Data = sortrows(gap600Data(:, 2:end), 1);
plot(gap600Data(:,1), gap600Data(:, 2), gap600Data(:,1), gap600Data(:, 3), gap600Data(:,1), gap600Data(:, 4));
legend('A-SPAIN', 'Janssen', 'S-SPAIN H');
subplot(3, 1, 2);
gap1200Data = data(data(:, 1) == 1200, :);
gap1200Data = sortrows(gap1200Data(:, 2:end), 1);
plot(gap1200Data(:,1), gap1200Data(:, 2), gap1200Data(:,1), gap1200Data(:, 3), gap1200Data(:,1), gap1200Data(:, 4));
subplot(3, 1, 3);
gap2000Data = data(data(:, 1) == 2000, :);
gap2000Data = sortrows(gap2000Data(:, 2:end), 1);
plot(gap2000Data(:,1), gap2000Data(:, 2), gap2000Data(:,1), gap2000Data(:, 3), gap2000Data(:,1), gap2000Data(:, 4));

figure(2);
aspainSNR = data(:, [1, 3]);
aspainSNR = sortrows(aspainSNR);
janssenSNR = data(:, [1, 4]);
janssenSNR = sortrows(janssenSNR);
sspainhSNR = data(:, [1, 5]);
sspainhSNR = sortrows(sspainhSNR);
boxchart(categorical(aspainSNR(:, 1), unique(aspainSNR(:, 1))), aspainSNR(:, 2));
hold on;
boxchart(categorical(janssenSNR(:, 1), unique(janssenSNR(:, 1))), janssenSNR(:, 2));
hold on;
boxchart(categorical(sspainhSNR(:, 1), unique(sspainhSNR(:, 1))), sspainhSNR(:, 2));
legend('A-SPAIN', 'Janssen', 'S-SPAIN H');
xlabel('gap size (samples)');
ylabel('SNR (dB)');
