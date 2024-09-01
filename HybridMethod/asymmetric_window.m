function asymmetricWin = asymmetric_window(winLen, leftLen, rightLen)
%ASYMMETRIC_WINDOW Summary of this function goes here
%   Detailed explanation goes here
symmetricWin = @(t) (1 + cos(pi * (1 + t))) / 2;
% napLog = @(x) (log((10^7) / x) / log((10^7) / (10^7 - 1)));
timeWrap = log2(0.5) / log2(symmetricWin(max(leftLen, rightLen) / (leftLen + rightLen)));

n = linspace(0, 1, winLen+2);
if leftLen > rightLen
  asymmetricWin = symmetricWin(n) .^ timeWrap;
else
  asymmetricWin = 1 - (1 - symmetricWin(n)) .^ timeWrap;
end
asymmetricWin = asymmetricWin(2:end-1);
end
