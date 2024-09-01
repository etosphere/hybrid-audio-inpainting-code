function snrValue = get_snr(originalSig, reconstructedSig)
%GET_SNR Summary of this function goes here
%   Detailed explanation goes here
snrValue = 20*log10(norm(originalSig)/norm(originalSig-reconstructedSig));
end
