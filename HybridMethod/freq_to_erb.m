function erb = freq_to_erb(freq)
%freq_to_erb Convert frequency to ERB scale
erb = freqtoaud(freq, 'erb');
end
