function gapped = gapping (signal, mask)  

% Input parameters
%       signal      clean input signal
%       mask        logical vector indicating the reliable samples

gapped = signal;
gapped(~mask) = 0;

end

