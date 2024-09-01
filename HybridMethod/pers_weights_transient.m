function [W] = pers_weights_transient(G)

X2 = abs(G).^2;
X2_freq = medfilt2(X2, [25, 1], 'symmetric');
W = X2_freq;

W = sqrt(W);
end


