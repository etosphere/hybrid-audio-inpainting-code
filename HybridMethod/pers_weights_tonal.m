function [W] = pers_weights_tonal(G)

X2 = abs(G).^2;
X2_time = medfilt2(X2, [3, 9], 'symmetric');
X2_freq = medfilt2(X2, [65, 1], 'symmetric');
W = max(X2_time - 0.15 * X2_freq, 1e-8);

W = sqrt(W);
end
