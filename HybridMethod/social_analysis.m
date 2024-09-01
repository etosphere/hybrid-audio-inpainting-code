function [x, residual] = social_analysis(F, insig, lambda, mask, shrinkParam, varargin)
%SOCIAL_ANALYSIS Analysis variant of social sparse decomposition.

% define initial parameters
maxIter = 100;
tolerance = 0.33; % relative error tolerance
isTransient = false;
if any(strcmp(varargin, 'tolerence'))
  toleranceIdx = find(strcmp(varargin, 'tolerence'));
  tolerance = varargin{toleranceIdx + 1};
end
if any(strcmp(varargin, 'transient'))
  isTransient = true;
end
if ~isfield(shrinkParam, 'neighbor')
  shrinkParam.shrinkage = 'PEW';
  shrinkParam.neighbor = ones(3, 5);
  shrinkParam.center = [2, 3];
end

F = frameaccel(F, length(insig)); % accelerate the frame
F.frana = @(insig)comp_dgtreal(insig,F.g,F.a,F.M,F.kv.lt,F.flags.do_timeinv);
F.frsyn = @(insig)comp_idgtreal(insig,F.g,F.a,F.M,F.kv.lt,F.flags.do_timeinv);
sigLen = F.L; % input signal length after postpadding
sig = postpad(insig, sigLen);
[~, frameFactor] = framebounds(F, sigLen); % normalization factor of frame operator (||DD^H||)
threshold = lambda / frameFactor;

switch shrinkParam.shrinkage
  case 'L'
    shrink = struct('type', 'l', 'lambda', threshold, 'mu', 1, 'glabel', 'time', 'neigh', 1, 'center', [1, 1], 'expo', 1, 'orth', 0);
  case 'WGL'
    shrink = struct('type', 'l', 'lambda', threshold, 'mu', 1, 'glabel', 'time', 'neigh', ones(3, 7), 'center', [2, 4], 'expo', 1, 'orth', 0);
  case 'EW'
    shrink = struct('type', 'l', 'lambda', threshold, 'mu', 1, 'glabel', 'time', 'neigh', 1, 'center', [1, 1], 'expo', 2, 'orth', 0);
  case 'PEW'
    shrink = struct('type', 'l', 'lambda', threshold, 'mu', 1, 'glabel', 'time', 'neigh', shrinkParam.neighbor, 'center', shrinkParam.center, 'expo', 2, 'orth', 0);
  otherwise
    error('Invalid shrinkage operator!')
end

sigMask = true(size(sig));
sigMask(~mask) = false;

iter = 1;
relativeErr = Inf;

x = sig;
u = F.frana(sig);
uSyn = sig;

tau = 1.5;
tauPrev = tau;
sigma = 1 / tau;
rho = 1;
shrink.lambda = shrink.lambda ./ sigma;

while (iter < maxIter && (relativeErr > tolerance || iter <= 4))
  xPrev = x;
  tempSig = x;
  
  tempSig(~sigMask) = 0;
  g1 = tempSig - sig;
  g = g1;
  v = u + sigma * F.frana(x-tau*(g + uSyn)); % v(i) = u(i) + σA(x(i) − τg(i) − τA∗u(i))
  uHalf = sigma * gen_thresh_gap(v./sigma, shrink, isTransient); % Moreau identity
  uHalf = v - uHalf;
  uSynHalf = F.frsyn(uHalf);

  % x,u,l updates // rho = (k-1)/(k+5);
  x = x - rho * tau * (g + uSynHalf);
  u = u + rho * (uHalf - u);
  uSyn = uSyn + rho * (uSynHalf - uSyn);

  % relativeErr = norm(xPrev - x);
  relativeErr = norm(xPrev(sigMask) - x(sigMask));
  iter = iter + 1;
end

residual = sig - x;
residual = residual(1:length(insig));
end
