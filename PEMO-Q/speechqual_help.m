% SPEECHQUAL    M. Hansen's method for estimation of perceived speech quality 
%
%   USAGE:
%     qc = speechqual(RefSig, TestSig, fs, [MaxLag, PauseCut, LevelAlign, weighting, assim, verbose])
%
%   INPUT:
%         RefSig : reference signal (either vector or wav-filename) (*) 
%        TestSig : processed (test) signal (either vector or wav-filename) (*)
%             fs : sample rate
%         MaxLag : parameter for automatic delay compensation. MaxLag defines the range of
%                  lags in ms ([-MaxLag, MaxLag]), over which the search for a possible delay
%                  between reference and test signal is performed. Set = 0 to skip the delay
%                  compensation(**). (optional, default = 0)
%       PauseCut : flag for removing silent intervals of at least 200 ms from signals
%                  (optional, default = 0)
%     LevelAlign : flag for compensating a possible overall level difference 
%                  (optional, default = 1)
%          assim : flag for partial assimilation of internal representations
%                  (Beerends-Berger-approach) (optional, default = 0)
%      weighting : flag for Hansen's band importance weighting (optional, default = 1)
%        verbose : flag for displaying parameter settings and progress messages 
%                  (optional, default = 0)
%
%   OUTPUT:
%             qc : objective speech quality measure
%
%   (Keep order of arguments, e.g. if you want to specify weighting, also specify LevelAlign.
%    You may type "[]" to adopt the default value.)
%
%   (*)  A signal amplitude of 1 (full scale in case of wav-files) is assumed to correspond to 100 dB SPL.
%   (**) If the delay is known, it is recommended to align the signals beforehand and skip
%        the automatic delay compensation.
%
%   DESCRIPTION:
%     SPEECHQUAL is an implementation of Hansen's method for speech quality
%     estimation. It predicts the perceived quality of a given test speech signal
%     relative to that of a reference signal, using the model of auditory perception
%     ("PEMO") by Dau et al. (1996).
%     After prealignment of the given pair of speech signals, the auditory model
%     transforms both signals into corresponding internal representations.
%     A set of linear weights is applied to the frequency bands of the two
%     representations to account for different assumed importances for the
%     perceived speech quality. (Can be disabled.)
%     Following an approach of Beerends/Berger, an optional subsequent 
%     assimilation procedure (not contained in Hansen's original method) can be applied
%     to reduce the differences between the internal representations sign-dependently,
%     putting more weight on positive than on negative deviations of the processed
%     signal (default: no assimilation). 
%     Finally, the linear cross correlation coefficient of the pair of internal
%     representations, "qc", is calculated and serves as a measure of the 
%     perceptual similarity. If the reference signal is of high quality,
%     qc may be interpreted as an objective measure of the speech quality (degradation)
%     of the test signal.
%
%   Ref.: Hansen and Kollmeier, J. Audio Eng. Soc., Vol. 48, No. 5, Mai 2000 
%
%   Copyright University of Oldenburg 2022
%   Author: Rainer Huber (rainer.huber@uni-oldenburg.de)
