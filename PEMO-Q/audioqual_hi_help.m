% AUDIOQUAL_HI    Estimation of perceived audio quality following the PEMO-Q method
%                 Extended version for normal hearing and hearing impaired listeners
%
%   USAGE:
%     [PSM, PSMt, ODG, PSM_inst] = audioqual_hi(RefSig, TestSig, fs, [HL, dt, modproc, MaxLag, PauseCut, LevelAlign, assim, verbose])
%
%   INPUT:
%         RefSig : reference signal (either vector or wav-filename) (*)
%        TestSig : test signal (either vector or wav-filename) (*) 
%             fs : sample rate (optional if signals are specified by wav-filenames)
%             HL : Hearing loss, i.e. audiogram data: 
%				   2xN matrix; first row: frequencies in Hz, second row: pure tone hearing threshold in dB HL
%                  set HL = 0 or HL = [] to indicate normal hearing (optional, default = 0)			
%             dt : time resolution (correlation interval length) in ms
%                  (optional, default = 10)
%        modproc : type of modulation processing in PEMO: 
%                  either filterbank ('fb') or lowpass ('lp') (default: 'lp')
%                  'fbN' with 1 <= N <= 10 applies a modulation filterbank with N filters
%                  (optional, default: N = 8)
%         MaxLag : parameter for automatic delay compensation. MaxLag defines the range of
%                  lags in ms ([-MaxLag, MaxLag]), over which the search for a possible delay
%                  between reference and test signal is performed. Set = 0 to skip the delay
%                  compensation(**). (optional, default = 0)
%       PauseCut : flag for removing silent intervals of at least 200 ms from signals
%                  (optional, default = 0)
%     LevelAlign : flag for compensating a possible overall level difference 
%                  (optional, default = 1)
%          assim : flag for partial assimilation of internal representations
%                  (Beerends-Berger-approach) (optional, default = 1)
%        verbose : flag for displaying parameter settings and progress messages 
%                  (optional, default = 0)
%
%   OUTPUT:
%            PSM : overall objective quality measure (overall correlation between internal representations)
%           PSMt : overall objective quality measure (= 5. percentile of "internal activity"-weighted PSM(t))
%            ODG : "Objective Difference Grade", i.e transformed PSMt-measure
%       PSM_inst : instantaneous objective quality measure (vector)
%                  (i.e. PSM_inst = PSM(t) with t = n*dt ms, n = 1,2,..., dt = 10 per default)
%
%   (Keep order of arguments, e.g. if you want to specify modproc, also specify dt.
%    You may type "[]" to adopt the default value.)
%
%   (*)  A signal amplitude of 1 (full scale in case of wav-files) is assumed to correspond to 100 dB SPL.
%   (**) If the delay is known, it is recommended to align the signals beforehand and skip
%        the automatic delay compensation.
%
%   DESCRIPTION:
%     AUDIOQUAL_HI is the implementation of the method for objective perceptual
%     assessment of audio quality "PEMO-Q" of Huber and Kollmeier (2006), slightly 
%     modified and extended for hearing impaired following the approach of 
%     Derleth et al. (2001). It predicts the perceived audio quality of a given test
%     signal relative to that of a reference signal, using the  model of auditory
%     perception ("PEMO") by Dau et al. (1996, 1997).
%     After prealignment of the given pair of audio signals, the auditory model
%     transforms both signals into corresponding internal representations.
%     Following an approach of Beerends/Berger, an optional subsequent 
%     assimilation procedure reduces the differences between the internal
%     representations sign-dependently, putting more weight on positive than on
%     negative deviations of the processed signal.
%     Finally, the linear cross correlation coefficient of the pair of internal
%     representations, PSM, is calculated and serves as a measure of the 
%     perceived similarity. If the reference signal is of high audio quality,
%     PSM may also be interpreted as an objective measure of the audio quality
%     (degradation) of the test signal. 
%     In addition, a time series of instantaneous audio quality is computed 
%     by consecutively correlating slices of the internal representations (length
%     of the slices = dt ms, 10 ms per default). The output vector PSM_inst contains
%     the sequence of PSM(n*dt).
%     The overall quality measure PSMt is calculated by taking the 5. percentile of
%     PSM_inst, after it is weighted by the moving average of the internal representation
%     of the test signal ("internal activity" weighting).
%     The third output argument, ODG, is obtained by mapping PSMt to a value that
%     corresponds to the subjective quality scale SDG ("Subjecitve Difference Grade")
%     defined in ITU-R BS.1116. The mapping function PSMt -> ODG was derived empirically
%     by fitting the results of several listening tests. THIS MAPPING ONLY HOLDS IF      
%     DEFAULT SETTINGS OF dt, LevelAlign, assim AND WIDE BAND SIGNALS ARE USED! 
%     (Regarding modproc, settings "lp", "fb" and "fb8" are possible.)
%     OTHERWISE, ODG WILL NOT BE COMPUTED! Moreover, be aware that this mapping might 
%     also not hold for applications other than audio codec evaluation. 
%     The values of SDG and ODG have the following meanings:
%     The quality impairment of the test signal relative to the refence signal is
%               0 : Imperceptible 
%              -1 : Perceptible but not annoying  
%              -2 : Slightly annoying 
%              -3 : Annoying
%              -4 : Very annoying
%
%   Copyright University of Oldenburg 2022
%   Author: Rainer Huber (rainer.huber@uni-oldenburg.de)
