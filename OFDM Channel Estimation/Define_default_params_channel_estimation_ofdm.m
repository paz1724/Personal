function [paramStruct] = Define_default_params_channel_estimation_ofdm()
																
%% Scenario
paramStruct.channelType    = 'AWGN';
paramStruct.tapTimeVec_ns  = 0;
paramStruct.tapPowerVec_dB = 0;
paramStruct.dopplerFreq    = 0;
paramStruct.Nsc            = 600;
paramStruct.totalSfnumber  = 5;
paramStruct.SNR_dB         = 20;
paramStruct.baseSeed       = 111;

%% Channel Estimation Methods
paramStruct.ceMethodDelay   = "MUSIC_Inverse"; % MUSIC_Inverse , MUSIC_Roots , ESPRIT , OMP , MVDR
paramStruct.ceMethodDoppler = "Taps_MUSIC";    % Taps_MUSIC , ACF_MUSIC , ACF_FFT

%% Channel Estimation Parameters
paramStruct.blockSizeScFactor   = 2; % For Sc covariance matrix spatial smoothing
paramStruct.blockSizeSymbFactor = 2; % For Symb covariance matrix spatial smoothing
paramStruct.lambdaFloorLevelMDL = 1e-5; 
paramStruct.modeMDL             = "MDL"; % MDL , Eig                 
paramStruct.MDLthFactor         = 0;

%% SubSpace
paramStruct.numPathsDelay         = 0; % 0 means automatic find the optimal number of paths using MDL
paramStruct.numPathsDoppler       = 0;
paramStruct.dtDesired_us          = 0.05;
paramStruct.dfDesired_hz          = 1;
paramStruct.maxRelevantT_us       = 7;
paramStruct.maxDopplerFreq_Hz     = 1500;
paramStruct.maxIterNewtonMethod   = 10; % Maximum number of iterations for Newton's method
paramStruct.toleranceNewtonMethod = 1e-3; % Tolerance for root convergence
paramStruct.numStdsTh             = 6;
paramStruct.eigMode               = "Matlab"; % Matlab , QR
paramStruct.bWhitenNoise          = false;
paramStruct.thOMP                 = 1e-5;
paramStruct.modeOMP               = "MVDR"; % MVDR , FFT

%% Booleans
paramStruct.bShowPhases = false;
paramStruct.bShowPower  = false;

%% General
paramStruct.overrideSeed     = [];
paramStruct.correlationType  = 'low';
paramStruct.cellId           = 1;
paramStruct.timingOffset_us  = 0;
paramStruct.freqOffset_Hz    = 0;
paramStruct.nTx              = 1;
paramStruct.numRx            = 1;
paramStruct.symbOffset       = 0;
paramStruct.subframeTime_sec = 1e-3;
paramStruct.nfftRx           = 128;
paramStruct.symbSpacing      = 1;
paramStruct.fsc              = 15e3; %[Hz]
paramStruct.NsymbInSf        = 14;
paramStruct.nNRSinSfPerTx    = 8;

%% Not Used
paramStruct.NREsNRSinSf = 16;
paramStruct.cp          = 'Normal';


end