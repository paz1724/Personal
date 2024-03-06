function [Hout,tap_mat_corr,corrInfo] = GenerateDynamicChannel(sParams)

channel       = sParams.channel;
NFFT          = sParams.NFFT;
Nsymb         = sParams.Nsymb;
fm            = sParams.fm;
scSpacing     = sParams.scSpacing;
SymbolSpacing = sParams.SymbolSpacing;
numRx         = sParams.numRx;
numTx         = sParams.numTx;
corr_name     = sParams.corr_name;
timingOffset  = sParams.timingOffset;
symbolOffset  = sParams.symbolOffset;
isMbsfn       = sParams.isMbsfn;

% channel :  either     channel name string ('AWGN' / 'ETU'  etc)
%                        or                 tap table: [ tap_time_in_ns   tap_power_in_dB ]
% BW in Hz  (used BW, e.g 18e6 for 20MHz)
% NFFT                           -  sets channel BW length  (e.g. 1024 for 10MHz)
% Nsymb                        -  number of symbols to generate
% fm                                -  doppler frequency [Hz]
% SymbolSpacing  -  Time distance between the generated symbols (counted in symbols)
% numRx, numTx        -   if not specified (or empty), 1x1 channel assumed
% corr_name = 'low' / 'med' / 'high'  or [Rx_corr_factor Tx_corr_factor]
% seed                          -   optional, for random number generation
% timing offset [us]    -   optional , for adding a constant phase slope to the channel
% symbolOffset - generate channel with offset in time of 'symbolOffset' symbols

% Defaults
if ~exist('numRx','var') || isempty(numRx), numRx=1;end
if ~exist('numTx','var') || isempty(numTx), numTx=1;end
if ~exist('timingOffset','var') || isempty(timingOffset), timingOffset = 0;end
if ~exist('symbolOffset','var')|| isempty(symbolOffset), symbolOffset=0;end

df = scSpacing; % subcarrier spacing
if strcmp(channel,'AWGN')
    linPhase = exp(-1j*2*pi*df*(timingOffset*1e-6)*((1:NFFT).'-NFFT/2));
    Hout = ones(NFFT,Nsymb,numRx,numTx).*repmat(linPhase,1,Nsymb,numRx,numTx);% add linear phase in F-domain due to timingOffset
%     Hout = ones(NFFT,Nsymb,numRx,numTx);
    tap_mat_corr = nan;
    corrInfo = nan;
    return
end
isTestCQI = strcmpi(channel,'TEST_CQI');


if exist('isMbsfn','var') && isMbsfn
    Ts = 1e-3/12;   % Average symbol time
else
    Ts = 1e-3/14;   % Average symbol time
end


% Channel PDP
if ischar(channel)
   tap_table = Get_tap_table(channel); %[ ns  dB]
else
    tap_table = channel;
end
num_taps = length(tap_table(:,2));
numChannels = numRx*numTx;
tap_time = tap_table(:,1)'*1e-9 + timingOffset*1e-6;
tap_sigma = 10.^(tap_table(:,2)/20);
tap_sigma = tap_sigma / sqrt(sum(abs(tap_sigma).^2)); % noramlize total power


% Correlation Matrix
if exist('corr_name','var') && ~isempty(corr_name)     
    [corr,corrInfo] = Calc_corr_matrix(numRx,numTx,corr_name);
    else
    corr = eye(numRx*numTx);
    corrInfo = [];
end
chol_corr = chol(corr);


% Time points in which the channel will be computed
n_reps = ceil(Nsymb/length(SymbolSpacing));                 
spacingVec  = repmat(SymbolSpacing,1,n_reps);
t = (cumsum(spacingVec) + symbolOffset)*Ts;
t = t(1:Nsymb);


if fm==0 || fm==inf   % extreme cases
    if fm==0            
        randomTaps = 1/sqrt(2) * repmat(rand(num_taps,1,numChannels)  +   1j*randn(num_taps,1,numChannels),1,Nsymb);              % Static channel        
    end
    if fm==inf
        randomTaps =  1/sqrt(2) * (randn(num_taps,Nsymb,numChannels) + 1j*randn(num_taps,Nsymb,numChannels));                      % Taps have no correlation between symbols in time
    end
    
    tap_mat = bsxfun(@times,tap_sigma, randomTaps);  % multiply each tap by its power

    tap_mat_corr = mtimes3d(chol_corr' ,permute(tap_mat,[ 3 1 2])); % create correlation between taps of different channels
    tap_mat_corr = permute(tap_mat_corr,[2 3 1]);                                                     % permute back to desired dimensions
 
else  % 0 < doppler < inf
    
    % Calculate scatterers to create doppler spectrum
    fs              = 1/(Ts*min(SymbolSpacing));
    Nsct            =  max(round(fm*Nsymb/fs),50);
    Nsct            = min(Nsct,500)* (1-isTestCQI) + isTestCQI*1;
    freqs           = fm*cos(2*pi*(1:Nsct)/Nsct).' ;
    FTdop           = exp(-1j*2*pi*freqs*t);                                   % Fourier transorm matrix
    
    phs             = rand(numChannels,num_taps,length(freqs))*2*pi;
    W               = bsxfun(@times,tap_sigma.' , exp(1j*phs)/sqrt(Nsct) );   % random scatteres per tap
    
    % Correlate the scatterers and re-order : numTaps X  Nsct X numChannels
    Wcorr           = mtimes3d(chol_corr' , W);
    Wcorr           = permute(Wcorr, [2 3 1]);
    
    % Fourier transform to time domain
    tap_mat_corr    = mtimes3d(Wcorr , FTdop);    
end

% Fourier transform the taps to frequency domain:
f = (-NFFT/2:NFFT/2-1)'*df;                             % OFDM subcarrier frequencies
FT = exp(-1j*2*pi*f*tap_time);                   % Fourier transform matrix

Hout = zeros(NFFT,Nsymb,numRx,numTx);
for tx=1:numTx
    for rx=1:numRx
        tap_mat = tap_mat_corr(:,:,(tx-1)*numRx + rx);
        Hout(:,:,rx,tx) = FT * tap_mat;
    end
end

