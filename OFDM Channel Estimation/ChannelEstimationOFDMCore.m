function [results] = ChannelEstimationOFDMCore(sParams)

%% CD
% cd C:\GitHub\Sandboxes\OFDM

%% Defaults
if ~exist('sParams','var')
    sParams = Define_default_params_channel_estimation_ofdm();
end

%% Channel Estimation
results = cChannelEstimationWrapper(sParams).Apply();

end
