classdef cChannelEstimationWrapper < handle
    properties
        sParams
    end

    %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
    %% ~~~~~~~~~~~~~~~~~~~~~~ Methods ~~~~~~~~~~~~~~~~~~~~~~ %%
    %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
    methods

        %% ~~~~~~~~~~~~~~~~~~~~~~ Constructor ~~~~~~~~~~~~~~~~~~~~~~ %%
        function obj = cChannelEstimationWrapper(sParams)
            
            obj.sParams = sParams;
        end

        %% ~~~~~~~~~~~~~~~~~~~~~~ Apply ~~~~~~~~~~~~~~~~~~~~~~ %%
        function [results] = Apply(obj)

            %% Update Derivative Parametrer
            obj.Update_derivative_params();

            %% Generate Channel
            rng('default');
            rng(obj.sParams.seed);

            %% cGenerateDynamicChannel
            oGDC = cGenerateDynamicChannel(obj.sParams);

            %% Apply Generate Dynamic Channel
            sChannel = oGDC.Apply();
            H_freq   = sChannel.Hrb;

            %% Generate_complex_noise
            noise = cStat.Generate_complex_noise(obj.sParams);

            %% Add noise to channel
            X = H_freq + noise;

            %% Show Doppler Phase
            obj.Show_channel_phase_and_power(X)

            %% Get_channel_time
            sChannelTime = oGDC.Get_channel_time(X);

            %% Plot_channel_time
            oGDC.Plot_channel_time(sChannelTime, obj.sParams);

            %% Show_image_channel_fast_slow_time
            if obj.sParams.dopplerFreq>0
                oGDC.Show_image_channel_fast_slow_time(sChannelTime);
            end

            %% Channel Estimation
            [sChannelEst, sResults] = obj.Channel_Estimation(oGDC, X);

            %% Results
            results.tauRMS_est           = mean(sChannelEst.tauRMS);
            results.tauRMS_actual        = mean(sResults.tauRMS_actual);
            results.dopplerSpread_est    = sChannelEst.dopplerSpread;
            results.dopplerSpread_actual = obj.sParams.dopplerFreq;
            results.Delay_RMSE           = sResults.sDelayRMSE.Linear;
            results.Delay_RMSE_dB        = sResults.sDelayRMSE.dB;
            results.Doppler_ErrNorm      = sResults.sDopplerErrNorm.Linear;
            results.Doppler_ErrNorm_dB   = sResults.sDopplerErrNorm.dB;

            %% Channel Estimation after removing strongest frequency
            % XInt = oGDC.Calc_H_Freq(sChannelEst.tapDelays, sChannelEst.tapPhasors, sChannelEst.dopplerSpread);
            % XHat = X-XInt;
            % 
            % [sChannelEst, sResults] = obj.Channel_Estimation(oGDC, XHat);

        end

        %% ~~~~~~~~~~~~~~~~~~~~ Update_derivative_params ~~~~~~~~~~~~~~~~~~~~ %%
        function [] = Update_derivative_params(obj)

            %% General Parameters
            obj.sParams.Nsymb            = obj.sParams.totalSfnumber*obj.sParams.NsymbInSf;
            dtDesired_sec                = obj.sParams.dtDesired_us*1e-6;
            obj.sParams.fsc              = obj.sParams.fsc;
            obj.sParams.Nfft             = 2^nextpow2(1/(dtDesired_sec*obj.sParams.fsc));
            obj.sParams.NfftSymbols      = 2^nextpow2(obj.sParams.Nsymb);
            obj.sParams.nfftRx           = min(obj.sParams.Nfft,2^(3+nextpow2(obj.sParams.Nsc)));
            obj.sParams.blockSizeSc      = max(1,floor(obj.sParams.Nsc/obj.sParams.blockSizeScFactor));
            obj.sParams.blockSizeSymb    = max(1,floor(obj.sParams.Nsymb/obj.sParams.blockSizeSymbFactor));
            obj.sParams.Tsymb            = 1/obj.sParams.fsc;
            obj.sParams.noiseSizeVec     = [obj.sParams.Nsc, obj.sParams.NsymbInSf*obj.sParams.totalSfnumber];

            %% Channel Seed
            obj.sParams.seed = obj.Calculate_seed();

        end

        %% ~~~~~~~~~~~~~~~~~~~~ Show_channel_phase_and_power ~~~~~~~~~~~~~~~~~~~~ %%
        function [] = Show_channel_phase_and_power(obj, X)

            %% Show Doppler Phase
            if obj.sParams.bShowPhases
                subcarriersToPlot = 1:4;
                symbolsToPlot = 1:7;
                cGenerateDynamicChannel.Show_channel_phase(X, subcarriersToPlot, symbolsToPlot);
            end

            %% Show Doppler Power
            if obj.sParams.bShowPower
                subcarriersToPlot = 1:4;
                symbolsToPlot = 1:7;
                cGenerateDynamicChannel.Show_channel_power(X, subcarriersToPlot, symbolsToPlot);
            end

        end

        %% ~~~~~~~~~~~~~~~~~~~~ Channel_Estimation ~~~~~~~~~~~~~~~~~~~~ %%
        function [sChannelEst, sResults] = Channel_Estimation(obj, oGDC, X)
            oCE                                 = cChannelEstimation(obj.sParams, oGDC);
            [sChannelEst, sResults] = oCE.Estimate_Channel(X);
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Calculate_seed ~~~~~~~~~~~~~~~~~~~~ %%
        function [seed] = Calculate_seed(obj)

            %% sSeed
            sSeed = obj.sParams;

            %% Defaults
            sSeed                  = cStruct.Set_default_value(sSeed,'baseSeed',222);
            sSeed                  = cStruct.Set_default_value(sSeed,'seedGap',11);
            sSeed                  = cStruct.Set_default_value(sSeed,'nexp',1);
            sSeed                  = cStruct.Set_default_value(sSeed,'nbatch',1);
            sSeed                  = cStruct.Set_default_value(sSeed,'Nbatches',1);

            sSeed                  = cStruct.Set_default_value(sSeed,'overrideSeedEnable',0);
            sSeed                  = cStruct.Set_default_value(sSeed,'overrideSeed',[]);

            %% Parameters
            baseSeed                = sSeed.baseSeed;
            seedGap                 = sSeed.seedGap;
            nexp                    = sSeed.nexp;
            nbatch                  = sSeed.nbatch;
            Nbatches                = sSeed.Nbatches;
            overrideSeedEnable	 	= sSeed.overrideSeedEnable;
            overrideSeed	        = sSeed.overrideSeed;

            %% Calculate Seed
            if overrideSeedEnable && ~isempty(overrideSeed)
                seed                = overrideSeed;
            else
                seed                = baseSeed+seedGap*(nexp-1)*Nbatches+seedGap*(nbatch-1);
            end

        end

    end

    %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
    %% ~~~~~~~~~~~~~~~~~~ Static Methods ~~~~~~~~~~~~~~~~~~~ %%
    %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
    methods(Static)

    end

end
