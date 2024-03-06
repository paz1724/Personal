classdef cChannelEstimation < handle
    properties

        %% Objects
        oGDC

        %% Properties
        sParams

    end

    %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
    %% ~~~~~~~~~~~~~~~~~~~~~~ Methods ~~~~~~~~~~~~~~~~~~~~~~ %%
    %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
    methods

        %% ~~~~~~~~~~~~~~~~~~~~~~ Constructor ~~~~~~~~~~~~~~~~~~~~~~ %%
        function obj = cChannelEstimation(sParams, oGDC)
            obj.oGDC = oGDC;
            obj      = cStruct.Update_param_struct_into_object(obj, sParams);

        end

        %% ~~~~~~~~~~~~~~~~~~~~ EstimateChannel ~~~~~~~~~~~~~~~~~~~~ %%
        function [sChannelEst, sResults] = Estimate_Channel(obj, X)
            % H_freq is the estimated channel in the frequency domain [NsubCarriers x Nsymbols]

            %% tau_grid and frequencies
            [obj.sParams.tau_grid, obj.sParams.freqVec_Hz] = obj.Generate_time_freq_grids_for_delay_spread(obj.oGDC.Nsc);
            [obj.sParams.freq_grid, obj.sParams.timeVec_sec] = obj.Generate_time_freq_grids_for_doppler_spread(obj.oGDC.Nsymb);

            %% Estimate Channel Delay
            tapDelays  = obj.Estimate_Channel_Delay(X);

            %% Estimate Taps Power And Phase
            tapPhasors = obj.Estimate_Taps_Power_And_Phase(X, tapDelays);

            %% Calc_delay_time_RMS
            [tauRMS, tauRMS_actual, sDelayRMSE] = obj.Calc_delay_time_RMS(tapDelays, tapPhasors);

            %% Estimate_Doppler_Spread
            dopplerSpread = obj.Estimate_Doppler_Spread(tapPhasors);
            sDopplerErrNorm = obj.Calc_Normalized_Error(dopplerSpread, obj.sParams.dopplerFreq);

            %% RMSE
            sResults.tauRMS_actual   = tauRMS_actual;
            sResults.sDelayRMSE      = sDelayRMSE;
            sResults.sDopplerErrNorm = sDopplerErrNorm;

            %% Estimated Channel Struct
            sChannelEst.tapDelays     = tapDelays;
            sChannelEst.tapPhasors    = tapPhasors;
            sChannelEst.tauRMS        = tauRMS;
            sChannelEst.dopplerSpread = dopplerSpread;

        end

        %% ~~~~~~~~~~~~~~~~~~~~ Estimate_Channel_Delay ~~~~~~~~~~~~~~~~~~~~ %%
        function [tapDelays] = Estimate_Channel_Delay(obj, X)
            % H_freq is the estimated channel in the frequency domain [NsubCarriers x Nsymbols]

            %% SubSpace Preparation
            baseTitleStr = vertcat("Spectrum for Channel Tap Delay Estimation",...
                "Channel: "+strrep(obj.oGDC.channel,'_',' '),...
                "SNR: "+obj.sParams.SNR_dB+" [dB]",...
                "N_{Sc}: "+obj.oGDC.Nsc);

            sParamsSubSpace.method              = obj.sParams.ceMethodDelay      ;
            sParamsSubSpace.blockSize           = obj.sParams.blockSizeSc        ;
            sParamsSubSpace.numPaths            = obj.sParams.numPathsDelay      ;
            sParamsSubSpace.lambdaFloorLevelMDL = obj.sParams.lambdaFloorLevelMDL;
            sParamsSubSpace.modeMDL             = obj.sParams.modeMDL;
            sParamsSubSpace.inputGrid           = obj.sParams.freqVec_Hz         ;
            sParamsSubSpace.outputGrid          = obj.sParams.tau_grid           ;
            sParamsSubSpace.tapType             = "Delay"                ;
            sParamsSubSpace.tapUnitsToShow      = "ns"                   ;
            sParamsSubSpace.tapFactorToShow     = 1e9                    ;
            sParamsSubSpace.baseTitleStr        = baseTitleStr           ;
            sParamsSubSpace.bShowPlots          = obj.sParams.plotEnable         ;

            sParamsSubSpace.maxIterNewtonMethod   = obj.sParams.maxIterNewtonMethod;
            sParamsSubSpace.toleranceNewtonMethod = obj.sParams.toleranceNewtonMethod;
            sParamsSubSpace.numStdsTh       = obj.sParams.numStdsTh;
            sParamsSubSpace.eigMode               = obj.sParams.eigMode;
            sParamsSubSpace.bWhitenNoise         = obj.sParams.bWhitenNoise;
            sParamsSubSpace.thOMP                 = obj.sParams.thOMP;
            sParamsSubSpace.modeOMP               = obj.sParams.modeOMP;
   
            %% Estimate Channel Delay Taps
            oCEDelay  = cSubSpace(sParamsSubSpace);
            tapDelays = oCEDelay.Apply(X);

        end

        %% ~~~~~~~~~~~~~~~~~~~~ Estimate Taps Power And Phase ~~~~~~~~~~~~~~~~~~~~ %%
        function tapPhasors = Estimate_Taps_Power_And_Phase(obj, X, tapDelays)

            %% Construct S matrix
            S = exp(-1i*2*pi*tapDelays*obj.sParams.freqVec_Hz).';

            %% Least Squares Solution for the phasors
            tapPhasors = (S'*S) \ S' * X;

        end

        %% ~~~~~~~~~~~~~~~~~~~~ Calc_delay_time_RMS ~~~~~~~~~~~~~~~~~~~~ %%
        function [tauRMS_est, tauRMS_actual, sRMSE] = Calc_delay_time_RMS(obj, tapDelays, tapPhasors)

            %% Tau RMS
            tauRMS_est = obj.Calc_tauRMS(tapPhasors, tapDelays);
            tauRMS_actual = obj.Calc_tauRMS(obj.oGDC.actualChannelPhasors, obj.oGDC.actualChannelDelays);

            sRMSE = obj.Calc_RMSE(tauRMS_est, tauRMS_actual);

        end

        %% ~~~~~~~~~~~~~~~~~~~~ Generate_time_freq_grids_for_delay_spread ~~~~~~~~~~~~~~~~~~~~ %%
        function [tau_grid_sec, freqVec_Hz] = Generate_time_freq_grids_for_delay_spread(obj, Nfreqs)
            % Define the search grid for possible tap delays
            tau_grid_sec = 1e-6*(0:obj.sParams.dtDesired_us:obj.sParams.maxRelevantT_us).';
            freqVec_Hz = (0:Nfreqs-1)*obj.oGDC.scSpacing;
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Generate_time_freq_grids_for_doppler_spread ~~~~~~~~~~~~~~~~~~~~ %%
        function [freq_grid_Hz, timeVec_sec] = Generate_time_freq_grids_for_doppler_spread(obj, Nsymb)
            % Define the search grid for possible tap delays
            freq_grid_Hz = (0:obj.sParams.dfDesired_hz:obj.sParams.maxDopplerFreq_Hz).';

            timeVec_sec = (0:Nsymb-1)*obj.oGDC.Tsymb;
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Calc_tauRMS ~~~~~~~~~~~~~~~~~~~~ %%
        function [tauRMS] = Calc_tauRMS(obj, tapPhasors, tapDelays)
            if numel(tapDelays)==1
                tauRMS = tapDelays;
            else
                tau2 = obj.Calc_tau2(tapPhasors, tapDelays);
                tau1 = obj.Calc_tau1(tapPhasors, tapDelays);

                tauRMS = sqrt(max(0,tau2-tau1.^2));
            end
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Estimate_Doppler_Spread ~~~~~~~~~~~~~~~~~~~~ %%
        function [dopplerSpread] = Estimate_Doppler_Spread(obj, tapPhasors)
            % H_freq is the estimated channel in the frequency domain [NsubCarriers x Nsymbols]

            switch obj.sParams.ceMethodDoppler
                case {"ACF_MUSIC","ACF_FFT"}
                    %% Estimate_Doppler_Spread_ACF
                    dopplerSpread = obj.Estimate_Doppler_Spread_ACF(tapPhasors);
                case "Taps_MUSIC"
                    %% Estimate_Doppler_Spread_Taps_MUSIC
                    dopplerSpread = obj.Estimate_Doppler_Spread_Taps_MUSIC(tapPhasors);
            end
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Estimate_Doppler_Spread_Taps_MUSIC ~~~~~~~~~~~~~~~~~~~~ %%
        function [dopplerSpread] = Estimate_Doppler_Spread_Taps_MUSIC(obj, tapPhasors)
            % H_freq is the estimated channel in the frequency domain [NsubCarriers x Nsymbols]

            %% SubSpace Preparation
            baseTitleStr = vertcat("Spectrum of ACF for Doppler Spread Estimation",...
                "Channel: "+strrep(obj.oGDC.channel,'_',' '),...
                "SNR: "+obj.sParams.SNR_dB+" [dB]",...
                "N_{Sc}: "+obj.oGDC.Nsc);

            sParamsSubSpace.blockSize  = obj.sParams.blockSizeSymb      ;
            sParamsSubSpace.inputGrid  = obj.sParams.timeVec_sec(1:size(tapPhasors,2)) ;
            sParamsSubSpace.outputGrid = obj.sParams.freq_grid          ;

            sParamsSubSpace.method                = obj.sParams.ceMethodDelay        ;
            sParamsSubSpace.numPaths              = obj.sParams.numPathsDoppler    ;
            sParamsSubSpace.lambdaFloorLevelMDL   = obj.sParams.lambdaFloorLevelMDL;
            sParamsSubSpace.modeMDL               = obj.sParams.modeMDL;
            sParamsSubSpace.maxIterNewtonMethod   = obj.sParams.maxIterNewtonMethod;
            sParamsSubSpace.toleranceNewtonMethod = obj.sParams.toleranceNewtonMethod;
            sParamsSubSpace.numStdsTh             = obj.sParams.numStdsTh;
            sParamsSubSpace.eigMode               = obj.sParams.eigMode;
            sParamsSubSpace.bWhitenNoise          = obj.sParams.bWhitenNoise;
            sParamsSubSpace.thOMP                 = obj.sParams.thOMP;
            sParamsSubSpace.modeOMP               = obj.sParams.modeOMP;
            
            sParamsSubSpace.tapType             = "Doppler"              ;
            sParamsSubSpace.tapUnitsToShow      = "Hz"                   ;
            sParamsSubSpace.tapFactorToShow     = 1                      ;
            sParamsSubSpace.baseTitleStr        = baseTitleStr           ;
            sParamsSubSpace.bShowPlots          = obj.sParams.plotEnable         ;

            %% Estimate Channel Delay Taps
            oSubSpace  = cSubSpace(sParamsSubSpace);
            tapDoppler = oSubSpace.Apply(tapPhasors.');

            %% Choose largest steering vector index for the doppler spread
            dopplerSpread   = tapDoppler(end);

        end

        %% ~~~~~~~~~~~~~~~~~~~~ Estimate_Doppler_Spread_ACF ~~~~~~~~~~~~~~~~~~~~ %%
        function [dopplerSpread] = Estimate_Doppler_Spread_ACF(obj, tapPhasors)

            %% Preprocessing Step: Compute the spatially smoothed covariance matrix
            Rxx_smoothed = cMath.Spatial_Smoothing(tapPhasors.', obj.sParams.blockSizeSymb);

            %% Autocorrelation Function
            Rhl = cMath.Calc_ACF_from_Covariance_Matrix(Rxx_smoothed);

            switch obj.sParams.ceMethodDoppler
                case "ACF_FFT"
                    %% Calc_ACF_PSD_By_FFT
                    dopplerSpread = obj.Estimate_Doppler_Spread_ACF_FFT(Rhl);

                case "ACF_MUSIC"
                    %% Estimate_Doppler_Spread_By_ACF_MUSIC
                    dopplerSpread = obj.Estimate_Doppler_Spread_ACF_MUSIC(Rhl);
            end

        end

        %% ~~~~~~~~~~~~~~~~~~~~ Estimate_Doppler_Spread_ACF_MUSIC ~~~~~~~~~~~~~~~~~~~~ %%
        function [dopplerSpread] = Estimate_Doppler_Spread_ACF_MUSIC(obj, Rhl)

            %% MUSIC Preparation
            baseTitleStr = vertcat("MUSIC Spectrum of ACF for Doppler Spread Estimation",...
                "Channel: "+strrep(obj.oGDC.channel,'_',' '),...
                "SNR: "+obj.sParams.SNR_dB+" [dB]",...
                "N_{Sc}: "+obj.oGDC.Nsc);

            sParamsSubSpace.blockSize  = obj.sParams.blockSizeSymb      ;
            sParamsSubSpace.inputGrid  = obj.sParams.timeVec_sec(1:numel(Rhl)) ;
            sParamsSubSpace.outputGrid = obj.sParams.freq_grid          ;

            sParamsSubSpace.method                = obj.sParams.ceMethodDelay        ;
            sParamsSubSpace.numPaths              = obj.sParams.numPathsDoppler    ;
            sParamsSubSpace.lambdaFloorLevelMDL   = obj.sParams.lambdaFloorLevelMDL;
            sParamsSubSpace.modeMDL               = obj.sParams.modeMDL;
            sParamsSubSpace.maxIterNewtonMethod   = obj.sParams.maxIterNewtonMethod;
            sParamsSubSpace.toleranceNewtonMethod = obj.sParams.toleranceNewtonMethod;
            sParamsSubSpace.numStdsTh       = obj.sParams.numStdsTh;
            sParamsSubSpace.eigMode               = obj.sParams.eigMode;
            sParamsSubSpace.bWhitenNoise         = obj.sParams.bWhitenNoise;
            sParamsSubSpace.thOMP                 = obj.sParams.thOMP;
            sParamsSubSpace.modeOMP               = obj.sParams.modeOMP;

            sParamsSubSpace.tapType             = "Doppler"              ;
            sParamsSubSpace.tapUnitsToShow      = "Hz"                   ;
            sParamsSubSpace.tapFactorToShow     = 1                      ;
            sParamsSubSpace.baseTitleStr        = baseTitleStr           ;
            sParamsSubSpace.bShowPlots          = obj.sParams.plotEnable         ;

            %% Estimate Channel Delay Taps
            oSubSpace  = cSubSpace(sParamsSubSpace);
            tapDoppler = oSubSpace.Apply(Rhl);

            %% Choose largest steering vector index for the doppler spread
            dopplerSpread   = tapDoppler(end);

        end
        
        %% ~~~~~~~~~~~~~~~~~~~~ Estimate_Doppler_Spread_ACF_FFT ~~~~~~~~~~~~~~~~~~~~ %%
        function [dopplerSpread] = Estimate_Doppler_Spread_ACF_FFT(obj, Rhl)

            %% Nfft & fVec
            Nfft        = obj.sParams.NfftSymbols*4;
            fs          = 1/obj.sParams.Tsymb;
            sPSD.fVecHz = linspace(-fs/2,fs/2,Nfft);

            %% H
            H_PSD         = abs(fftshift(fft(Rhl,Nfft)));
            H_PSD_dB      = 20*log10(H_PSD);
            sPSD.H_PSD_dB = H_PSD_dB-max(H_PSD_dB);
            [~, maxInd_H] = max(sPSD.H_PSD_dB);
            fd_H_Hz       = sPSD.fVecHz(maxInd_H);

            %% J0
            lags           = (0:numel(Rhl)-1).';
            J0             = besselj(0, 2 * pi * obj.sParams.dopplerFreq * lags * obj.sParams.Tsymb);
            J0_PSD         = abs(fftshift(fft(J0,Nfft)));
            J0_PSD_dB      = 20*log10(J0_PSD);
            sPSD.J0_PSD_dB = J0_PSD_dB-max(J0_PSD_dB);
            [~, maxInd_J0] = max(sPSD.J0_PSD_dB);
            fd_J0_Hz       = sPSD.fVecHz(maxInd_J0);

            %% sDoppler
            sDoppler.fd_H_Hz  = fd_H_Hz;
            sDoppler.fd_J0_Hz = fd_J0_Hz;

            %% Plot_Channel_Taps_PSD
            if obj.sParams.plotEnable
                obj.Plot_Channel_Taps_PSD(sPSD, sDoppler);
            end
            
            %% Doppler Spread
            dopplerSpread = abs(sDoppler.fd_H_Hz);

        end

    end

    %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
    %% ~~~~~~~~~~~~~~~~~~ Static Methods ~~~~~~~~~~~~~~~~~~~ %%
    %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
    methods(Static)

        %% ~~~~~~~~~~~~~~~~~~~~ Plot_Channel_Taps_PSD ~~~~~~~~~~~~~~~~~~~~ %%
        function [] = Plot_Channel_Taps_PSD(sPSD, sDoppler)

            %% Extract PSD and Doppler Structs
            H_PSD_dB  = sPSD.H_PSD_dB ;
            J0_PSD_dB = sPSD.J0_PSD_dB;

            fVecHz    = sPSD.fVecHz   ;

            fd_H_Hz  = sDoppler.fd_H_Hz;
            fd_J0_Hz = sDoppler.fd_J0_Hz;

            %% Plot Signal
            sPlot.title     = {...
                "Channel Taps PSD",...
                "H: f_d: "+fd_H_Hz+" [Hz]",...
                "J0: f_d: "+fd_J0_Hz+" [Hz]",...
                };
            sPlot.xlabel    = 'Frequency [Hz]';
            sPlot.ylabel    = 'Doppler FFT [dB]';
            sPlot.legend    = {"H","J0"};
            signalToPlot    = horzcat(H_PSD_dB,  J0_PSD_dB);
            cPlot.Plot_signal(fVecHz, signalToPlot ,sPlot);
            hold on;
            plot(fd_H_Hz,0,'bx','MarkerSize',17,'LineWidth',2);
            plot(fd_J0_Hz,0,'rx','MarkerSize',17,'LineWidth',2);

        end

        %% ~~~~~~~~~~~~~~~~~~~~ Calc_tau2 ~~~~~~~~~~~~~~~~~~~~ %%
        function [tau2] = Calc_tau2(tapPhasors, tapDelays)
            tau2 = abs(tapPhasors.').^2*tapDelays.^2./sum(abs(tapPhasors.').^2,2);
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Calc_tau1 ~~~~~~~~~~~~~~~~~~~~ %%
        function [tau1] = Calc_tau1(tapPhasors, tapDelays)
            tau1 = abs(tapPhasors.').^2*tapDelays./sum(abs(tapPhasors.').^2,2);
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Calc_RMSE ~~~~~~~~~~~~~~~~~~~~ %%
        function [sRMSE] = Calc_RMSE(estVal, actualVal)
            err = estVal-actualVal;
            sRMSE.Linear = sqrt(mean(abs(err).^2));
            sRMSE.dB = 10*log10(sRMSE.Linear);
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Calc_Normalized_Error ~~~~~~~~~~~~~~~~~~~~ %%
        function [sErrNorm] = Calc_Normalized_Error(estVal, actualVal)
            errNorm = abs(estVal-actualVal)/abs(actualVal);
            sErrNorm.Linear = errNorm;
            sErrNorm.dB = 10*log10(sErrNorm.Linear);
        end

    end

end
