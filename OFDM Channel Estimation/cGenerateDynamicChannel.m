classdef cGenerateDynamicChannel < handle
    properties

        %% GDC Properties
        Nsc
        NsymbInSf
        subframeTime_sec
        totalSfnumber
        fsc
        channel
        tapTimeVec_ns
        tapPowerVec_dB
        channelType
        NFFT
        Nsymb
        fm
        scSpacing
        SymbolSpacing
        numRx
        numTx
        corr_name
        timingOffset_us
        freqOffset_Hz
        symbolOffset
        Tsymb
        bShowPlots

        actualChannelDelays
        actualChannelPhasors

    end

    %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
    %% ~~~~~~~~~~~~~~~~~~~~~~ Methods ~~~~~~~~~~~~~~~~~~~~~~ %%
    %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
    methods

        %% ~~~~~~~~~~~~~~~~~~~~~~ Constructor ~~~~~~~~~~~~~~~~~~~~~~ %%
        function obj = cGenerateDynamicChannel(sParams)

            %% Particle Filter
            obj.Nsc              = sParams.Nsc;
            obj.NsymbInSf        = sParams.NsymbInSf;
            obj.subframeTime_sec = sParams.subframeTime_sec;
            obj.totalSfnumber    = sParams.totalSfnumber;
            obj.fsc              = sParams.fsc;
            obj.channel          = sParams.channelType;
            obj.tapTimeVec_ns    = sParams.tapTimeVec_ns;
            obj.tapPowerVec_dB   = sParams.tapPowerVec_dB;
            obj.NFFT             = sParams.Nfft;
            obj.Nsymb            = sParams.Nsymb;
            obj.fm               = sParams.dopplerFreq;
            obj.scSpacing        = sParams.fsc;
            obj.SymbolSpacing    = sParams.symbSpacing;
            obj.numRx            = sParams.numRx;
            obj.numTx            = sParams.nTx;
            obj.corr_name        = sParams.correlationType;
            obj.timingOffset_us  = sParams.timingOffset_us;
            obj.freqOffset_Hz    = sParams.freqOffset_Hz;
            obj.symbolOffset     = sParams.symbOffset;
            obj.Tsymb            = sParams.Tsymb;
            obj.bShowPlots       = sParams.plotEnable;

        end

        %% ~~~~~~~~~~~~~~~~~~~~ Apply ~~~~~~~~~~~~~~~~~~~~ %%
        function [sChannel] = Apply(obj)

            switch obj.channel
                case 'AWGN'
                    sChannel = obj.Calc_channel_AWGN();
                otherwise
                    sChannel = obj.Calc_channel_Rayleigh();
            end

            obj.actualChannelDelays  = sChannel.tapDelays;
            obj.actualChannelPhasors = sChannel.tapPhasors;
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Calc_channel_AWGN ~~~~~~~~~~~~~~~~~~~~ %%
        function sChannel = Calc_channel_AWGN(obj)
            linPhaseSc    = exp(-1j*2*pi*obj.scSpacing*(obj.timingOffset_us*1e-6)*((1:obj.NFFT).'-obj.NFFT/2));
            linPhaseSymb = exp(-1j*2*pi*obj.Tsymb*obj.freqOffset_Hz*(0:obj.Nsymb-1));
            sChannel.Hout = ones(obj.NFFT, obj.Nsymb, obj.numRx, obj.numTx).*linPhaseSc.*linPhaseSymb;
            % repmat(linPhaseSc,1, obj.Nsymb, obj.numRx, obj.numTx).*...
            % repmat(linPhaseSymb,obj.NFFT, 1, obj.numRx, obj.numTx);
            sChannel.Hrb  = sChannel.Hout(1:obj.Nsc,:,:,:);

            % add linear phase in F-domain due to obj.timingOffset
            % Hout = ones(obj.NFFT,obj.Nsymb,obj.numRx,obj.numTx);
            sChannel.tapDelays  = 0;
            sChannel.tapPhasors = 1;
            sChannel.corrInfo   = nan;
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Calc_channel_Rayleigh ~~~~~~~~~~~~~~~~~~~~ %%
        function sChannel = Calc_channel_Rayleigh(obj)

            % channel : either channel name string ('AWGN' / 'ETU'  etc)
            % or tap table: [ tap_time_in_ns   tap_power_in_dB ]
            % BW in Hz (used BW, e.g 18e6 for 20MHz)
            % NFFT               -  sets channel BW length  (e.g. 1024 for 10MHz)
            % Nsymb              -  number of symbols to generate
            % fm                 -  doppler frequency [Hz]
            % SymbolSpacing      -  Time distance between the generated symbols (counted in symbols)
            % numRx, obj.numTx   -  if not specified (or empty), 1x1 obj.channel assumed
            % corr_name          -  'low' / 'med' / 'high'  or [Rx_corr_factor Tx_corr_factor]
            % seed               -   optional, for random number generation
            % timing offset [us] -   optional, for adding a constant phase slope to the obj.channel
            % symbolOffset       - generate channel with offset in time of 'obj.symbolOffset' symbols

            %% Sample Time
            Ts = obj.subframeTime_sec/obj.NsymbInSf;   % Average symbol time

            %% Channel PDP
            if ischar(obj.channel)
                tap_table = obj.Get_channel_tap_table(); %[ ns  dB]
            else
                tap_table = obj.channel;
            end
            num_taps    = length(tap_table(:,2));
            numChannels = obj.numRx*obj.numTx;
            tapDelays    = tap_table(:,1)'*1e-9 + obj.timingOffset_us*1e-6;
            tap_sigma   = 10.^(tap_table(:,2)/20);
            tap_sigma   = tap_sigma / sqrt(sum(abs(tap_sigma).^2)); % noramlize total power

            %% Correlation Matrix
            if exist('obj.corr_name','var') && ~isempty(obj.corr_name)
                [corr, corrInfo] = Calc_corr_matrix(obj.numRx, obj.numTx, obj.corr_name);
            else
                corr     = eye(obj.numRx * obj.numTx);
                corrInfo = [];
            end
            chol_corr = chol(corr);

            %% Time points in which the obj.channel will be computed
            n_reps     = ceil(obj.Nsymb/length(obj.SymbolSpacing));
            spacingVec = repmat(obj.SymbolSpacing,1,n_reps);
            t          = (cumsum(spacingVec) + obj.symbolOffset)*Ts;
            t          = t(1:obj.Nsymb);

            if obj.fm==0 || obj.fm==inf   % extreme cases
                if obj.fm==0
                    randomTaps = 1/sqrt(2) * repmat(rand(num_taps,1,numChannels)  +   1j*randn(num_taps,1,numChannels),1,obj.Nsymb);              % Static obj.channel
                end
                if obj.fm==inf
                    randomTaps =  1/sqrt(2) * (randn(num_taps,obj.Nsymb,numChannels) + 1j*randn(num_taps,obj.Nsymb,numChannels));                      % Taps have no correlation between symbols in time
                end

                tap_mat = bsxfun(@times,tap_sigma, randomTaps);  % multiply each tap by its power

                tapPhasors = mtimes3d(chol_corr' ,permute(tap_mat,[ 3 1 2])); % create correlation between taps of different channels
                tapPhasors = permute(tapPhasors,[2 3 1]);                                                     % permute back to desired dimensions

            else  % 0 < doppler < inf

                %% Calculate scatterers to create doppler spectrum
                fs    = 1/(Ts*min(obj.SymbolSpacing));
                Nsct  = max(round(obj.fm*obj.Nsymb/fs),50);
                Nsct  = min(Nsct,500);
                freqs = obj.fm*cos(2*pi*(1:Nsct)/Nsct).' ;
                FTdop = exp(-1j*2*pi*freqs*t);                                   % Fourier transorm matrix

                phs   = rand(numChannels,num_taps,length(freqs))*2*pi;
                W     = bsxfun(@times,tap_sigma.' , exp(1j*phs)/sqrt(Nsct) );   % random scatteres per tap

                %% Correlate the scatterers and re-order : numTaps X  Nsct X numChannels
                Wcorr = mtimes3d(chol_corr' , W);
                Wcorr = permute(Wcorr, [2 3 1]);

                %% Fourier transform to time domain
                tapPhasors = mtimes3d(Wcorr , FTdop);
            end

            %% Calc_H_Freq
            [Hrb, Hout] = obj.Calc_H_Freq(tapDelays, tapPhasors, obj.freqOffset_Hz);

            %% Out Struct
            sChannel.Hout = Hout;
            sChannel.Hrb = Hrb;
            sChannel.tapDelays  = tapDelays(:);
            sChannel.tapPhasors = tapPhasors;
            sChannel.corrInfo = corrInfo;
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Calc_H_Freq ~~~~~~~~~~~~~~~~~~~~ %%
        function [Hrb, Hout] = Calc_H_Freq(obj, tapDelays, tapPhasors, freqOffset_Hz)

            %% Fourier transform the taps to frequency domain:
            f = (-obj.NFFT/2:obj.NFFT/2-1)'*obj.scSpacing;  % OFDM subcarrier frequencies
            FT = exp(-1j*2*pi*f*tapDelays(:).');       % Fourier transform matrix

            Hout = zeros(obj.NFFT,obj.Nsymb,obj.numRx,obj.numTx);
            for tx=1:obj.numTx
                for rx=1:obj.numRx
                    tap_mat = tapPhasors(:,:,(tx-1)*obj.numRx + rx);
                    Hout(:,:,rx,tx) = FT * tap_mat;
                end
            end

            %% Frequency Offset
            linPhaseSymb = exp(-1j*2*pi*obj.Tsymb*freqOffset_Hz*(0:obj.Nsymb-1));
            Hout = Hout.*linPhaseSymb;

            %% Take relevant subcarriers
            Hout = squeeze(Hout);
            Hout = reshape(Hout,[obj.NFFT,obj.NsymbInSf,obj.totalSfnumber,obj.numTx]);
            Hrb  = Hout(1:obj.Nsc,:,:,:); % cut 1 RB
            Hrb  = reshape(Hrb,size(Hrb,1),[]);
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Get_channel_tap_table ~~~~~~~~~~~~~~~~~~~~ %%
        function tap_table = Get_channel_tap_table(obj)

            % returns a channel's tap table in the format:
            % tap table: [ tap_time_in_ns   tap_power_in_dB ]
            %  TS 36.101  tables B.2.1
            switch obj.channel
                case {'AWGN','EYE','FLAT'}
                    tap_table = [0  0];
                    tap_table(:,1) = tap_table(:,1);
                case {'EVEH_A','EVA'}
                    tap_table = [0 30 150 310 370 710 1090 1730 2510
                        0 -1.5 -1.4 -3.6 -0.6 -9.1 -7 -12 -16.9]';
                case {'EPED_A','EPA'}
                    tap_table = [...
                        0   30  70  90  110     190     410
                        0   -1  -2  -3  -8      -17.2   -20.8]';
                case 'NEW'
                    tap_table = [...
                        0,  7000,   14000,  2000;...
                        0, 7,     -3,      -2]';
                case 'ETU'
                    tap_table = [...
                        0,  50, 120,    200,    230,    500,    1600,   2300,   5000;...
                        -1, -1, -1,     0,      0,      0,     -3,     -5,      -7]';
                case 'ETU_like'
                    tap_table = [...
                        0,  50, 120,    200,    230,    500,    1600,   2300,   5000;...
                        -1, -1, -1,     0,      0,      0,     12,     15,      -7]';
                case 'Two_Taps'
                    tap_table = [...
                        0 4000 ; ...
                        0 0 ]';
                case 'Custom'
                    tap_table = [...
                        obj.tapTimeVec_ns; ...
                        obj.tapPowerVec_dB ]';
            end

        end

        %% ~~~~~~~~~~~~~~~~~~~~ Get_channel_time ~~~~~~~~~~~~~~~~~~~~ %%
        function [sChannelTime] = Get_channel_time(obj, X)

            if ischar(obj.channel)
                obj.channelType = obj.channel;
                maxRelevantT_us = 10;
            else
                maxRelevantT_us = max(obj.channel(:,1))+5;
                obj.channelType = 'Custom';
            end

            hTime        = (obj.NFFT/obj.Nsc)*ifft(X,obj.NFFT,1);
            hTimeAvg     = mean(hTime,[2,3]);
            hTimeAvg_dB  = db(hTimeAvg);
            hTime        = reshape(hTime,obj.NFFT,[]);
            hTime_dB     = db(hTime);
            tVec         = (0:obj.NFFT-1)./(obj.NFFT*obj.fsc);
            tVec_us      = tVec*1e6;
            relevantInds = tVec_us<=maxRelevantT_us;
            tVec_us      = tVec_us(relevantInds);
            hTime_dB     = hTime_dB(relevantInds,:);
            hTimeAvg_dB  = hTimeAvg_dB(relevantInds,:);

            %% Channel Time Struct
            sChannelTime.tVec_us     = tVec_us;
            sChannelTime.hTime_dB    = hTime_dB;
            sChannelTime.hTimeAvg_dB = hTimeAvg_dB;

        end

        %% ~~~~~~~~~~~~~~~~~~~~ %% Plot_channel_time ~~~~~~~~~~~~~~~~~~~~ %%
        function [] = Plot_channel_time(obj, sChannelTime, sParams)

            if obj.bShowPlots
                tVec_us     = sChannelTime.tVec_us;
                hTimeAvg_dB = sChannelTime.hTimeAvg_dB;

                Create_figure(sParams.graphAxes);
                plot(tVec_us,hTimeAvg_dB);
                xlabel('Time [us]');
                ylabel('Channel Impulse Response [dB]');
                title({'Time Domain',"Channel Impulse Response Vs Time",...
                    "Channel: "+strrep(obj.channelType,'_',' ')+...
                    ",  Doppler: "+obj.fm+" [Hz]"+...
                    ",  SNR: "+sParams.SNR_dB+" [dB]"});
                if isempty(sParams.graphAxes)
                    plotbrowser on;
                end
            end
        end

        %% ~~~~~~~~~~~~~~~~~~~~ %% Show_image_channel_fast_slow_time ~~~~~~~~~~~~~~~~~~~~ %%
        function [] = Show_image_channel_fast_slow_time(obj, sChannelTime)

            if obj.bShowPlots
                %% Channel Time Struct
                tVec_us     = sChannelTime.tVec_us;
                hTime_dB    = sChannelTime.hTime_dB;

                figure;
                Tdecimated = 50;
                T = length(tVec_us);
                M = round(T/Tdecimated);
                tDec_us = tVec_us(1:M:end);

                imagesc(hTime_dB); colorbar;
                set(gca, 'XTick', 1:size(hTime_dB,2), 'XTickLabel', 1:size(hTime_dB,2),'XTickLabelRotation',90) % 10 ticks
                set(gca, 'YTick', 1:M:size(hTime_dB,1), 'YTickLabel', tDec_us) % 20 ticks
                ylabel('Time [us]');
                xlabel('Symbols');
                title({'Time Domain',"Channel Impulse Response Vs Time",...
                    "Channel: "+obj.channelType+",  Doppler: "+obj.fm+" [Hz]"});
                plotbrowser on;
            end
        end

    end

    %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
    %% ~~~~~~~~~~~~~~~~~~ Static Methods ~~~~~~~~~~~~~~~~~~~ %%
    %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
    methods(Static)

        %% ~~~~~~~~~~~~~~~~~~~~ %% Debug: Show_channel_phase ~~~~~~~~~~~~~~~~~~~~ %%
        function [] = Show_channel_phase(Hrb, subcarriersToPlot, symbolsToPlot)

            sLegend.commonString = 'Subcarriers To Plot';
            sLegend.values       = subcarriersToPlot;
            sLegend.units        = [];

            H = Hrb(subcarriersToPlot,:,:);
            H = reshape(H,size(H,1),[]).';
            figure;
            subplot(2,2,1);
            plot(angle(H));
            grid on; grid minor;
            title({'Frequency Domain','Channel Phase Vs Symbol'});
            ylabel('Channel Phase [rad]');
            xlabel('Symbol Index');
            cPlot.Click_Legend_List(sLegend);

            subplot(2,2,2);
            plot(unwrap(angle(H)));
            grid on; grid minor;
            title({'Frequency Domain','Channel Unwrapped Phase Vs Symbol'});
            ylabel('Channel Phase [rad]');
            xlabel('Symbol Index');
            cPlot.Click_Legend_List(sLegend);

            sLegend.commonString = 'Symbols To Plot';
            sLegend.values       = symbolsToPlot;
            sLegend.units        = [];

            H = reshape(Hrb,size(H,1),[]).';
            H = H(:,symbolsToPlot);
            subplot(2,2,3);
            plot(angle(H));
            grid on; grid minor;
            title({'Frequency Domain','Channel Phase Vs Sub-Carrier'});
            ylabel('Channel Phase [rad]');
            xlabel('Sub-Crrier Index');
            cPlot.Click_Legend_List(sLegend);

            subplot(2,2,4);
            plot(unwrap(angle(H)));
            grid on; grid minor;
            title({'Frequency Domain','Channel Unwrapped Phase Vs Sub-Crrier'});
            ylabel('Channel Phase [unwrapped rad]');
            xlabel('Sub-Crrier Index');
            cPlot.Click_Legend_List(sLegend);

            plotbrowser('on');
        end

        %% ~~~~~~~~~~~~~~~~~~~~ %% Debug: Show_channel_power ~~~~~~~~~~~~~~~~~~~~ %%
        function [] = Show_channel_power(Hrb, subcarriersToPlot, symbolsToPlot)

            H = Hrb(subcarriersToPlot,:,:);
            H = reshape(H,size(H,1),[]).';
            figure;
            subplot(1,2,1);
            plot(20*log10(abs(H)));
            grid on; grid minor;
            title({'Frequency Domain','Channel Power Vs Symbol'});
            ylabel('Channel Power [dB]');
            xlabel('Symbol Index');

            sLegend.commonString = 'Subcarriers To Plot';
            sLegend.values       = subcarriersToPlot;
            sLegend.units        = [];
            cPlot.Click_Legend_List(sLegend);

            H = reshape(Hrb,size(H,1),[]).';
            H = H(:,symbolsToPlot);
            subplot(1,2,2);
            plot(20*log10(abs(H)));
            grid on; grid minor;
            title({'Frequency Domain','Channel Power Vs Sub-Carrier'});
            ylabel('Channel Power [dB]');
            xlabel('Sub-Crrier Index');

            sLegend.commonString = 'Symbols To Plot';
            sLegend.values       = symbolsToPlot;
            sLegend.units        = [];
            cPlot.Click_Legend_List(sLegend);

            plotbrowser('on');
        end

    end

end
