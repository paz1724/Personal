classdef cModifiedDFT < handle
    properties

        %% Objects
        oEWA                    % object of RMSE

        %% MDFT properties
        deltaOmega              % Omega resolution after the Modified DFT multiplication
        searchAreaAroundInitW   % Area of indices to serch around the init w
        searchAreaCurr          % Current search area for omega that will grow if peak not found
        maxIter                 % Maximum iteration for omega expansion area search
        bKeepSearchingForPeak   % Boolean for keep searhcing for omega peak
        PAPRThdB                % PAPR Threshold for PAPR of omega power vector
        kVecPrev                % Previous omega samples already scanned for complexity optimization

        % Number of hypothesizes: J = (searchAreaAroundInitW+1) * 2 / deltaOmega
    end

    methods

        %% ~~~~~~~~~~~~~~~~~~~~~~ Constructor ~~~~~~~~~~~~~~~~~~~~~~ %%
        function obj = cModifiedDFT(oEWA)

            %% Objects
            obj.oEWA = oEWA;

            %% Initial Values for properties
            obj.deltaOmega            = 1/3;
            obj.searchAreaAroundInitW = 20;
            obj.bKeepSearchingForPeak = true;
            obj.PAPRThdB              = 12;
            obj.maxIter               = 200;
            obj.kVecPrev              = [];

        end

        %% ~~~~~~~~~~~~~~~~~~~~ Apply ~~~~~~~~~~~~~~~~~~~~ %%
        function [sPerfOut] = Apply(obj, y, t)

            %% Phasor Samples
            z = exp(1i*2*pi*y);

            %% Omega_estimation --> O(N*J)
            w_est = obj.Omega_estimation(z, t);

            %% Theta0_estimation --> O(N)
            theta0_est = obj.Theta0_estimation(z, t, w_est);

            %% y prediction based on the estimates --> O(N)
            yPredBest  = mod(w_est.*t.'+theta0_est,1);

            %%  Calc RMSE
            sRMSE = obj.oEWA.Calc_RMSE(w_est, theta0_est, t, y, yPredBest, 'Modified DFT');

            %% Out Struct
            sPerfOut.w_est          = w_est;
            sPerfOut.theta0_est     = theta0_est;
            sPerfOut.w_RMSE         = sRMSE.w;
            sPerfOut.theta0_RMSE    = sRMSE.theta;
            sPerfOut.w_RMSE_dB      = sRMSE.w_dB;
            sPerfOut.theta0_RMSE_dB = sRMSE.theta_dB;

        end

        %% ~~~~~~~~~~~~~~~~~~~~ Omega_estimation ~~~~~~~~~~~~~~~~~~~~ %%
        function [w_est] = Omega_estimation(obj, z, t)

            iter = 0;
            obj.searchAreaCurr = obj.searchAreaAroundInitW;
            while obj.bKeepSearchingForPeak && iter<obj.maxIter
                iter = iter+1;
                %% Generate_modified_DFT
                [MDFTmat, omegaVec] = obj.Generate_modified_DFT(t);

                %% Apply Modified DFT on the recieved signal --> O(N*J)
                ZVec = MDFTmat.'* z;

                %% Frequency Power  --> O(J)
                PVec = abs(ZVec).^2;

                %% PAPR Stop Condition
                PAPRdB = 10*log10(max(PVec)/mean(PVec));
                if obj.oEWA.bShowPlots
                    disp("Iteration: "+iter+" , Omega PAPR: "+PAPRdB+" , Threshold: "+obj.PAPRThdB+" , Search Area: "+obj.searchAreaCurr);
                end
                if PAPRdB>obj.PAPRThdB
                    obj.bKeepSearchingForPeak = false;
                else % Increase the search area of w by multiplication factor
                    obj.searchAreaCurr = obj.searchAreaCurr+obj.searchAreaAroundInitW;
                end
                
            end

            %% Weighted_least_squares_parabolic_interpolation --> O(J)
            w_est = obj.oEWA.Weighted_least_squares_parabolic_interpolation(PVec, omegaVec, obj.deltaOmega);

            %% Plot omega Estimation
            if obj.oEWA.bShowPlots
                obj.Plot_omega_estimation(PVec, w_est, omegaVec);
            end
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Theta0_estimation ~~~~~~~~~~~~~~~~~~~~ %%
        function [theta0_est] = Theta0_estimation(obj, z, t, w_est)

            %% y prediction w/o theta0 based on w_est --> O(N)
            zOnlyTheta0  = z.*exp(-1i*2*pi*w_est.*t);
            theta0EstVec = angle(zOnlyTheta0)/(2*pi);
            sHist        = obj.oEWA.Get_hist_statistics(theta0EstVec, obj.oEWA.histBinsDecFactorInitGuess);

            %% Weighted_least_squares_parabolic_interpolation --> O(N)
            bins        = sHist.bins;
            delta       = sHist.delta;
            f           = sHist.f;
            theta0_est  = obj.oEWA.Weighted_least_squares_parabolic_interpolation(f, bins, delta);

            %% Plot_theta0_estimation
            if obj.oEWA.bShowPlots
                obj.Plot_theta0_estimation(t, theta0EstVec, theta0_est, obj.oEWA.theta0, bins, f);
            end

        end

        %% ~~~~~~~~~~~~~~~~~~~~ Generate_modified_DFT ~~~~~~~~~~~~~~~~~~~~ %%
        function [MDFTmat, omegaVec] = Generate_modified_DFT(obj, t)

            %% New Search area samples vector
            kStart      = floor((obj.oEWA.w_init-obj.searchAreaCurr)/obj.deltaOmega);
            kEnd        = ceil((obj.oEWA.w_init+obj.searchAreaCurr)/obj.deltaOmega);
            kVec        = kStart:kEnd;

            %% Remove omega samples we have already scanned
            kVec         = setdiff(kVec,obj.kVecPrev);

            %% Add the new omega samples to the kVecPrev
            obj.kVecPrev = horzcat(obj.kVecPrev,kVec);

            %% Omega samples vector
            omegaVec    = obj.deltaOmega.*kVec;

            %% Phase matrix
            phaseMat    = 2*pi*t.*omegaVec;

            %% MDFT Matrix
            MDFTmat     = exp(-1i*phaseMat);

        end

        %% ~~~~~~~~~~~~~~~~~~~~ Show_estimation_errors_vs_iterations ~~~~~~~~~~~~~~~~~~~~ %%
        function [] = Show_estimation_errors_vs_iterations(obj, wEstVec, theta0EstVec, minCostVec)
            k               = length(wEstVec);
            %% Plot Final Estimation
            sPlot.title     = 'Estimated Parameters vs time';
            sPlot.xlabel    = 'Iterations';
            sPlot.ylabel    = 'Estimated Parameters';
            sPlot.legend    = {'\omega',"\theta_0","Minimum Cost"};
            sPlot.yline     = [obj.oEWA.w,obj.oEWA.theta0,0];
            signalToPlot    = horzcat(wEstVec,theta0EstVec, minCostVec);
            obj.oEWA.Plot_multiple_signals(1:k, signalToPlot ,sPlot);

        end

        %% ~~~~~~~~~~~~~~~~~~~~ Plot_omega_estimation ~~~~~~~~~~~~~~~~~~~~ %%
        function [] = Plot_omega_estimation(obj, PVec, w_est, omegaVec)
            sPlot.title     = {'Signal after MDFT vs \omega',...
                "\omega Resolution = "+obj.deltaOmega+ " [rad]/[sec]",...
                "True \omega = "+obj.oEWA.w,...
                "Estimated \omega = "+w_est,...
                };
            sPlot.xlabel    = '\omega [rad/sec]';
            sPlot.ylabel    = 'Signal after MDFT [dB]';
            signalToPlot    = 10*log10(PVec);
            obj.oEWA.Plot_signal(omegaVec, signalToPlot ,sPlot);
            xline(obj.oEWA.w,'--g','LineWidth',1.5);
            xline(w_est,'--m','LineWidth',1.5);
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Plot_theta0_estimation ~~~~~~~~~~~~~~~~~~~~ %%
        function [] = Plot_theta0_estimation(obj, t, theta0EstVec, theta0_est, theta0, bins, f)
            
            %% Modulo Protection
            theta0_est    = mod(theta0_est,1);
            bins          = mod(bins,1);
            binsMat       = [bins,f];
            sortedBinsMat = sortrows(binsMat,1);

            %% Signal Plot
            figure;
            subplot(1,2,1);
            plot(t,theta0EstVec,'-o');
            grid minor;
            sPlot.title  = "\theta_0 vs Time [sec]";
            sPlot.xlabel = 'Time [sec]';
            sPlot.ylabel = '\theta_0 [rad/(2*\pi)]';
            title(sPlot.title);
            xlabel(sPlot.xlabel);
            ylabel(sPlot.ylabel);
            yline(theta0_est,'--g','LineWidth',1.5);
            yline(theta0,'--m','LineWidth',1.5);
            legendStr = {"\theta_0 samples","Estimated \theta_0","Actual \theta_0"};
            cPlot.Click_Legend(legendStr);

            %% Histogram Plot
            subplot(1,2,2);
            plot(sortedBinsMat(:,1),sortedBinsMat(:,2),'-o');
            titleStr = {"Histogram of \theta_0 values",...
                "Estimated \theta_0: "+theta0_est+" [rad/(2*\pi)]",...
                "Actual value of \theta_0: "+obj.oEWA.theta0+" [rad/(2*\pi)]"};
            legendStr = {"Histogram","Estimated \theta_0","Actual \theta_0"};
            xline(theta0_est,'--g','LineWidth',1.5);
            xline(obj.oEWA.theta0,'--m','LineWidth',1.5);
            grid minor;
            title(titleStr);
            xlabel("w [rad/sec]");
            ylabel("Histogram");
            cPlot.Click_Legend(legendStr);

            plotbrowser('on');
        end

    end

    methods(Static)

    end
end