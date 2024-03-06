classdef cEstimateWheelAngle < handle
    properties
        w                           % omega: the slope of theta
        theta0                      % the initial phase of theta
    
        sigma_n                     % the std of the received signal's noise
        pUniform                    % the probability of uniform sample in the received signal
        N                           % number of received samples
    
        w_init                      % Initial guess for w
        theta0_init                 % Initial guess for theta0
        w_init_std                  % Initial std guess for w

        outlierPercentTh            % The percentage Threshold in the CDF of the differential of theta to regard as an outlier
        histBinsDecFactorInitGuess  % Decimation factor for hist bins division for initil guess for w
        bShowPlots                  % Boolean for showing debug plots

    end

    methods

        %% ~~~~~~~~~~~~~~~~~~~~~~ Constructor ~~~~~~~~~~~~~~~~~~~~~~ %%
        function obj = cEstimateWheelAngle(w, theta0, sigma_n, pUniform, N, bShowPlots)
            obj.w                          = w;
            obj.theta0                     = theta0;

            obj.sigma_n                    = sigma_n;
            obj.pUniform                   = pUniform;
            obj.N                          = N;

            obj.outlierPercentTh           = 70;
            obj.histBinsDecFactorInitGuess = 30;
            obj.bShowPlots                 = bShowPlots;

        end

        %% ~~~~~~~~~~~~~~~~~~~~~~ Generate Input Samples ~~~~~~~~~~~~~~~~~~~~~~ %%
        function [y, t] = GenerateInputSamples(obj)
            %% Generate uniform-spaced time samples
            t                = vertcat(0,sort(rand(obj.N-1,1)));

            %% Generate theta samples
            nVec             = obj.sigma_n*randn(obj.N,1);
            theta            = mod(obj.w*t+obj.theta0+nVec,1);

            %% Received Signal: y
            numPointsUniform = round(obj.N*obj.pUniform);
            uniformInds      = randperm(obj.N, numPointsUniform);
            y                = theta;
            y(uniformInds)   = rand(numPointsUniform,1);

            %% Plot Signal
            if obj.bShowPlots
                sPlot.title     = {
                "Signal vs Time",...
                "w: "+obj.w,...
                "\theta_0: "+obj.theta0,...
                "\sigma_n: "+obj.sigma_n,...
                "pUniform: "+obj.pUniform,...
                "N Samples: "+obj.N,...
                };
                sPlot.xlabel    = 'time [sec]';
                sPlot.ylabel    = 'Signal [rad*\pi]';
                sPlot.legend    = {'\theta','y'};
                obj.Plot_signal(t, horzcat(theta,y), sPlot);
            end
        end

        %% ~~~~~~~~~~~~~~~~~~~~~~ Generate Init Guess: W , theta0 ~~~~~~~~~~~~~~~~~~~~~~ %%
        function [] = GenerateInitGuess(obj, y, t)

            %% Differential
            tDiffVec            = t(2:end)-t(1:end-1);
            diffY               = (y(2:end)-y(1:end-1))./tDiffVec;

            %% Get Statistics of the differential
            sStatistics         = obj.Get_statistics(diffY);
            
            %% Remove Outliers higher than 70% CDF
            diffYForPlot              = diffY;
            outlierInds               = abs(diffYForPlot)>sStatistics.cdfPercentOutlier;
            diffYForPlot(outlierInds) = [];

            %% Take max bin in histogram of differentials for initial guess
            sHist           = obj.Get_hist_statistics(diffYForPlot, obj.histBinsDecFactorInitGuess);
            obj.w_init      = obj.Weighted_least_squares_parabolic_interpolation(sHist.f, sHist.bins, sHist.delta);
            obj.w_init_std  = sStatistics.std;
            obj.theta0_init = y(1);

            %% Plot differential Y
            if obj.bShowPlots
                obj.Plot_differential_Y(diffY, sStatistics, t, sHist);
            end
        end

        %% ~~~~~~~~~~~~~~~~~~~~~~ Calc_RMSE ~~~~~~~~~~~~~~~~~~~~~~ %%
        function [sRMSE] = Calc_RMSE(obj, w_est, theta0_est, t, y, yPredBest, method)
            %%  RMSE
            sRMSE = obj.Calc_RMSE_Error(w_est, theta0_est);

            %% Show results
            if obj.bShowPlots
                RMSE_w_title      = "w: True: " + obj.w + " , Est: " + w_est + " , RMSE: = "+sRMSE.w_dB+" [dB]";
                RMSE_theta0_title = "\theta_0: True: " + obj.theta0 + " , Est: " + theta0_est + " , RMSE:  = "+sRMSE.theta_dB+" [dB]";

                %% Plot Final Estimation
                sPlot.title     = {method,'Signal vs Time',RMSE_w_title,RMSE_theta0_title};
                sPlot.xlabel    = 'time [sec]';
                sPlot.ylabel    = 'signal [rad*\pi]';
                sPlot.legend    = {'y_{meas}',"y_{pred}"};
                signalToPlot    = horzcat(y,yPredBest.');
                obj.Plot_signal(t, signalToPlot ,sPlot);
            end
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Cost_function ~~~~~~~~~~~~~~~~~~~~ %%
        function J = Cost_function(obj, w, theta0, y, t)
            yPred        = mod(w .* t + theta0, 1);
            dist         = (y - yPred).^2;
            distReshaped = reshape(dist,size(dist,1),[]);
            J            = obj.Non_outlier_average(distReshaped.');
            J            = reshape(J,size(dist,2),size(dist,3));

            %% Plot Final Estimation
            bPlot = false;
            if bPlot
                sPlot.title     = 'Signal vs Time';
                sPlot.xlabel    = 'time [sec]';
                sPlot.ylabel    = 'signal [rad*\pi]';
                sPlot.legend    = {'y_{meas}',"y_{pred}"};
                signalToPlot    = horzcat(y,yPred(:,2,24));
                obj.Plot_signal(t, signalToPlot ,sPlot);
            end

        end

        %% ~~~~~~~~~~~~~~~~~~~~ Non_outlier_average ~~~~~~~~~~~~~~~~~~~~ %%
        function noOutlierDist = Non_outlier_average(obj, dist)
            distCdfPercent = obj.Get_percentage_value_from_vector(dist.',obj.outlierPercentTh).';
            noOutlierDist = dist;
            noOutlierDist(bsxfun(@gt,noOutlierDist,distCdfPercent)) = nan;
            noOutlierDist = mean(noOutlierDist,2,"omitnan");
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Get_statistics ~~~~~~~~~~~~~~~~~~~~ %%
        function [sStatistics] = Get_statistics(obj, signal)
            sStatistics.cdfPercentOutlier = obj.Get_percentage_value_from_vector(abs(signal),obj.outlierPercentTh);
            sStatistics.median            = median(signal);
            % sStatistics.std               = std(signal);
            validInds                     = abs(signal)<sStatistics.cdfPercentOutlier;
            sStatistics.std               = std(signal(validInds));
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Plot_differential_Y ~~~~~~~~~~~~~~~~~~~~ %%
        function [] = Plot_differential_Y(obj, diffY, sStatistics, t, sHist)
            diffYForPlot = diffY;
            outlierInds = abs(diffYForPlot)>sStatistics.cdfPercentOutlier;
            diffYForPlot(outlierInds) = nan;
            if isnan(sStatistics.std)
                sStatistics.std = 0;
            end

            titleStr        = {...
                'd(\theta)/dt vs time',...
                "Median slope = " + sStatistics.median,...
                "Std slope    = " + sStatistics.std,...
                "CDF "+obj.outlierPercentTh+"% = " + sStatistics.cdfPercentOutlier,...
                };

            sPlot.title  = titleStr;
            sPlot.xlabel = 'time [sec]';
            sPlot.ylabel = 'd(\theta)/dt [\pi*rad/sec]';

            figure;
            %% Plot Signal
            subplot(1,2,1);
            plot(t(2:end),diffYForPlot,'-o');
            grid minor;
            title(sPlot.title);
            xlabel(sPlot.xlabel);
            ylabel(sPlot.ylabel);
            yline(obj.w_init,'--g','LineWidth',1.5);
            yline(obj.w,'--m','LineWidth',1.5);
            yline(sStatistics.median+sStatistics.std,'--c','LineWidth',1.5);
            yline(sStatistics.median-sStatistics.std,'--c','LineWidth',1.5);
            legendStr = {"\omega samples",...
                "Initial Guess for w",...
                "Actual Value of w",...
                "Upper Std",...
                "Lower Std",...
                };
            cPlot.Click_Legend(legendStr);
            
            %% Plot Histogram
            subplot(1,2,2);
            bins = sHist.bins;
            f    = sHist.f;
            plot(bins,f,'-o');
            legendStr = {"Histogram",
                "Initial Guess for w",
                "Actual Value of w"};
            xline(obj.w_init,'--g','LineWidth',1.5);
            xline(obj.w,'--m','LineWidth',1.5);
            grid minor;
            title({"Histogram of slope values", ...
                "Initial Guess for w: "+obj.w_init+" [rad/sec]", ...
                "Actual value of w: "+obj.w+" [rad/sec]"});
            xlabel("w [rad/sec]");
            ylabel("Histogram");
            cPlot.Click_Legend(legendStr);

            plotbrowser('on');

        end

        %% ~~~~~~~~~~~~~~~~~~~~ Calc_RMSE_Error ~~~~~~~~~~~~~~~~~~~~ %%
        function [sRMSE] = Calc_RMSE_Error(obj, w_est, theta0_est)
            wErr             = w_est-obj.w;
            thetaErr         = theta0_est-obj.theta0;
            thetaErr         = obj.Protection_from_modulu_for_theta(thetaErr);

            M               = length(wErr);
            wErr            = wErr(max(1,ceil(M*0.2)):M);
            thetaErr        = thetaErr(max(1,ceil(M*0.2)):M);

            sRMSE.w        = sqrt(mean(abs(wErr).^2));
            sRMSE.theta    = sqrt(mean(abs(thetaErr).^2));
            sRMSE.w_dB     = 10*log10(sRMSE.w);
            sRMSE.theta_dB = 10*log10(sRMSE.theta);
            
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Show_estimation_errors_vs_iterations ~~~~~~~~~~~~~~~~~~~~ %%
        function [] = Show_estimation_errors_vs_iterations(obj, t, sObjectives)

            %% sObjectives
            wEstVec         = sObjectives.wEstVec;
            theta0EstVec    = sObjectives.theta0EstVec;
            estErr          = sObjectives.estErr;

            %% Calculate RMSE
            sRMSE = obj.Calc_RMSE_Error(wEstVec, theta0EstVec);

            %% Plot Final Estimation
            sPlot.title     = 'Estimated Parameters vs time';
            sPlot.xlabel    = 'time [sec]';
            sPlot.ylabel    = ["w" , "\theta_0", "Estimation Error"];
            sPlot.legend    = {...
                "w RMSE: "+sRMSE.w_dB+ " [dB]", ...
                "\theta_0 RMSE: "+sRMSE.theta_dB+ " [dB]", ...
                "Estimation Error"};
            sPlot.yline     = [obj.w,obj.theta0,0];
            signalToPlot    = horzcat(wEstVec, theta0EstVec, estErr);
            obj.Plot_multiple_signals(t, signalToPlot ,sPlot);

        end

    end

    methods(Static)

        %% ~~~~~~~~~~~~~~~~~~~~ Get_hist_statistics ~~~~~~~~~~~~~~~~~~~~ %%
        function [sHist] = Get_hist_statistics(signal, decFactor)
            [histbins,edges] = histcounts(signal, ceil(length(signal)/decFactor));
            bins             = (edges(1:end-1)+edges(2:end)).'/2;
            f                = histbins.'/sum(histbins);
            [~,maxInd]       = max(f);
            maxHist          = bins(maxInd);
            delta            = mean(diff(bins));

            %% Outpuit Struct
            sHist.maxHist = maxHist;
            sHist.bins    = bins;
            sHist.delta   = delta;
            sHist.f       = f;
            sHist.maxInd  = maxInd;

        end

        %% ~~~~~~~~~~~~~~~~~~~~ Plot_signal ~~~~~~~~~~~~~~~~~~~~ %%
        function [] = Plot_signal(xVec, signal, sPlot)

            %% Default
            sPlot = Set_default_value(sPlot,'title','Signal Vs Time');
            sPlot = Set_default_value(sPlot,'xlabel','Time');
            sPlot = Set_default_value(sPlot,'ylabel','Signal');
            sPlot = Set_default_value(sPlot,'legend',[]);

            %% Plot Struct
            titleStr     = sPlot.title;
            xlabelStr    = sPlot.xlabel;
            ylabelStr    = sPlot.ylabel;
            legendStr    = sPlot.legend;

            %% Plot
            figure;
            plot(xVec,signal,'-o');
            grid minor;
            title(titleStr);
            xlabel(xlabelStr);
            ylabel(ylabelStr);
            cPlot.Click_Legend(legendStr);
            plotbrowser('on');
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Plot_multiple_signals ~~~~~~~~~~~~~~~~~~~~ %%
        function [] = Plot_multiple_signals(xVec, signals, sPlot)

            %% Default
            sPlot = Set_default_value(sPlot,'title','Signal Vs Time');
            sPlot = Set_default_value(sPlot,'xlabel','Time');
            sPlot = Set_default_value(sPlot,'ylabel','Signal');
            sPlot = Set_default_value(sPlot,'legend',[]);
            sPlot = Set_default_value(sPlot,'yline',[]);

            %% Plot Struct
            titleStr     = sPlot.title;
            xlabelStr    = sPlot.xlabel;
            ylabelStr    = sPlot.ylabel;
            legendStr    = sPlot.legend;
            ylineVec     = sPlot.yline;

            %% Plot
            figure;
            h = [];
            M = size(signals,2);
            for i=1:M
                h(i) = subplot(M,1,i);
                plot(xVec,signals(:,i),'-o');
                grid minor;
                title({titleStr,legendStr{i}});
                xlabel(xlabelStr);
                ylabel(ylabelStr(i));
                if ~isempty(ylineVec)
                    yline(ylineVec(i),'--g','LineWidth',2);
                end
            end
            plotbrowser('on');
            linkaxes(h,'x');

        end
        
        %% ~~~~~~~~~~~~~~~~~~~~ Set_default_value ~~~~~~~~~~~~~~~~~~~~ %%
        function [paramStruct] = Set_default_value(paramStruct,fieldname,defaultVal)

            if ~isfield(paramStruct,fieldname) || isempty(paramStruct.(fieldname))
                paramStruct.(fieldname) 	= defaultVal;
            end
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Get_percentage_value_from_vector ~~~~~~~~~~~~~~~~~~~~ %%
        function [percentValVec] = Get_percentage_value_from_vector(inputMat,percentage)

            %%
            percentValVec           = zeros(1,size(inputMat,2));
            for i = 1:size(inputMat,2)
                inputVec          	= inputMat(:,i);
                percentVal          = cEstimateWheelAngle.Get_cdf_percent(inputVec,percentage);
                percentValVec(i)    = percentVal;
            end
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Get_cdf_percent ~~~~~~~~~~~~~~~~~~~~ %%
        function [percentVal] = Get_cdf_percent(inputVec,percentage)

            L              = length(inputVec);
            indVec         = 1:L;
            sortedInputVec = sort(inputVec);
            effIndex       = percentage*L/100;
            percentVal     = cEstimateWheelAngle.Linear_interpolation_extrapolation(indVec, sortedInputVec, effIndex);
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Linear_interpolation_extrapolation ~~~~~~~~~~~~~~~~~~~~ %%
        function [relYVec] = Linear_interpolation_extrapolation(xVec, yVec, x0Vec)

            relYVec = zeros(size(x0Vec));
            for xInd = 1:length(x0Vec)
                x0            = x0Vec(xInd);
                relY          = cEstimateWheelAngle.Linear_interpolation_extrapolation_core(xVec, yVec, x0);
                relYVec(xInd) = relY;
            end
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Linear_interpolation_extrapolation_core ~~~~~~~~~~~~~~~~~~~~ %%
        function [relY] = Linear_interpolation_extrapolation_core(xVec, yVec, x0)

            %% Divide to 2 halves to determine if Increasing or Decreasing
            L                   = length(xVec);
            if L==1
                relY = yVec;
            else
                firstHalfAvg      	= mean(xVec(1:floor(L/2)));
                lastHalfAvg        	= mean(xVec(floor(L/2)+1:L));

                %% Cross Indices
                if firstHalfAvg<lastHalfAvg
                    %% Increasing
                    indBefore         	= find(xVec < x0, 1, 'last');
                    indAfter            = find(xVec >= x0, 1, 'first');
                else
                    %% Decreasing
                    indBefore         	= find(xVec > x0, 1, 'last');
                    indAfter            = find(xVec <= x0, 1, 'first');
                end

                %% x y Terms
                y1                  = xVec(indBefore);
                x1                  = yVec(indBefore);
                y2                  = xVec(indAfter);
                x2                  = yVec(indAfter);

                %% Cases
                if ~isempty(y1) && ~isempty(y2)
                    if y2-y1~=0
                        relY        = x1 + (x2 - x1) / (y2 - y1) * (x0 - y1);
                    else
                        relY        = (x1 + x2)/2;
                    end
                elseif ~isempty(y1) && isempty(y2)
                    y2              = y1;
                    x2              = x1;
                    y1              = xVec(indBefore-1);
                    x1              = yVec(indBefore-1);
                    if ~isempty(y1) && ~isempty(y2)
                        if y2-y1~=0
                            relY  	= x1 + (x2 - x1) / (y2 - y1) * (x0 - y1);
                        else
                            relY 	= (x1 + x2)/2;
                        end
                    else
                        relY        = nan;
                    end
                elseif isempty(y1) && ~isempty(y2)
                    y1              = y2;
                    x1              = x2;
                    y2              = xVec(indAfter+1);
                    x2              = yVec(indAfter+1);
                    if ~isempty(y1) && ~isempty(y2)
                        if y2-y1~=0
                            relY  	= x1 + (x2 - x1) / (y2 - y1) * (x0 - y1);
                        else
                            relY  	= (x1 + x2)/2;
                        end
                    else
                        relY        = nan;
                    end
                else
                    relY            = nan;
                end
            end
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Protection_from_modulu_for_theta ~~~~~~~~~~~~~~~~~~~~ %%
        function [thetaErr] = Protection_from_modulu_for_theta(thetaErr)
            if thetaErr>0.5
                thetaErr = thetaErr-1;
            elseif thetaErr<-0.5
                thetaErr = thetaErr+1;
            end
        end
    
        %% Get_polinom_matrix_and_near_peak_z_indices --> O(J)
        function [X, indP] = Get_polinom_matrix_and_near_peak_P_indices(PVec)
            [~,maxIndP] = max(PVec);

            if maxIndP==1
                sideMode = "Right";
            elseif maxIndP==length(PVec)
                sideMode = "Left";
            else
                sideMode = "Center";
            end

            switch sideMode
                case "Center"
                    indP = maxIndP-1:maxIndP+1;
                    X = (-1:1).'.^(0:2);
                case "Left"
                    indP = maxIndP-2:maxIndP;
                    X = (-2:0).'.^(0:2);
                case "Right"
                    indP = maxIndP:maxIndP+2;
                    X = (0:2).'.^(0:2);
            end
        end

        %% Weighted Least Squares Parabolic Interpolation --> O(J)
        function [val_est] = Weighted_least_squares_parabolic_interpolation(PVec, xVec, delta)

            %% Get_polinom_matrix_and_near_peak_P_indices --> O(J)
            [X, indP] = cEstimateWheelAngle.Get_polinom_matrix_and_near_peak_P_indices(PVec);

            %% P Vector
            valP = PVec(indP);

            %% Weights are the P values arounf the max
            W = diag(valP);

            %% Weighted Least Squares Solution --> O(1)
            alphaVec = (X'*W*X)\X'*W*valP;

            %% Parabolic Interpolation --> O(1)
            maxX = -alphaVec(2)./(2*alphaVec(3));

            %% Transform maxX to maxW --> O(1)
            maxIndP = indP(2);
            val_est = xVec(maxIndP)+maxX.*delta;

        end
    end
        
end