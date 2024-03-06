classdef cGridSearch < handle
    properties

        %% Objects
        oEWA                    % object of RMSE

        %% NORMSEOS properties
        N_scan_w                % Scan number of points over w
        N_scan_theta0           % Scan number of points over theta0
        maxIterationsMF         % Max iterations until convergence
        minCostDiffTh           % Threshold of target cost difference convergence
        minCostTh               % Threshold of target cost convergence
        w_init_area_scale       % Scale of Initial Area of scan of w regard w_est
        w_init_area             % Initial Area of scan of w
        centerJumpFactor_w      % Factor of search area to jump the center to when the minimum is in the edge
        zoomFactor_w            % Zoom factor for w search area when the minimum is inside the area
        zoomFactor_theta0       % Zoom factor for theta0 search area when the minimum is inside the area
    end

    methods

        %% ~~~~~~~~~~~~~~~~~~~~~~ Constructor ~~~~~~~~~~~~~~~~~~~~~~ %%
        function obj = cGridSearch(oEWA)

            %% Objects
            obj.oEWA             = oEWA;

            %% Initial Values for properties
            obj.N_scan_w           = 16;        % Scan number of points over w -->Jw
            obj.N_scan_theta0      = 32;        % Scan number of points over theta0 -->Jt
            obj.maxIterationsMF    = 10;        % Max iterations until convergence
            obj.minCostDiffTh      = 1e-5;      % Threshold of target cost difference convergence
            obj.minCostTh          = 1e-3;      % Threshold of target cost convergence
            obj.w_init_area_scale  = 0.2;       % Scale of Initial Area of scan of w regard w_est
            obj.w_init_area        = 20;        % Initial Area of scan of w
            obj.centerJumpFactor_w = 1.5;       % Factor of search area to jump the center to when the minimum is in the edge
            obj.zoomFactor_w       = 0.4;       % Zoom factor for w search area when the minimum is inside the area
            obj.zoomFactor_theta0  = 0.8;       % Zoom factor for theta0 search area when the minimum is inside the area

        end

        %% ~~~~~~~~~~~~~~~~~~~~ Apply ~~~~~~~~~~~~~~~~~~~~ %%
        function [sPerfOut] = Apply(obj, y, t)

            %% Init
            w_est           = obj.oEWA.w_init;
            w_center        = w_est;
            w_area          = obj.w_init_area;% max([obj.oEWA.w_init_std/2, w_est*obj.w_init_area_scale,obj.w_init_area]);
            theta0_est      = obj.oEWA.theta0_init;
            theta0_center   = theta0_est;
            theta0_area     = 1;
            minCostPrev     = inf;
            bKeepSearching  = true;
            iter            = 0;

            %% Init Estimation Vectors
            wEstVec         = zeros(obj.maxIterationsMF,1);
            theta0EstVec    = zeros(obj.maxIterationsMF,1);
            minCostVec      = zeros(obj.maxIterationsMF,1);

            %% Search Space
            w_pred_vec = linspace( ...
                w_center-w_area, ...
                w_center+w_area, ...
                obj.N_scan_w);

            theta0_pred_vec = linspace( ...
                0, ...
                1, ...
                obj.N_scan_theta0);
            theta0_pred_vec = permute(theta0_pred_vec,[1 3 2]);

            %% Iteration Loop --> O(Niter*N*J) (J=jw*Jt)
            while bKeepSearching
                iter = iter+1;

                %% Show iterations
                if obj.oEWA.bShowPlots
                    disp("Iteration: "+iter);
                end

                %% Cost_function --> O(N*J)
                costValueMat = obj.oEWA.Cost_function(w_pred_vec, theta0_pred_vec, y, t);

                %% Find Minimum Value and Index --> O(J)
                minCost                  = min(costValueMat(:));
                minCostDiff              = abs(minCostPrev-minCost);

                %% Find minimum cost
                [minInd.w,minInd.theta0] = ind2sub(size(costValueMat),find(costValueMat == minCost,1,'first'));
                minCostPrev              = minCost;
                w_est                    = w_pred_vec(minInd.w);
                theta0_est               = theta0_pred_vec(minInd.theta0);

                %% Show Image for Debug
                if obj.oEWA.bShowPlots
                    obj.Debug_grid_scan_image(costValueMat, theta0_pred_vec, w_pred_vec, minInd, minCost, minCostDiff, iter, y, t);
                end

                %% Stop Condition or update search space
                if (minCostDiff<obj.minCostDiffTh  && minCost<obj.minCostTh) || iter==obj.maxIterationsMF
                    bKeepSearching = false;
                else
                    %% Adjust w search space
                    [w_pred_vec, w_center, w_area] = obj.Adjust_search_space_w( ...
                        w_pred_vec, ...
                        w_center, ...
                        w_area, ...
                        minInd.w);

                    %% Adjust theta0 search space
                    [theta0_pred_vec, theta0_center, theta0_area] = obj.Adjust_search_space_theta0( ...
                        theta0_pred_vec, ...
                        theta0_center, ...
                        theta0_area, ...
                        minInd.theta0);

                end

                %% Update Estimation Vectors
                wEstVec(iter)       = w_est;
                theta0EstVec(iter)  = theta0_est;
                minCostVec(iter)    = minCost;

            end

            %% Truncate estimation vectors up to current iteration
            wEstVec         = wEstVec(1:iter);
            theta0EstVec    = theta0EstVec(1:iter);
            minCostVec      = minCostVec(1:iter);

            %% Show_estimation_errors_vs_iterations
            if obj.oEWA.bShowPlots
                obj.Show_estimation_errors_vs_iterations(wEstVec, theta0EstVec, minCostVec);
            end
            
            %% y prediction based on the estimates
            yPredBest  = mod(w_est.*t.'+theta0_est,1);

            %%  Calc RMSE
            sRMSE = obj.oEWA.Calc_RMSE(w_est, theta0_est, t, y, yPredBest, 'No Outlier RMSE Optimization Search');
 
            %% Out Struct
            sPerfOut.w_est          = w_est;
            sPerfOut.theta0_est     = theta0_est;
            sPerfOut.w_RMSE         = sRMSE.w;
            sPerfOut.theta0_RMSE    = sRMSE.theta;
            sPerfOut.w_RMSE_dB      = sRMSE.w_dB;
            sPerfOut.theta0_RMSE_dB = sRMSE.theta_dB;

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

        %% ~~~~~~~~~~~~~~~~~~~~ Debug_grid_scan_image ~~~~~~~~~~~~~~~~~~~~ %%
        function [] = Debug_grid_scan_image(obj, costValueMat, theta0_pred_vec, w_pred_vec, minInd, minCost, minCostDiff, iter, y, t)
            figure;
            J = 10*log10(costValueMat);
            imagesc(J); colorbar;
            set(gca, 'XTick', 1:size(J,2), 'XTickLabel', theta0_pred_vec(:),'XTickLabelRotation',90) % 10 ticks
            set(gca, 'YTick', 1:size(J,1), 'YTickLabel', w_pred_vec) % 20 ticks
            xlabel('\theta_0');
            ylabel('\omega');
            title({"Iteration: "+iter,...
                "Cost Value [dB] Vs w , \theta_0",...
                "\omega^* = "+obj.oEWA.w+" , \theta_0^* = "+obj.oEWA.theta0,...
                "\omega_{Est} = "+w_pred_vec(minInd.w)+" , \theta_{Est} = "+theta0_pred_vec(minInd.theta0),...
                "Minimum Cost: "+10*log10(minCost)+" [dB]",...
                "Minimum Cost Diff: "+10*log10(minCostDiff)+" [dB]",...
                });

            hold on;
            plot(minInd.theta0,minInd.w,'ro');
            plotbrowser('on');

            %% Plot Final Estimation
            bPlot = false;
            if bPlot
                theta0Ind   = 31;
                wInd        = 45;
                obj.Debug_cost_matrix(wInd, theta0Ind, t, y, w_pred_vec, theta0_pred_vec);
            end

        end

        %% ~~~~~~~~~~~~~~~~~~~~ Adjust_search_space_w ~~~~~~~~~~~~~~~~~~~~ %%
        function [predVec, centerVal, searchArea] = Adjust_search_space_w(obj, predVec, centerVal, searchArea, minInd)

            if minInd==1 
                predVec = predVec - obj.centerJumpFactor_w*searchArea;
                centerVal = centerVal - obj.centerJumpFactor_w*searchArea;
            elseif minInd==obj.N_scan_w
                predVec = predVec + obj.centerJumpFactor_w*searchArea;
                centerVal = centerVal + obj.centerJumpFactor_w*searchArea;
            else
                centerVal      = predVec(minInd);
                searchArea     = searchArea*obj.zoomFactor_w;
                predVec = linspace( ...
                    centerVal-searchArea, ...
                    centerVal+searchArea, ...
                    obj.N_scan_w);  
            end
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Adjust_search_space_theta0 ~~~~~~~~~~~~~~~~~~~~ %%
        function [predVec, centerVal, searchArea] = Adjust_search_space_theta0(obj, predVec, centerVal, searchArea, minInd)

            if minInd>1 && minInd<obj.N_scan_w
                centerVal      = predVec(minInd);
                searchArea     = searchArea*obj.zoomFactor_theta0;
                predVec = linspace( ...
                    max(0,centerVal-searchArea), ...
                    min(1,centerVal+searchArea), ...
                    obj.N_scan_w);  
                predVec = permute(predVec,[1,3,2]);
            end
        end

    end

    methods(Static)

        %% ~~~~~~~~~~~~~~~~~~~~ Debug_cost_matrix ~~~~~~~~~~~~~~~~~~~~ %%
        function [] = Debug_cost_matrix(wInd, theta0Ind, t, y, w_pred_vec, theta0_pred_vec)
            w = w_pred_vec(wInd);
            theta0 = theta0_pred_vec(theta0Ind);
            yPred   = mod(w .* t + theta0, 1);
            sPlot.title     = 'Signal vs Time';
            sPlot.xlabel    = 'time [sec]';
            sPlot.ylabel    = 'signal [rad*\pi]';
            sPlot.legend    = {'y_{meas}',"y_{pred}"};
            signalToPlot    = horzcat(y,yPred);
            RMSE.Plot_signal(t, signalToPlot ,sPlot);
        end

    end
end