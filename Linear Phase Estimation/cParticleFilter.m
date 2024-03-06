classdef cParticleFilter < handle
    properties
        
        %% Objects
        oEWA                    % object of cEstimateWheelAngle

        %% PF Properties
        J                       % Number of particles
        wInitStdScale           % Scale of initial Std of w regards to w_init
        w_init_area             % Initial search area for w
        w_init_std              % Initial std for w
        theta0_init_std         % Initial std for theta0
        w_pred_noise_std        % Prediction noise std for w
        theta0_pred_noise_std   % Prediction noise std for theta0
        meas_noise_std          % Measurements noise std
        alphaIIR_w              % IIR factor for adding noise to the particles for w
        alphaIIR_theta0         % IIR factor for adding noise to the particles for theta0
        min_meas_noise_sigma    % minimum value for the measurement noise sigma
        numSamplesWindow        % Number of samples from the past iterations to process in each iter
        estMethod               % Estimation method: Best weight or average all weights
        numMaxIterations        % Number of maximum iterations to converge
        DminTh                  % Minimum Distance for stop condition

    end

    methods

        %% ~~~~~~~~~~~~~~~~~~~~~~ Constructor ~~~~~~~~~~~~~~~~~~~~~~ %%
        function obj = cParticleFilter(oEWA)
            
            %% Objects
            obj.oEWA                    = oEWA;

            %% Particle Filter
            obj.J                       = 200;           % Number of particles
            obj.wInitStdScale           = 0.2;           % Scale of initial Std of w regards to w_init
            obj.w_init_area             = 20;            % Initial std for w
            obj.w_init_std              = obj.oEWA.w_init*obj.wInitStdScale; %max([obj.oEWA.w_init_std/2, obj.oEWA.w_init*obj.wInitStdScale,obj.w_init_area]);
            obj.theta0_init_std         = 0.2;           % Initial std for theta0
            obj.w_pred_noise_std        = 0.3;           % Prediction noise std for w
            obj.theta0_pred_noise_std   = 0.05;          % Prediction noise std for theta0
            obj.meas_noise_std          = oEWA.sigma_n;  % Measurements noise std
            obj.alphaIIR_w              = 0.8;           % IIR factor for adding noise to the particles w
            obj.alphaIIR_theta0         = 0.9;           % IIR factor for adding noise to the particles theta0
            obj.min_meas_noise_sigma    = 1e-3;          % minimum value for the measurement noise sigma
            obj.numSamplesWindow        = 100;           % M: Number of samples from the past iterations to process in each iter
            obj.estMethod               = "Best";        % Estimation method: Best weight or average all weights
            obj.numMaxIterations        = 17;            % Number of maximum iterations to converge
            obj.DminTh                  = 1e-3;          % Minimum Distance for stop condition

        end
 
        %% ~~~~~~~~~~~~~~~~~~~~ Apply ~~~~~~~~~~~~~~~~~~~~ %%
        function [sPerfOut] = Apply(obj, y, t)
            % Particle Filter to estimate w and theta0

            %% Initialize particles
            new_particles = obj.Init_particles();

            %% Init Estimation Vectors
            wEstVec         = zeros(obj.numMaxIterations,1);
            theta0EstVec    = zeros(obj.numMaxIterations,1);
            DminVec         = zeros(obj.numMaxIterations,1);

            %% Init weights
            weights = ones(obj.J,1)/obj.J;

            %% Iteration Loop --> O(N*J*M)
            for iter = 1:obj.numMaxIterations
                %% Show iterations
                if obj.oEWA.bShowPlots
                    disp("Iter = "+iter);
                end

                %% Add gaussian noise to particles --> O(J)
                particles = obj.Add_gaussian_noise_to_particles(new_particles);

                %% Prediction: based on the model over hypothesizes (particles)
                yPred = obj.Prediction_based_on_particles(particles, t);

                %% Generate weights --> O(J*M)
                [weights, bValidWeights, Dmin] = obj.Generate_weights(yPred, y, weights);

                %% Debug: Plot 2 best hypothesizes
                bPlot2BestHypo = false;
                if bPlot2BestHypo
                    obj.Plot_2_best_hypothesizes(y, yPred, weights, t);
                end

                %% Resample particles: Low Variance Resampling --> O(J)
                if bValidWeights
                    new_particles = obj.Low_variance_resampling(particles, weights);
                else
                    new_particles = particles;
                end

                %% Estimate_w_theta0 --> O(J)
                [sPerfOut, sPlot]    = obj.Estimate_w_theta0(particles, weights, t);

                %% Update Estimation Vectors
                wEstVec(iter)       = sPerfOut.w_est;
                theta0EstVec(iter)  = sPerfOut.theta0_est;
                DminVec(iter)  = Dmin;
                
                %% Stop Condition
                if Dmin<obj.DminTh
                    wEstVec         = wEstVec(1:iter);
                    theta0EstVec    = theta0EstVec(1:iter);
                    DminVec         = DminVec(1:iter);
                    break;
                end
            end

            %% Show results
            if obj.oEWA.bShowPlots
                %% Show_estimation_errors_vs_iterations
                sObjectives.wEstVec         = wEstVec;
                sObjectives.theta0EstVec    = theta0EstVec;
                sObjectives.estErr          = DminVec;
                obj.oEWA.Show_estimation_errors_vs_iterations(1:iter, sObjectives);

                %% Plot Final Estimation
                obj.Plot_final_estimation(sPlot, t, y)
            end
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Init_particles ~~~~~~~~~~~~~~~~~~~~ %%
        function [particles] = Init_particles(obj)
            particles = zeros(obj.J, 2);
            particles(:, 1) = obj.oEWA.w_init + obj.w_init_std * randn(obj.J, 1); % Initialize w around w_init
            particles(:, 2) = mod(obj.oEWA.theta0_init + obj.theta0_init_std * randn(obj.J, 1),1); % Initialize w around w_init
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Add_gaussian_noise_to_particles ~~~~~~~~~~~~~~~~~~~~ %%
        function [particles] = Add_gaussian_noise_to_particles(obj, particles)

            %% Add gaussian noise to the particles
            newParticlesW = particles(:, 1) + obj.w_pred_noise_std * randn(obj.J, 1);
            newParticlesTheta0 = particles(:, 2) + obj.theta0_pred_noise_std * randn(obj.J, 1);

            %% IIR Filtering
            particles(:, 1) = obj.alphaIIR_w * particles(:, 1) + (1-obj.alphaIIR_w) * newParticlesW;
            particles(:, 2) = mod(obj.alphaIIR_theta0 * particles(:, 2) + (1-obj.alphaIIR_theta0) * newParticlesTheta0,1);

        end

        %% ~~~~~~~~~~~~~~~~~~~~ Add_gaussian_noise_to_particles ~~~~~~~~~~~~~~~~~~~~ %%
        function [yPred] = Prediction_based_on_particles(obj, particles, tVec)
            yPred = mod(particles(:, 1).*tVec.'+particles(:, 2),1);
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Low_variance_resampling ~~~~~~~~~~~~~~~~~~~~ %%
        function [newWeights, bValidWeights, Dmin] = Generate_weights(obj, yPred, yMeas, weights)

            %% Distance from measurement
            d = yPred-yMeas.';

            %% Set lower limit for the measurement noise sigma
            meas_noise_sigma = max(obj.min_meas_noise_sigma,obj.meas_noise_std);

            %% Mahalanobis distance-based gaussian metric
            D_mahalanobis   = abs(d/meas_noise_sigma).^2;
            D               = obj.oEWA.Non_outlier_average(D_mahalanobis);
            newWeights      = 1/(2*pi*meas_noise_sigma^2).*exp(-0.5*D);

            %% Add regularization constant for stability
            newWeights = newWeights+eps;

            %% Normalize weights
            newWeights = newWeights./sum(newWeights);
            
            %% Check weight validity
            bValidWeights = max(newWeights)>1/obj.J+eps;
            if ~bValidWeights
                newWeights = weights;
            end

            %% Minimum Distance for stop condition
            Dmin = min(D);

        end

        %% ~~~~~~~~~~~~~~~~~~~~ Low_variance_resampling ~~~~~~~~~~~~~~~~~~~~ %%
        function [new_particles] = Low_variance_resampling(obj, particles, weights)
            step = 1/obj.J;
            r = step*rand();
            c = weights(1);
            i = 1;
            new_particles = zeros(obj.J,2);
            for j=1:obj.J
                U = r+(j-1)*step;
                while U>c && i<obj.J-1
                    i = i+1;
                    c = c+weights(i);
                end
                new_particles(j,:) = particles(i,:);
            end
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Estimate_w_theta0 ~~~~~~~~~~~~~~~~~~~~ %%
        function [sPerfOut, sPlot] = Estimate_w_theta0(obj, particles, weights, t)

            [sortedWeights,sortedInds] = sort(weights,'descend');

            switch obj.estMethod
                case "Best"
                    %% Take the highest weight for best estimate
                    w_est           = particles(sortedInds(1),1);
                    theta0_est      = particles(sortedInds(1),2);
                case "Average"
                    %% Average all weights for each of the chosen particles has even opportunity
                    w_est           = mean(particles(:,1));
                    theta0_est      = mean(particles(:,2));
                    
            end

            %% Prediction
            yPredBest       = mod(w_est.*t.'+theta0_est,1);

            %%  RMSE
            sRMSE = obj.oEWA.Calc_RMSE_Error(w_est, theta0_est);

            RMSE_w_title      = "w: True: " + obj.oEWA.w + " , Est: " + w_est + " , RMSE: = "+sRMSE.w_dB+" [dB]";
            RMSE_theta0_title = "\theta_0: True: " + obj.oEWA.theta0 + " , Est: " + theta0_est + " , RMSE:  = "+sRMSE.theta_dB+" [dB]";

            %% Out Struct
            sPerfOut.w_est          = w_est;
            sPerfOut.theta0_est     = theta0_est;
            sPerfOut.w_RMSE         = sRMSE.w;
            sPerfOut.theta0_RMSE    = sRMSE.theta;
            sPerfOut.w_RMSE_dB      = sRMSE.w_dB;
            sPerfOut.theta0_RMSE_dB = sRMSE.theta_dB;

            sPlot.maxSortedInd = sortedInds(1);
            sPlot.yPredBest  = yPredBest;
            sPlot.RMSE_w_title = RMSE_w_title;
            sPlot.RMSE_theta0_title = RMSE_theta0_title;

        end

        %% ~~~~~~~~~~~~~~~~~~~~ Plot Final Estimation ~~~~~~~~~~~~~~~~~~~~ %%
        function [] = Plot_final_estimation(obj, sPlot, t, y)
            sPlot.title     = {
                "Particle Filter",...
                "Signal vs Time",...
                sPlot.RMSE_w_title,...
                sPlot.RMSE_theta0_title,...
                "w: "+obj.oEWA.w,...
                "\theta_0: "+obj.oEWA.theta0,...
                "\sigma_n: "+obj.oEWA.sigma_n,...
                "pUniform: "+obj.oEWA.pUniform,...
                "N Samples: "+obj.oEWA.N,...
                "N Particles: "+obj.J,...
                };
            sPlot.xlabel    = 'time [sec]';
            sPlot.ylabel    = 'signal [rad/\pi]';
            sPlot.legend    = {'y_{meas}',"y_{pred,1}: "+sPlot.maxSortedInd};
            signalToPlot    = horzcat(y,sPlot.yPredBest.');
            obj.oEWA.Plot_signal(t, signalToPlot ,sPlot);
        end

        

    end

    methods(Static)

        %% ~~~~~~~~~~~~~~~~~~~~ %% Debug: Plot 2 best hypothesizes ~~~~~~~~~~~~~~~~~~~~ %%
        function [] = Plot_2_best_hypothesizes(yMeas, yPred, weights, tVec)
            [~,sortedInds]  = sort(weights,'descend');
            sPlot.title     = 'signal vs time';
            sPlot.xlabel    = 'time [sec]';
            sPlot.ylabel    = 'signal [rad/\pi]';
            sPlot.legend    = {'y_{meas}','y_{pred}^(1)','y_{pred}^(2)'};
            signal          = horzcat(yMeas,yPred(sortedInds(1),:).',yPred(sortedInds(2),:).');
            cEstimateWheelAngle.Plot_signal(tVec, signal, sPlot);
        end

    end
end