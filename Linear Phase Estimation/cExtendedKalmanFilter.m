classdef cExtendedKalmanFilter < handle
    properties
        
        %% Objects
        oEWA                    % object of cEstimateWheelAngle

        %% EKF Properties
        Q                       % Process noise covariance matrix Q
        R                       % Measurement noise covariance matrix R
        P_init                  % Initial state covariance
        delta                   % Delta of w and theta0 for the Jacobian
        regularizationConst     % Regularization const for EKF

    end

    methods

        %% ~~~~~~~~~~~~~~~~~~~~~~ Constructor ~~~~~~~~~~~~~~~~~~~~~~ %%
        function obj = cExtendedKalmanFilter(oEWA)
            
            %% Objects
            obj.oEWA                    = oEWA;

            %% Extended Kalman Filter
            obj.Q                   = diag([1,1e-3]); % Process noise covariance matrix Q
            obj.R                   = 0.2 + oEWA.sigma_n^2;      % Measurement noise covariance matrix R
            obj.P_init              = diag([1e5,0.05]); % Initial state covariance
            obj.delta               = 0.05;               % Delta of w and theta0 for the Jacobian

        end

        %% ~~~~~~~~~~~~~~~~~~~~ Apply ~~~~~~~~~~~~~~~~~~~~ %%
        function [sPerfOut] = Apply(obj, y, t)
            % Particle Filter to estimate w and theta0
            % Extended Kalman Filter to estimate parameters w and theta0 from noisy signal

            %% Initialization
            x_est           = [obj.oEWA.w_init; obj.oEWA.theta0_init];  % Initial state estimate
            P               = obj.P_init;  % Initial state covariance
            z_est           = zeros(obj.oEWA.N,1);
            x_est_mat       = zeros(obj.oEWA.N,2);
            meas_err_vec    = zeros(obj.oEWA.N,1);

            %% Init Estimation Vectors
            wEstVec         = zeros(obj.oEWA.N,1);
            theta0EstVec    = zeros(obj.oEWA.N,1);

            %% Time iteration loop --> O(N)
            for k = 1:obj.oEWA.N
                %% Show iterations
                if obj.oEWA.bShowPlots
                    disp("Iteration: "+k);
                end

                %% State Prediction (identity since w and theta0 are constants)
                x_pred = x_est;
                P_pred = P + obj.Q;

                %% Measurement Update --> O(1)
                % Numerically approximate the Jacobian of the measurement function at x_pred
                H_k = obj.Compute_jacobian(x_pred, t(k));

                %% Kalman Gain --> O(1)
                K = P_pred * H_k' / (H_k * P_pred * H_k' + obj.R);

                %% Output prediction --> O(1)
                z_pred = mod(x_pred(1) * t(k) + x_pred(2), 1);  % Predicted measurement

                %% Protection from measurement error --> O(1)
                meas_err = y(k) - z_pred;
                meas_err_vec(k) = meas_err;

                if abs(meas_err)>0.2
                    K = 0; % We dont believe in the measurements
                end

                %% Update state estimate and covariance --> O(1)
                x_est  = x_pred + K * meas_err;  % State update
                P      = (eye(2) - K * H_k) * P_pred;  % Covariance update

                %% Estimated Signal
                z_est(k)       = z_pred;
                x_est_mat(k,:) = x_est.';

                %% Update Estimation Vectors
                w_est               = x_est(1);
                theta0_est          = x_est(2);
                wEstVec(k)       = w_est;
                theta0EstVec(k)  = theta0_est;

                %% Plot
                bPlotSignal    = false;
                if bPlotSignal
                    obj.Debug_plot_EKF(y, z_est, t, x_est_mat, meas_err_vec, k);
                end
            end

            %% Show results
            if obj.oEWA.bShowPlots
                %% Show_estimation_errors_vs_iterations
                sObjectives.wEstVec         = wEstVec;
                sObjectives.theta0EstVec    = theta0EstVec;
                sObjectives.estErr          = meas_err_vec;
                obj.oEWA.Show_estimation_errors_vs_iterations(t, sObjectives);
            end

            %% y Prediction
            yPredBest  = mod(w_est.*t.'+theta0_est,1);

            %%  RMSE
            sRMSE = obj.oEWA.Calc_RMSE(w_est, theta0_est, t, y, yPredBest, 'Extended Kalman Filter');

            %% Out Struct
            sPerfOut.w_est          = w_est;
            sPerfOut.theta0_est     = theta0_est;
            sPerfOut.w_RMSE         = sRMSE.w;
            sPerfOut.theta0_RMSE    = sRMSE.theta;
            sPerfOut.w_RMSE_dB      = sRMSE.w_dB;
            sPerfOut.theta0_RMSE_dB = sRMSE.theta_dB;

        end

        %% ~~~~~~~~~~~~~~~~~~~~ Debug_plot_EKF ~~~~~~~~~~~~~~~~~~~~ %%
        function [] = Debug_plot_EKF(obj, y, z_est, t, x_est_mat, meas_err_vec, k)
            %% Plot Final Estimation
            sPlot.title     = 'Signal vs Time';
            sPlot.xlabel    = 'Time [sec]';
            sPlot.ylabel    = 'Signal [rad/\pi]';
            sPlot.legend    = {'y_{meas}',"y_{pred}"};
            signalToPlot    = horzcat(y(1:k),z_est(1:k));
            obj.oEWA.Plot_signal(t(1:k), signalToPlot ,sPlot);

            sPlot.title     = 'Estimated Parameters vs time';
            sPlot.xlabel    = 'time [sec]';
            sPlot.ylabel    = 'Estimated Parameters';
            sPlot.legend    = {'w',"\theta_0","Estimation Error"};
            sPlot.yline     = [obj.oEWA.w,obj.oEWA.theta0,0];
            signalToPlot    = horzcat(x_est_mat(1:k,:),meas_err_vec(1:k));
            obj.oEWA.Plot_multiple_signals(t(1:k), signalToPlot ,sPlot);

        end

        %% ~~~~~~~~~~~~~~~~~~~~ Compute_jacobian ~~~~~~~~~~~~~~~~~~~~ %%
        function H = Compute_jacobian(obj, x, t)
            % Numerically approximate the Jacobian matrix
            h1_n       = mod((x(1) + obj.delta) * t + x(2), 1);
            h1_n_1     = mod(x(1) * t + x(2), 1);

            h2_n       = mod(x(1) * t + (x(2) + obj.delta), 1);
            h2_n_1     = mod(x(1) * t + x(2), 1);

            dh_dw      = (h1_n - h1_n_1) / obj.delta;
            dh_dtheta0 = (h2_n - h2_n_1) / obj.delta;

            H          = [dh_dw, dh_dtheta0];
        end

    end

    methods(Static)

    end
end