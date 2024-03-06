function [sEWA] = Estimate_wheel_angle()

close all;

%% Inputs
methodList = "MDFT";
% methodList = ["MDFT","EKF","PF","GS"];
noiseMode  = "On";
N          = 1000;
w          = 2400;
theta0     = 0.17;
bShowPlots = true;
Nstats     = 1;
seed       = 1111;
rng(seed);
bSaveResults = true;

%% When running statistics do no show plots
if Nstats>1
    bShowPlots = false;
else
    bSaveResults = false;
end

%% Noise Mode
switch noiseMode
    case "On"
        sigma_n  = 0.03;
        pUniform = 0.05;
    case "Off"
        sigma_n  = 0;
        pUniform = 0;
    case "Only"
        sigma_n  = 0.03;
        pUniform = 0;
end

%% Construct object of: Estimate Angular Velocity
oEWA = cEstimateWheelAngle(w, theta0, sigma_n, pUniform, N, bShowPlots);

%% Statistics Loop
for n = 1:Nstats

    %% Display Statistics
    if ~bShowPlots
        disp("Stat: "+n);
    end

    %% Generate Input Samples
    [y, t] = oEWA.GenerateInputSamples();

    %% Generate Init Guess W, theta0
    oEWA.GenerateInitGuess(y, t)

    %% Estimation Methods for w & theta0
    for method = methodList
        switch method
            case "MDFT"
                %% Modified DFT
                oMDFT = cModifiedDFT(oEWA);
                sEstOut = oMDFT.Apply(y, t);

            case "EKF"
                %% Extended Kalman Filter
                oEKF = cExtendedKalmanFilter(oEWA);
                sEstOut = oEKF.Apply(y, t);

            case "PF"
                %% Particle Filter
                oPF = cParticleFilter(oEWA);
                sEstOut = oPF.Apply(y, t);

            case "GS"
                %% Grid Search
                oGS = cGridSearch(oEWA);
                sEstOut = oGS.Apply(y, t);

        end

        %% Aggregate Results
        sEst.(method).w(n,1) = sEstOut.w_est;
        sEst.(method).theta0(n,1) = sEstOut.theta0_est;

        sRMSE.(method).w(n,1) = sEstOut.w_RMSE;
        sRMSE.(method).theta0(n,1) = sEstOut.theta0_RMSE;

    end
end

%% Average RMSE
for method = methodList
    sRMSE_Average.(method).w = 10*log10(mean(sRMSE.(method).w));
    sRMSE_Average.(method).theta0 = 10*log10(mean(sRMSE.(method).theta0));
end

%% Output Struct
sEWA.sEst          = sEst;
sEWA.sRMSE         = sRMSE;
sEWA.sRMSE_Average = sRMSE_Average;

%% Save Results
if bSaveResults
    filenameToSave = "C:\GitHub\Sandboxes\Interviews\Kayhut\Performance_Results.mat";
    disp("Saving Results to: "+filenameToSave+" ...");
    save(filenameToSave,'sEWA','-v7.3');
    disp("Saved Results to: "+filenameToSave+" !!!");
end

end
