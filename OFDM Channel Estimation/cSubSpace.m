classdef cSubSpace < handle
    properties

        %% cSubSpace Properties
        sParams
        inputSpacing
        peakTh

    end

    %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
    %% ~~~~~~~~~~~~~~~~~~~~~~ Methods ~~~~~~~~~~~~~~~~~~~~~~ %%
    %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
    methods

        %% ~~~~~~~~~~~~~~~~~~~~~~ Constructor ~~~~~~~~~~~~~~~~~~~~~~ %%
        function obj = cSubSpace(sParams)

            %% Update Defulats
            sParams = cSubSpace.Update_defualts(sParams);
            obj = cStruct.Update_param_struct_into_object(obj, sParams);
            obj.sParams.inputSpacing = mean(diff(obj.sParams.inputGrid));

        end

        %% ~~~~~~~~~~~~~~~~~~~~ Apply ~~~~~~~~~~~~~~~~~~~~ %%
        function [chosenTaps] = Apply(obj, X)

            %% Preprocessing Step: Compute the spatially smoothed covariance matrix
            % Adjust the size of subblocks for spatial smoothing (L) as needed
            Rxx_fb = cMath.Spatial_Smoothing(X, obj.sParams.blockSize);

            %% Real_covariance_matrix_for_Roots
            Rxx_fb = obj.Real_covariance_matrix_for_Roots(Rxx_fb);

            %% Eigen Decomposition
            [E, D] = cMath.Eigen_Decomposition(Rxx_fb, obj.sParams.eigMode);

            %% Estimate_number_of_taps
            p  = size(D,1);
            N  = size(X,2); % The inner dimension of the auto-correlation
            obj.Estimate_number_of_taps(D, p, N)

            %% Noise Whitening
            sStat = cMath.Signal_And_Noise_Decomposition(E, D, obj.sParams);

            %% Display Number of Paths
            switch obj.sParams.method
                case "MUSIC_Inverse"
                    chosenTaps = obj.Estimate_Taps_By_Inverse(sStat.En);
                case "MUSIC_Roots"
                    chosenTaps = obj.Estimate_Taps_By_Roots(sStat.En);
                case "ESPRIT"
                    chosenTaps = obj.Estimate_Taps_By_ESPRIT(sStat.Es);
                case "OMP"
                    chosenTaps = obj.Estimate_Taps_By_OMP(X, sStat.Rss);
                case "MVDR"
                    chosenTaps = obj.Estimate_Taps_By_MVDR(X, sStat.Rss);
            end

        end
    

        %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
        %% ~~~~~~~~~~~~~~~~~~~~~~~~ Common Processses ~~~~~~~~~~~~~~~~~~~~~~~ %%
        %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%

        %% ~~~~~~~~~~~~~~~~~~~~ Estimate_number_of_taps ~~~~~~~~~~~~~~~~~~~~ %%
        function [] = Estimate_number_of_taps(obj, D, p, N)
            if obj.sParams.numPaths==0 % 0 means automatic find the optimal number of paths using MDL
                
                %% MDL: Minimum_Description_Length
                obj.sParams.numPaths = cMath.Minimum_Description_Length(D, p, N, obj.sParams);

            end
            obj.sParams.numPaths = min(obj.sParams.numPaths,obj.sParams.blockSize-1);
            if obj.sParams.bShowPlots
                disp("Number of Paths: "+obj.sParams.numPaths);
            end
        end
 
        %% ~~~~~~~~~~~~~~~~~~~~ Find_spectrum_peak_indices ~~~~~~~~~~~~~~~~~~~~ %%
        function [peakInds, peakVals] = Find_spectrum_peak_indices(obj, PSpec)
            PSpecForPeakDetection = vertcat(0, abs(PSpec)); % Adjust for zero indexing
            [peakVals, peakInds] = findpeaks(PSpecForPeakDetection, ...
                'SortStr', 'descend', 'NPeaks', obj.sParams.numPaths);
            peakInds       = peakInds - 1; % Adjusting peak locations
            thPercent       = 100*(1-2*obj.sParams.numPaths/numel(PSpecForPeakDetection));
            thPercent       = min(thPercent,95);
            percentVal      = cMath.Get_cdf_percent(PSpecForPeakDetection,thPercent);
            medVal          = median(PSpecForPeakDetection);
            obj.peakTh     = percentVal+medVal;%medVal+noiseVarPostPG*obj.sParams.numStdsTh;
            validInds      = peakVals>obj.peakTh;
            if ~any(validInds)
                obj.peakTh = 0;
                validInds  = peakVals>obj.peakTh;
            end
            peakMat        = horzcat(peakInds, peakVals);
            sortedPeakMat  = sortrows(peakMat(validInds,:),1);
            peakInds       = sortedPeakMat(:,1);
            peakVals       = sortedPeakMat(:,2);
            obj.sParams.numPaths = numel(peakInds);
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Get_tap_delays ~~~~~~~~~~~~~~~~~~~~ %%
        function [chosenTaps] = Get_tap_delays(obj, peakInds)
            chosenTaps = obj.sParams.outputGrid(peakInds);
        end


        %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
        %% ~~~~~~~~~~~~~~~~~~~~~~~~~~ MUSIC Inverse ~~~~~~~~~~~~~~~~~~~~~~~~~ %%
        %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%

        %% ~~~~~~~~~~~~~~~~~~~~ Estimate_Taps_By_Inverse ~~~~~~~~~~~~~~~~~~~~ %%
        function [chosenTaps] = Estimate_Taps_By_Inverse(obj, En)

            %% Step 1: SubSpace algorithm for estimating tap delays
            PSpecVec = obj.Calc_P_MUSIC(En);

            %% Step 2: Identify the peaks in the MUSIC spectrum to estimate the tap delays
            peakInds = obj.Find_spectrum_peak_indices(PSpecVec);

            %% Step 3: Tap Delays
            chosenTaps = obj.Get_tap_delays(peakInds);

            %% Plot the MUSIC spectrum
            if obj.sParams.bShowPlots
                obj.Plot_P_spectrum(PSpecVec, chosenTaps);
            end

        end

        %% ~~~~~~~~~~~~~~~~~~~~ Calc_P_MUSIC ~~~~~~~~~~~~~~~~~~~~ %%
        function [PSpecVec] = Calc_P_MUSIC(obj, En)
            PSpecVec = zeros(size(obj.sParams.outputGrid));
            for i = 1:numel(obj.sParams.outputGrid)
                outputVal = obj.sParams.outputGrid(i);
                % Adjusted steering vector construction
                inputVec = obj.sParams.inputGrid(1:size(En,1));
                steeringVec = exp(-1i*2*pi*inputVec*outputVal).'; % Steering vector
                % Adjusted MUSIC spectrum calculation
                PSpecVec(i) = 1 / (steeringVec' * (En * En') * steeringVec); 
            end
        end

        %% ~~~~~~~~~~~~~~~~~~~~ %% Plot_P_spectrum ~~~~~~~~~~~~~~~~~~~~ %%
        function [] = Plot_P_spectrum(obj, PSpecVec, chosenTaps)
            
            PSpecPowerdB    = 10*log10(abs(PSpecVec));
            PSpecPowerRoots = interp1(obj.sParams.outputGrid,PSpecPowerdB,chosenTaps);
            peakTh_dB       = 10*log10(obj.peakTh);
            
            figure;
            plot(obj.sParams.tapFactorToShow*obj.sParams.outputGrid, PSpecPowerdB, '-o');
            xlabel(obj.sParams.tapType+" Taps ["+obj.sParams.tapUnitsToShow+"]");
            ylabel('Sub-Space Spectrum [dB]');
            
            [~,maxInd] = max(PSpecPowerRoots);
            titleStr = vertcat(obj.sParams.baseTitleStr,...
                "Number of Paths: "+obj.sParams.numPaths,...
                obj.sParams.tapType+" Tap: "+obj.sParams.tapFactorToShow*chosenTaps(:)+" ["+obj.sParams.tapUnitsToShow+"] , Power: "+PSpecPowerRoots+" [dB]",...
                "Strongest Tap",...
                obj.sParams.tapType+" Tap: "+obj.sParams.tapFactorToShow*chosenTaps(maxInd)+" ["+obj.sParams.tapUnitsToShow+"] , Power: "+PSpecPowerRoots(maxInd)+" [dB]",...
                "Largest Tap",...
                obj.sParams.tapType+" Tap: "+obj.sParams.tapFactorToShow*chosenTaps(end)+" ["+obj.sParams.tapUnitsToShow+"] , Power: "+PSpecPowerRoots(end)+" [dB]");
            title(titleStr);
            grid minor;
            hold on;
            plot(obj.sParams.tapFactorToShow*chosenTaps, PSpecPowerRoots, 'rx', 'MarkerSize', 8, 'LineWidth',2);
            if ~isempty(peakTh_dB)
                yline(peakTh_dB,'g--','LineWidth',2);
                legendStr = {'Sub-Space Spectrum', 'Estimated Taps','Minimum Tap Power'};
            else
                legendStr = {'Sub-Space Spectrum', 'Estimated Taps'};
            end

            %% ClickLegendCb
            cPlot.Click_Legend(legendStr);
            plotbrowser('on');
        end


        %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
        %% ~~~~~~~~~~~~~~~~~~~~~~~~~~ MUSIC Roots ~~~~~~~~~~~~~~~~~~~~~~~~~ %%
        %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%

        %% ~~~~~~~~~~~~~~~~~~~~ Estimate_Taps_By_Roots ~~~~~~~~~~~~~~~~~~~~ %%
        function [chosenTaps] = Estimate_Taps_By_Roots(obj, Q_N)

            %% Step 2: Define the FB root-MUSIC polynomial and find its roots
            f_mu_fb = obj.Generate_steering_vector_function(Q_N);

            %% Step 3: Find Roots Using Newtons Method
            initialGuesses = exp(-1i*2*pi*obj.sParams.inputSpacing*obj.sParams.outputGrid); % Adjust range and count as necessary
            roots          = obj.Find_Roots_Using_Newtons_Method(f_mu_fb, initialGuesses);

            %% Step 4: Convert roots to angles and tap delays
            chosenTaps = sort(-angle(roots)/(2*pi*obj.sParams.inputSpacing));

            %% Debug: Plot the MUSIC spectrum
            if obj.sParams.bShowPlots
                %% Generate PSpec Vector
                PSpecVec = zeros(numel(initialGuesses),1);
                for i = 1:numel(initialGuesses)
                    PSpecVec(i) = f_mu_fb(initialGuesses(i));
                end
                %% Plot_P_spectrum
                obj.Plot_P_spectrum(PSpecVec, chosenTaps);
            end

        end

        %% ~~~~~~~~~~~~~~~~~~~~ Real_covariance_matrix_for_Roots ~~~~~~~~~~~~~~~~~~~~ %%
        function [Rxx_re] = Real_covariance_matrix_for_Roots(obj, Rxx_fb)
            if obj.sParams.method=="MUSIC_Roots"
                %% Step 0: Number of elements in the array
                N = size(Rxx_fb, 1);

                %% Step 1: Construct U for transforming Rxx_fb to Rxx_re
                obj.sParams.U = obj.Construct_U(N);

                %% Step 2: Transform to real matrix
                Rxx_re = real(obj.sParams.U' * Rxx_fb * obj.sParams.U); % Transformation to real-valued matrix
            else
                Rxx_re = Rxx_fb;
            end
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Generate_steering_vector_function ~~~~~~~~~~~~~~~~~~~~ %%
        function [f_mu_fb] = Generate_steering_vector_function(obj, Q_N)
        % The steering vector function
            QQH     = Q_N * Q_N';
            kVec    = (0:(obj.sParams.blockSize-1)).';
            s_z     = @(z) obj.sParams.U'*z.^kVec; % Define steering vector as a function of z
            f_mu_fb = @(z) s_z(1./z)' * QQH * s_z(z); % Define the FB root-MUSIC polynomial
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Find_Roots_Using_Newtons_Method ~~~~~~~~~~~~~~~~~~~~ %%
        function rootsNM = Find_Roots_Using_Newtons_Method(obj, f_mu_fb, initialGuesses)
            % f_mu_fb: the projection function of the steering hypo vectors and the null space eigenvectors
            % initialGuesses: Vector of initial guesses for the roots

            %% Initialize roots array
            rootsNM = zeros(numel(initialGuesses),1); 
            %% Loop over initial guesses
            for i = 1:length(initialGuesses)

                %% Current guess
                z = initialGuesses(i); 

                %% Iterations until convergence
                for iter = 1:obj.sParams.maxIterNewtonMethod
                    %% Numerical differentiation for derivative
                    delta = 1e-6; % Small change
                    f_prime = (f_mu_fb(z + delta) - f_mu_fb(z)) / delta; % Numerical derivative

                    %% Newton's update
                    z_new = z - f_mu_fb(z) / f_prime;

                    %% Check for convergence
                    if abs(z_new - z) < obj.sParams.toleranceNewtonMethod
                        rootsNM(i) = z_new;
                        break; % Converged
                    end

                    %% Prepare for next iteration
                    z = z_new;

                end
            end

            %% Filter roots to keep only those that make physical sense, e.g., inside the unit circle
            rootsNM = rootsNM(abs(rootsNM) <= 1 & rootsNM ~= 0); % Remove non-converged and out-of-bound roots

            %% Unique roots with tolerance
            rootsAbs                = abs(rootsNM);
            rootsPhase              = angle(rootsNM);
            A                       = horzcat(rootsAbs,rootsPhase);
            [uniqueVals, uniqueInd] = obj.Unique_Rows_Tol(A, obj.sParams.toleranceNewtonMethod);
            rootsNM                   = rootsNM(uniqueInd);

            %% Update number of paths
            obj.sParams.numPaths            = min(numel(rootsNM),obj.sParams.numPaths);

             %% Generate PSpec Vector
             PSpecVec = zeros(numel(rootsNM),1);
             for i = 1:numel(rootsNM)
                 PSpecVec(i) = abs(f_mu_fb(rootsNM(i)));
             end

            %% Choose the numPaths closest roots to unit circle
            distFromUnitCircle = 1-abs(rootsNM);
            rootsDistMat       = horzcat((1:numel(rootsNM)).',rootsNM, PSpecVec, distFromUnitCircle);
            sortedRoots        = sortrows(rootsDistMat,[4 3]);
            rootsNM              = sortedRoots(1:obj.sParams.numPaths,2);

        end


        %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
        %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ESPRIT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
        %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%

        %% ~~~~~~~~~~~~~~~~~~~~ Estimate_Taps_By_ESPRIT ~~~~~~~~~~~~~~~~~~~~ %%
        function [chosenTaps] = Estimate_Taps_By_ESPRIT(obj, Es)
            % Rxx_fb is the estimated smoothed foreward backward covariance matrix of the 
            % channel frequency response X

            %% Step 1: ESPRIT algorithm: Least Squares Estimation of the Phase Shift
            S1  = Es(1:end-1, :);
            S2  = Es(2:end, :);
            Psi = S1 \ S2;

            %% Step 2: Eigenvalues of Psi give the phase shifts
            switch obj.sParams.eigMode
                case "Matlab"
                    phi = sort(-angle(eig(Psi)));
                case "QR"
                    [~,D_Psi] = cMath.Eig_QR(Psi);
                    phi = sort(-angle(diag(D_Psi)));
            end

            %% Step 3: Convert phase shifts to time delays
            chosenTaps = phi / (2 * pi * obj.sParams.inputSpacing);

        end


        %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
        %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ OMP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
        %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%

        %% ~~~~~~~~~~~~~~~~~~~~ Estimate_Taps_By_OMP ~~~~~~~~~~~~~~~~~~~~ %%
        function [chosenTaps] = Estimate_Taps_By_OMP(obj, X, Rss)
            
            %% Covariance Martix Size
            Nfb            = size(Rss,1);

            inputVec       = obj.sParams.inputGrid(1:Nfb);
            F_OverCompDict = exp(-1i*2*pi*inputVec.*obj.sParams.outputGrid).';
            W_t            = zeros(size(Rss,1),obj.sParams.numPaths);
            
            XGal           = X(1:Nfb,:);
            r_t            = XGal;
            LTotal         = zeros(size(F_OverCompDict,2),obj.sParams.numPaths);
            maxLVec        = zeros(1,obj.sParams.numPaths);
            J_t            = zeros(1,obj.sParams.numPaths);

            %% Generate_MVDR_W: MVDR Beamformer Matrix
            switch obj.sParams.modeOMP
                case "MVDR"
                    W = obj.Generate_MVDR_W(F_OverCompDict, Rss);
                case "FFT"
                    W = F_OverCompDict;
            end
            for t = 1:obj.sParams.numPaths

                %% MVDR_Projection
                PSpec = obj.MVDR_Projection(W, r_t);

                %% Projection of atoms
                L = abs(PSpec);
                
                %% argmax over the projection of atoms
                [maxL,j_t] = max(L);
                
                %% Stop Condition
                if maxL<obj.sParams.thOMP
                    LTotal               = LTotal(:,1:t-1) ;
                    maxLVec              = maxLVec(:,1:t-1);
                    J_t                  = J_t(:,1:t-1)    ;
                    obj.sParams.numPaths = t-1;
                    break;
                end

                %% Plot Signal
                if true && obj.sParams.bShowPlots
                    obj.Plot_OMP_Per_Iteration(j_t, L, t, maxL)
                end

                %% Aggregate atoms
                LTotal(:,t)  = L;
                maxLVec(:,t) = maxL;
                J_t(:,t)     = j_t;
                W_t(:,t)     = W(:,j_t).';

                %% Least Squares for optimal coeffiecients Fn
                Wt       = W_t(:,1:t);
                aLS      = Wt \ XGal;
                XHat_t   = Wt * aLS;

                %% Subtract Residual
                r_t      = XGal-XHat_t;
            end

            %% Chosen Taps
            chosenTaps = obj.sParams.outputGrid(J_t);

            %% Show_OMP_Spec_Image
            if obj.sParams.bShowPlots
                obj.Show_OMP_Spec_Image(LTotal, J_t, chosenTaps, maxLVec);
            end

            
        end

        %% ~~~~~~~~~~~~~~~~~~~~ %% Show_OMP_Spec_Image ~~~~~~~~~~~~~~~~~~~~ %%
        function [] = Show_OMP_Spec_Image(obj, LTotal, J_t, chosenTaps, maxLVec)

            Piters = numel(J_t);
            PSpecPowerRoots = maxLVec(:);
            [~,maxInd] = max(PSpecPowerRoots);

            %% Channel Time Struct
            figure;
            gridVec    = obj.sParams.tapFactorToShow*obj.sParams.outputGrid;
            Tdecimated = 60;
            T          = length(gridVec);
            M          = round(T/Tdecimated);
            gridDecVec    = gridVec(1:M:end);

            imagesc(LTotal(:,1:Piters)); colorbar;
            hold on;
            for k = 1:numel(J_t)
                plot(k, J_t(k),'marker','x','LineWidth',2,'MarkerSize',10,'Color','red')
            end
            set(gca, 'XTick', 1:obj.sParams.numPaths, 'XTickLabel', 1:obj.sParams.numPaths,'XTickLabelRotation',90) % 10 ticks
            set(gca, 'YTick', 1:M:numel(gridVec), 'YTickLabel', gridDecVec) % 20 ticks
                       
            ylabel(obj.sParams.tapType+" ["+obj.sParams.tapUnitsToShow+"]");
            xlabel('Iterations (# Paths)');
            title({"OMP Spectrum Projection",obj.sParams.tapType+" Estimation"});
            titleStr = vertcat(obj.sParams.baseTitleStr,...
                "Number of Paths: "+obj.sParams.numPaths,...
                obj.sParams.tapType+" Tap: "+obj.sParams.tapFactorToShow*chosenTaps(:)+" ["+obj.sParams.tapUnitsToShow+"] , Power: "+PSpecPowerRoots,...
                "Strongest Tap",...
                obj.sParams.tapType+" Tap: "+obj.sParams.tapFactorToShow*chosenTaps(maxInd)+" ["+obj.sParams.tapUnitsToShow+"] , Power: "+PSpecPowerRoots(maxInd),...
                "Largest Tap",...
                obj.sParams.tapType+" Tap: "+obj.sParams.tapFactorToShow*chosenTaps(end)+" ["+obj.sParams.tapUnitsToShow+"] , Power: "+PSpecPowerRoots(end));
            title(vertcat("OMP Spectrum Projection",titleStr));
            plotbrowser on;
        end

        %% ~~~~~~~~~~~~~~~~~~~~ %% Plot_OMP_Per_Iteration ~~~~~~~~~~~~~~~~~~~~ %%
        function [] = Plot_OMP_Per_Iteration(obj, j_t, L, t, maxL)
            maxTap = obj.sParams.outputGrid(j_t)*obj.sParams.tapFactorToShow;
            sPlot.title     = {...
                "Projection of Atoms over Signal",...
                "Iteration: "+t,...
                "Max "+obj.sParams.tapType+": "+maxTap+" ["+obj.sParams.tapUnitsToShow+"]",...
                };
            sPlot.xlabel    = obj.sParams.tapType+" ["+obj.sParams.tapUnitsToShow+"]";
            sPlot.ylabel    = 'Projection of Atoms over Signal';
            signalToPlot    = L;
            cPlot.Plot_signal(obj.sParams.outputGrid*obj.sParams.tapFactorToShow, signalToPlot ,sPlot);
            hold on;
            plot(maxTap,maxL,'gx','MarkerSize',17,'LineWidth',2);
        end


        %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
        %% ~~~~~~~~~~~~~~~~~~~~~~~ Modified MVDR ~~~~~~~~~~~~~~~~~~~~~~~~ %%
        %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%

        %% ~~~~~~~~~~~~~~~~~~~~ Estimate_Taps_By_MVDR ~~~~~~~~~~~~~~~~~~~ %%
        function [chosenTaps] = Estimate_Taps_By_MVDR(obj, X, Rss)

            %% Covariance Martix Size
            Nfb = size(Rss,1);

            %% Step 1: Steering Matrix
            inputVec    = obj.sParams.inputGrid(1:Nfb);
            steeringMat = exp(-1i*2*pi*inputVec.*obj.sParams.outputGrid).';

            %% Generate_MVDR_W: MVDR Beamformer Matrix
            W = obj.Generate_MVDR_W(steeringMat, Rss);

            %% Step 2: MVDR beamforming for each potential delay tap
            PSpecVec = obj.MVDR_Projection(W, X);

            %% Step 3: Find_spectrum_peak_indices
            [peakInds, peakVals] = obj.Find_spectrum_peak_indices(PSpecVec);

            %% Step 4: Tap Delays
            chosenTaps = obj.Get_tap_delays(peakInds);

            %% Plot the MUSIC spectrum
            if obj.sParams.bShowPlots
                obj.Plot_P_spectrum(PSpecVec, chosenTaps);
            end

        end


    end

    %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
    %% ~~~~~~~~~~~~~~~~~~ Static Methods ~~~~~~~~~~~~~~~~~~~ %%
    %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
    methods(Static)

        %% ~~~~~~~~~~~~~~~~~~~~ Apply ~~~~~~~~~~~~~~~~~~~~ %%
        function [sParams] = Update_defualts(sParams)
            sParams = cStruct.Set_default_value(sParams,'eigMode',"Matlab");
            sParams = cStruct.Set_default_value(sParams,'bWhitenNoise',false);
            sParams = cStruct.Set_default_value(sParams,'method',"Inverse");
            sParams = cStruct.Set_default_value(sParams,'numPaths',0);
            sParams = cStruct.Set_default_value(sParams,'lambdaFloorLevelMDL',1e-4);
            sParams = cStruct.Set_default_value(sParams,'modeMDL',"MDL");
            sParams = cStruct.Set_default_value(sParams,'maxIterNewtonMethod',100);
            sParams = cStruct.Set_default_value(sParams,'toleranceNewtonMethod',1e-5);
            sParams = cStruct.Set_default_value(sParams,'numStdsTh',0.1);
            sParams = cStruct.Set_default_value(sParams,'thOMP',1e-6); 
            sParams = cStruct.Set_default_value(sParams,'modeOMP',"MVDR"); 
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Construct_U ~~~~~~~~~~~~~~~~~~~~ %%
        function [U] = Construct_U(N)
            % Determine the size of Rxx_fb to decide on the U matrix
            % N: Number of elements in the array
            % Construct the unitary transformation matrix U
            if mod(N, 2) == 0
                % For even N
                I = eye(N/2);
                J = fliplr(eye(N/2));
                U = (1/sqrt(2)) * [I, 1i*I; J, -1i*J];
            else
                % For odd N, adjust the construction accordingly
                M = (N-1)/2;
                I = eye(M);
                J = fliplr(eye(M));
                upper = [I, zeros(M,1), 1i*I];
                middle = [zeros(1,M),sqrt(2),zeros(1,M)];
                lower = [J, zeros(M,1), -1i*J];
                U = (1/sqrt(2)) * [upper; middle; lower];
            end

        end

        %% ~~~~~~~~~~~~~~~~~~~~ MVDR_Projection ~~~~~~~~~~~~~~~~~~~~ %%
        function [PSpecVec] = MVDR_Projection(W, X)

            %% Covariance Martix Size
            Nfb = size(W,1);
            PSpecVec = mean(W' * X(1:Nfb,:),2); % Apply beamforming
            
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Generate_MVDR_W ~~~~~~~~~~~~~~~~~~~~ %%
        function [W] = Generate_MVDR_W(steeringMat, Rss)
            %% Covariance Martix Size
            N = size(steeringMat,2);
            W = zeros(size(steeringMat));
            for k = 1:N
                v0          = steeringMat(:,k);
                w           = cSubSpace.Generate_MVDR_w_BF(v0, Rss);
                W(:,k)      = w;
            end
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Generate_MVDR_w_BF ~~~~~~~~~~~~~~~~~~~~ %%
        function [w] = Generate_MVDR_w_BF(v0, Rss)
            alpha       = 1 / (v0' * (Rss \ v0));
            w           = alpha * (Rss \ v0); % MVDR weights
        end

        %% ~~~~~~~~~~~~~~~~~~~~ Unique Rows Tol ~~~~~~~~~~~~~~~~~~~~ %%
        function [uniqueA, uniqueInd] = Unique_Rows_Tol(A, tol, wrapAroundVal_Dim)

            %% Defaults
            % tol
            if ~exist('tol','var')
                tol = 1e-3;
            end

            % wrapAroundVal_Dim
            if ~exist('wrapAroundVal_Dim','var')
                wrapAroundVal_Dim = zeros(1,size(A,2));
            end

            %%
            maskNoWrapDim = wrapAroundVal_Dim==0;
            maskWrapDim   = wrapAroundVal_Dim > 0;
            D             = vecnorm(permute(A(:,maskNoWrapDim), [1 3 2]) - permute(A(:,maskNoWrapDim), [3 1 2]), 2, 3);
            D2            = D.^2;
            indWrapDim    = find(maskWrapDim);
            for indDim = indWrapDim
                valWrap = wrapAroundVal_Dim(indDim);
                distDim = abs(A(:,indDim) - A(:,indDim)');
                distDim = min(distDim, valWrap - distDim);
                D2      = D2+distDim.^2;
            end
            isNeighbour = D2 < tol^2;
            mask        = true(size(isNeighbour,1),1);
            uniqueInd   = zeros(0,1);
            for ii = 1:size(isNeighbour,1)
                ind2add          = find(mask,1,'first');
                thisConn         = isNeighbour(ind2add,:);
                mask(thisConn)   = false;
                uniqueInd(end+1,1) = ind2add;
                if nnz(mask) == 0
                    break;
                end
            end
            uniqueA = A(uniqueInd,:);

        end


    end

end
