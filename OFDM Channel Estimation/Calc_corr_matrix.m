function [corr corrInfo] = Calc_corr_matrix(numRx,numTx,corr_name)

[eNB_params UE_params gamma nearestStandardCorrName] = parse_correlation(corr_name);
alpha = eNB_params.alpha;
beta = UE_params.beta ;
device_params(1).corr_factor = alpha;
device_params(2).corr_factor = beta;

device_params(1).setup = eNB_params.setup;
device_params(2).setup = UE_params.setup;

device_params(1).num_antennas = numTx;
device_params(2).num_antennas = numRx;

if strcmpi(device_params(2).setup, 'std') && ~isempty(strfind(corr_name, 'step_2')) && (numRx== 8 || (numRx== 4 && numTx == 4 && sysParams.frames(1).subframe(2).TrCH.PDSCH.numLayers == 4))
    device_params(2).setup = 'step_2';
end

for i_device = 1:2
    R_device{i_device} = Get_correlation_matrix(device_params(i_device));
end

%% test the number of cross poles
if length(device_params(1).setup) > length('diamond') && strcmp(device_params(1).setup(length('diamond') + 1), '_')
    device_params(1).num_of_cross_poles_constellations = str2double(device_params(1).setup(length('diamond') + 2:end));
elseif length(device_params(2).setup) > length('diamond') && strcmp(device_params(2).setup(length('diamond') + 1), '_')
    device_params(1).num_of_cross_poles_constellations = 0;
    device_params(2).num_of_cross_poles_constellations = str2double(device_params(2).setup(length('diamond') + 2:end));
elseif length(device_params(1).setup) > length('step') && strcmp(device_params(1).setup(length('step') + 1), '_')
    device_params(1).num_of_cross_poles_constellations = str2double(device_params(1).setup(length('step') + 4:end));
else
    device_params(1).num_of_cross_poles_constellations = 1;
end
if length(device_params(2).setup) > length('diamond') && strcmp(device_params(2).setup(length('diamond') + 1), '_')
    device_params(2).num_of_cross_poles_constellations = str2double(device_params(2).setup(length('diamond') + 2:end));
elseif length(device_params(2).setup) > length('step') && strcmp(device_params(2).setup(length('step') + 1), '_')
    device_params(2).num_of_cross_poles_constellations = str2double(device_params(2).setup(length('step') + 4:end));
else
    device_params(2).num_of_cross_poles_constellations = 1;
end
num_cross_poles = (device_params(1).num_of_cross_poles_constellations == 2) + (device_params(2).num_of_cross_poles_constellations == 2);
if strcmpi(device_params(2).setup, 'CP') || strcmpi(device_params(1).setup, 'CP')
    num_cross_poles = floor(device_params(1).num_antennas/2) + floor(device_params(2).num_antennas/2);
end

%% create the correlation with the num cross poles
if num_cross_poles > 0
    if num_cross_poles == 1
        R_gamma = [1 gamma ; gamma' 1];
        if (device_params(1).num_of_cross_poles_constellations == 2)
            corr_mat(:,:) = kron(kron(R_gamma,R_device{1}),R_device{2});
        else
            corr_mat(:,:) = kron(R_device{1},kron(R_gamma,R_device{2}));
        end
    elseif num_cross_poles >= 2
        if all(size(R_device{1}) == 1) || all(size(R_device{2}) == 1)
            disp('the correlation config doesnt consider alfa or beta due to single cross-pol numRx = 2')
        end
        R_gamma = [1 0 -gamma 0; ...
            0 1 0 gamma; ...
            -gamma 0 1 0; ...
            0 gamma 0 1;];

        num_Tx = device_params(1).num_antennas;
        num_Rx = device_params(2).num_antennas;
        P = zeros(num_Tx .* num_Rx);
        for ii = 1:num_Rx
            for jj = 1:(num_Tx/2)
                P((jj-1)*num_Rx + ii , 2*(jj-1)*num_Rx + ii  ) = 1;
            end
            for jj = (num_Tx/2+1):num_Tx
                P((jj-1)*num_Rx + ii , 2*(jj-num_Tx/2)*num_Rx - num_Rx + ii  ) = 1;
            end
        end
        corr_mat(:,:) = P * kron(kron(R_device{1}, R_gamma),R_device{2}) * P.';
    else
        error('the cross pole config does not support more than 2 cross poles');
    end
else
    corr_mat(:,:) = kron(R_device{1}, R_device{2});
end
%% test the eigenvalues of the correlation

% Asserting that corr_mat is Positive Semi-Definite after round-off to 4-digit precision
corr_mat_eigs = eig(round(corr_mat(:,:)*1e4)/1e4); %;eig(round(corr_mat(:,:)*1e4)/1e4) ;eig(corr_mat(:,:))
if any(corr_mat_eigs <= 0)
    error('Correlation matrix is not Positive Semi-Definite. Refer to 36.101.B.2.3.2-3 for methods to insure its PSD property after round-off to 4 digits precision');
end

corrInfo.alpha=alpha;
corrInfo.beta=beta;
corr = corr_mat;
corrInfo.corrMat=corr_mat;

end
