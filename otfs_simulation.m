% OTFS Simulation Framework
clc
clear
% Set parameters for OTFS simulation
% number of Doppler bins (time slots)
param.N=16;
% number of delay bins (subcarriers)
param.M=64;
% carrier frequency
param.fc=4e9;
% speed of light
param.c=299792458;
% subcarrier spacing
param.delta_f=15e3;
% block duration
param.T=1/param.delta_f;
% delay resolution
param.delay_resolution = 1/(param.M*param.delta_f);
% Doppler resolution
param.Doppler_resolution = 1/(param.N*param.T);
% normalized DFT matrix
Fn=dftmtx(param.N);
param.Fn=Fn/norm(Fn);
% modulation size (QAM modulation)
param.mod_size=4;
% channel model
param.model='EVA';
% maximum user equipment speed in km/h(for standard channel model)
param.max_UE_speed_kmh=100; % maximum user equipment speed in km/h
% channel parameters generation mode
param.channel_param_generation_mode='std'; % 'std': standard channel parameters, 'syn': synthetic channel parameters
% modulation method
param.otfs_modulation_method=1;% 1: Method 1, 2: Method 2, 3: Method 3
% OTFS type
param.otfs_type='OTFS'; % OTFS, RZP-OTFS, RCP-OTFS, CP-OTFS, ZP-OTFS
% response calculation method(for otfs_type='OTFS')
param.res_calc_method=2; % 1: TDL model, 2: time-domain channel matrix, 3: delay-time channel matrix, 4: delay-Doppler channel matrix
% noise power(dB)
param.SNR_dB=Inf;
% OTFS demodulation method
param.otfs_demod_method=3; % 1: Method 1, 2: Method 2, 3: Method 3
% OTFS detection method
param.otfs_detection_method='LMMSE_time'; % LMMSE_delay-Doppler, LMMSE_time
% row-column permutation matrix (Eq. (4.33))
param.P=zeros(param.N*param.M,param.N*param.M);
for j=1:param.N
    for i=1:param.M
        E=zeros(param.M,param.N);
        E(i,j)=1;
        param.P((j-1)*param.M+1:j*param.M,(i-1)*param.N+1:i*param.N)=E;
    end
end

% SNR range and Monte Carlo parameters
SNR_dB_range = 0:1:25;  % From 0dB to 25dB with step size 1dB
num_monte_carlo = 20;    % Number of Monte Carlo simulations

% Define OTFS types for comparison
otfs_types = {'OTFS', 'RZP-OTFS', 'RCP-OTFS', 'CP-OTFS', 'ZP-OTFS'};
BER = zeros(length(otfs_types), length(SNR_dB_range));  % Store average BER for each type at each SNR
current_snr_bers = zeros(length(otfs_types), num_monte_carlo);  % Store BER results for current SNR

% Simulate for each SNR value
for snr_idx = 1:length(SNR_dB_range)
    % Update current SNR value
    param.SNR_dB = SNR_dB_range(snr_idx);
    fprintf('\nSimulating SNR = %d dB:\n', param.SNR_dB);
    
    % Monte Carlo loop
    for mc = 1:num_monte_carlo
        fprintf('Monte Carlo iteration %d/%d:\n', mc, num_monte_carlo);
        [param.taps, param.l_i, param.k_i, param.g_i] = generate_channel_params(param);
        [param.noise, param.sigma_w_2] = generate_noise(param);
        % Simulate for each OTFS type
        for type_idx = 1:length(otfs_types)
            % Set current OTFS type
            param.otfs_type = otfs_types{type_idx};
            
            % OTFS Transmitter
            [tx_info_bits, tx_info_symbols,s]=otfs_transmitter(param);

            % OTFS Channel
            [r,G,H]=otfs_channel(s,tx_info_symbols, param);

            % OTFS Receiver
            [rx_info_bits, rx_info_symbols]=otfs_receiver(r, G, H, param);

            % Calculate BER for current iteration
            current_snr_bers(type_idx, mc) = sum(tx_info_bits ~= rx_info_bits) / length(tx_info_bits);
            
            % Display current progress
            fprintf('  %s: Current BER = %e\n', param.otfs_type, current_snr_bers(type_idx, mc));
        end
    end
    
    % Calculate average BER for each OTFS type at current SNR
    for type_idx = 1:length(otfs_types)
        BER(type_idx, snr_idx) = mean(current_snr_bers(type_idx, :));
        fprintf('==> %s: Average BER at SNR = %d dB is %e\n', ...
                otfs_types{type_idx}, SNR_dB_range(snr_idx), BER(type_idx, snr_idx));
    end
end

% Plotting the BER vs SNR curves
figure;
markers = {'o-', 's-', 'd-', '^-', 'v-'};  % Different marker styles
colors = {'b', 'r', 'g', 'm', 'c'};        % Different colors

hold on;
for i = 1:length(otfs_types)
    semilogy(SNR_dB_range, BER(i,:), markers{i}, 'LineWidth', 2, ...
             'MarkerSize', 8, 'Color', colors{i}, 'DisplayName', otfs_types{i});
end
hold off;

grid on;
xlabel('SNR (dB)');
ylabel('Average Bit Error Rate (BER)');
title('Comparison of Different OTFS Types Performance');
ylim([1e-4 1]);  % Set y-axis range
legend('Location', 'southwest');
ylabel('Average Bit Error Rate (BER)');
title('Comparison of Different OTFS Types Performance');
ylim([1e-4 1]);  % Set y-axis range
legend('Location', 'southwest');
