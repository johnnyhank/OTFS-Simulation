function [noise, sigma_w_2] = generate_noise(param)
    % Add AWGN noise to the received signal
    % calculate average QAM symbol energy
    M = param.M; % number of subcarriers
    N = param.N; % number of time slots
    mod_size = param.mod_size; % modulation size
    Es = mean(abs(qammod(0:mod_size-1,mod_size).^ 2));
    % SNR=Es/noise power
    SNR=10.^(param.SNR_dB/10);
    % noise power
    sigma_w_2=Es/SNR;
    % generate Gaussian noise samples with variance=sigma_w_2
    noise = sqrt(sigma_w_2/2)*(randn(N*M,1) + 1i*randn(N*M,1));
end