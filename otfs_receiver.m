function [rx_info_bits, rx_info_symbols] = otfs_receiver(r, G, H, param)
N=param.N; % number of Doppler bins (time slots)
M=param.M; % number of delay bins (subcarriers)
mod_size=param.mod_size; % modulation size
Fn=param.Fn; % normalized DFT matrix
P=param.P; % row-column permutation matrix (Eq. (4.33))
% % Add AWGN noise to the received signal
% % calculate average QAM symbol energy
% Es = mean(abs(qammod(0:mod_size-1,mod_size).^ 2));
% % SNR=Es/noise power
% SNR=10.^(param.SNR_dB/10);
% % noise power
% sigma_w_2=Es/SNR;
% % generate Gaussian noise samples with variance=sigma_w_2
% noise = sqrt(sigma_w_2/2)*(randn(N*M,1) + 1i*randn(N*M,1));
% % add AWGN to the received signal
sigma_w_2 = param.sigma_w_2; % noise power
noise= param.noise; % pre-generated noise
r=r+noise;

% OTFS demodulation
switch param.otfs_demod_method
    case 1
        % Method 1 (Eqs. (4.24) and (4.27))
        Y_tilda=reshape(r,M,N);
        Y=Y_tilda*Fn;
    case 2
        % Method 2 (Eq. (4.35))
        y=kron(eye(M),Fn)*(P.')*r;
        Y=reshape(y,N,M).';
    case 3
        % Method 3 (Eq. (4.35))
        y=(P.')*kron(Fn,eye(M))*r;
        Y=reshape(y,N,M).';
    otherwise
        error('Unsupported OTFS demodulation method');
end

% LMMSE detection
switch param.otfs_detection_method
    case 'LMMSE_delay-Doppler'
        % Add condition number check
        cond_number = cond(H'*H);
        if cond_number > 1e10
            warning('Delay-Doppler channel matrix condition number is too high: %e', cond_number);
        end
        % vectorize Y
        y=reshape(Y.',N*M,1);
        % estimated delay-Doppler matrix (Eq. (6.18))
        x_hat=(H'*H+sigma_w_2)^(-1)*(H'*y);
        % the following two lies are added to match the transmitter code of OTFS frame generation 
        X_hat = reshape(x_hat, N, M).';
        rx_info_symbols = X_hat(:);
        % QAM demodulation
        % rx_info_bits =qamdemod(x_hat,mod_size,'gray','OutputType','bit')
        rx_info_bits = qamdemod(rx_info_symbols, mod_size, 'gray', 'OutputType', 'bit');
    case 'LMMSE_time'
        % Add before LMMSE detection
        cond_number = cond(G'*G);
        if cond_number > 1e10
            warning('Channel matrix condition number is too high: %e', cond_number);
        end
        % estimated time domain samples (Eq. (6.19))
        s_hat=(G'*G+sigma_w_2)^(-1)*(G'*r);
        % MxN estimated delay-Doppler symbols (using Method 1 in code 12)
        X_hat_tilda=reshape(s_hat,M,N);
        X_hat=X_hat_tilda*Fn;
        x_hat=reshape(X_hat.',N*M,1);
        % the following two lines are added to match the transmitter code of OTFS frame generation
        X_hat = reshape(x_hat, N, M).';
        rx_info_symbols = X_hat(:);
        % QAM demodulation
        % rx_info_bits=qamdemod(x_hat,mod_size,'gray','OutputType','bit')
        rx_info_bits = qamdemod(rx_info_symbols, mod_size, 'gray', 'OutputType', 'bit');
    otherwise
        error('Unsupported OTFS detection method');
end

