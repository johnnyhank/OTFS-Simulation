function [tx_info_bits, tx_info_symbols, s,P] = otfs_transmitter(param)
    M=param.M;
    N=param.N;
    % number of information symbols in one frame
    N_syms_per_frame=N*M;
    % number of information bits in one frame
    N_bits_per_frame=N_syms_per_frame*log2(param.mod_size);
    % generate random bits
    tx_info_bits=randi([0,1],N_bits_per_frame,1);
    % QAM modulation
    tx_info_symbols=qammod(tx_info_bits,param.mod_size,'gray','InputType','bit');
    % Generate the MxN OTFS delay-Doppler frame
    X=reshape(tx_info_symbols,M,N);
    % Vectorized OTFS frame information symbols
    x=reshape(X.',N*M,1);
    Im=eye(M);
    P=param.P; % row-column permutation matrix (Eq. (4.33))
    Fn=param.Fn; % normalized DFT matrix
    switch param.otfs_modulation_method
        case 1
            % Method 1 (Eqs. (4.17) and (4.20))
            X_tilda=X*Fn';
            % s=reshape(X_tilda,1,N*M);
            s=reshape(X_tilda,N*M,1);
        case 2
            % Method 2 (Eq. (4.35))
            s=P*kron(Im,Fn')*x;
        case 3
            % Method 3 (Eq. (4.35))
            s=kron(Fn',Im)*P*x;
        otherwise
            error('Unsupported OTFS modulation method selected.');
    end
end