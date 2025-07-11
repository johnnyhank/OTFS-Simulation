function [taps, l_i, k_i, g_i] = generate_channel_params(param)
    M=param.M; % number of delay bins (subcarriers)
    N=param.N; % number of Doppler bins (time slots)
    P=param.P; % row-column permutation matrix (Eq. (4.33))
    Im=eye(M);
    Fn=param.Fn; % normalized DFT matrix
    max_cond_number = 1e8; % maximum condition number for the channel matrices H and G
    max_attempts = 10;  % max number of attempts 
    for attempt = 1:max_attempts
        switch param.channel_param_generation_mode
            case 'std'
                [taps, l_i, k_i, g_i] = generate_std_channel_params(param);
            case 'syn'
                taps=6;
                l_max=4;
                k_max=4;
                [taps, l_i, k_i, g_i] = generate_syn_channel_params(taps, l_max, k_max);
        end
        
        z=exp(1i*2*pi/N/M);
        delay_spread=max(l_i);
        % Generate discrete-time baseband channel gs in TDL form (Eq. (2.22))
        gs=zeros(delay_spread+1,N*M);
        for q=0:N*M-1
            for i=1:taps
                gs(l_i(i)+1,q+1)=gs(l_i(i)+1,q+1)+g_i(i)*z^(k_i(i)*(q-l_i(i)));
            end
        end
        % Generate discrete-time baseband channel matrix (Eq. (4.38))
        G=zeros(N*M,N*M);
        for q=0:N*M-1
            for ell=0:delay_spread
                if(q>=ell)
                    G(q+1,q-ell+1)=gs(ell+1,q+1);
                end
            end
        end
        
        % generate delay-Doppler channel matrix (Eq. (6.1))
        H=kron(Im,Fn)*(P'*G*P)*kron(Im,Fn');
        
        % check condition numbers
        cond_G = cond(G'*G);
        cond_H = cond(H'*H);
        
        if cond_G < max_cond_number && cond_H < max_cond_number
            break;
        end
        
        if attempt == max_attempts
            warning('max_attempts reached, using the last valid channel matrix');
        end
    end
end