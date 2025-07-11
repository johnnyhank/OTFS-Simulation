function [r,G,H] = otfs_channel(s, tx_info_symbols, param)
    M=param.M; % number of delay bins (subcarriers)
    N=param.N; % number of Doppler bins (time slots)
    P=param.P; % row-column permutation matrix (Eq. (4.33))
    Im=eye(M);
    Fn=param.Fn; % normalized DFT matrix
    X=reshape(tx_info_symbols,M,N);
    X_tilda=reshape(s,M,N); 
    taps=param.taps; % number of taps in the channel
    l_i=param.l_i; % normalized delays of the taps
    k_i=param.k_i; % normalized Doppler shifts of the taps
    g_i=param.g_i; % complex gains of the taps
    % max_cond_number = 1e8; % maximum condition number for the channel matrices H and G
    % max_attempts = 10;  % max number of attempts 
    % for attempt = 1:max_attempts
    %     switch param.channel_param_generation_mode
    %         case 'std'
    %             [taps, l_i, k_i, g_i] = generate_std_channel_params(param);
    %         case 'syn'
    %             taps=6;
    %             l_max=4;
    %             k_max=4;
    %             [l_i, k_i, g_i] = generate_syn_channel_params(taps, l_max, k_max);
    %     end
        
    %     z=exp(1i*2*pi/N/M);
    %     delay_spread=max(l_i);
    %     % Generate discrete-time baseband channel gs in TDL form (Eq. (2.22))
    %     gs=zeros(delay_spread+1,N*M);
    %     for q=0:N*M-1
    %         for i=1:taps
    %             gs(l_i(i)+1,q+1)=gs(l_i(i)+1,q+1)+g_i(i)*z^(k_i(i)*(q-l_i(i)));
    %         end
    %     end
    %     % Generate discrete-time baseband channel matrix (Eq. (4.38))
    %     G=zeros(N*M,N*M);
    %     for q=0:N*M-1
    %         for ell=0:delay_spread
    %             if(q>=ell)
    %                 G(q+1,q-ell+1)=gs(ell+1,q+1);
    %             end
    %         end
    %     end
        
    %     % generate delay-Doppler channel matrix (Eq. (6.1))
    %     H=kron(Im,Fn)*(P'*G*P)*kron(Im,Fn');
        
    %     % check condition numbers
    %     cond_G = cond(G'*G);
    %     cond_H = cond(H'*H);
        
    %     if cond_G < max_cond_number && cond_H < max_cond_number
    %         break;
    %     end
        
    %     if attempt == max_attempts
    %         warning('max_attempts reached, using the last valid channel matrix');
    %     end
    % end
    z=exp(1i*2*pi/N/M);
    delay_spread=max(l_i);
    switch param.otfs_type
        case 'OTFS'
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
            % generate delay-time channel matrix (Eq. (4.55))
            H_tilda=P*G*P.';
            % generate delay-Doppler channel matrix (Eq. (6.1))
            H=kron(Im,Fn)*(P'*G*P)*kron(Im,Fn');
            % Generate r by passing the Tx signal through the channel
            switch param.res_calc_method
                case 1
                    % Method 1: Using the TDL model (Eq. (4.36))
                    r=zeros(N*M,1);
                        for q=0:N*M-1
                            for ell=0:(delay_spread-1)
                                if(q>=ell)
                                    r(q+1)=r(q+1)+gs(ell+1,q+1)*s(q-ell+1);
                                end
                            end
                        end
                case 2
                    % Method 2: Using the time-domain channel matrix (G) (Eq. (4.37))
                    r=G*s;
                case 3
                    % Method 3: Using the delay-time channel matrix (H_tilda) (Eq. (4.54))
                    x_tilda=reshape(X_tilda.',N*M,1);
                    y_tilda=H_tilda*x_tilda;
                    r=P*y_tilda;
                case 4
                    % Method 4: Using the delay-Doppler channel matrix (H) (Eq. (4.59))
                    x=reshape(X.',N*M,1);
                    y=H*x;
                    r=P*kron(Im,Fn')*y;
            otherwise
                    error('Unsupported response calculation method selected.');
            end
        case 'RZP-OTFS'
            % Generate discrete-time baseband channel in TDL form (Eq. (2.22))
            gs=zeros(delay_spread+1,N*M);
            for q=0:N*M-1
                for i=1:taps
                    gs(l_i(i)+1,q+1)=gs(l_i(i)+1,q+1)+g_i(i)*z^(k_i(i)*(q-l_i(i)));
                end
            end
            % Generate discrete-time baseband channel matrix (Eq. (4.38))
            G_rzp=zeros(N*M,N*M);
            for q=0:N*M-1
                for ell=0:delay_spread
                    if(q>=ell)
                        G_rzp(q+1,q-ell+1)=gs(ell+1,q+1);
                    end
                end
            end
            % generate received signal after discarding CP
            % fprintf('s的维度为: [%d, %d]\n', size(s,1), size(s,2));
            % fprintf('G_rzp的维度为: [%d, %d]\n', size(G_rzp,1), size(G_rzp,2));
            r=G_rzp*s;
            G =G_rzp; % Assign G_rzp to G for consistency
            % generate delay-Doppler channel matrix (Eq. (6.1))
            H=kron(Im,Fn)*(P'*G*P)*kron(Im,Fn');
            % fprintf('r的维度为: [%d, %d]\n', size(r,1), size(r,2));
        case 'RCP-OTFS'
            % Generate discrete-time baseband channel in TDL form (Eq. (2.22))
            gs=zeros(delay_spread+1,N*M);
            for q=0:N*M-1
                for i=1:taps
                    gs(l_i(i)+1,q+1)=gs(l_i(i)+1,q+1)+g_i(i)*z^(k_i(i)*(q-l_i(i)));
                end
            end
            % Generate discrete-time baseband channel matrix (Eq. (4.83))
            G_rcp=zeros(N*M,N*M);
            for q=0:N*M-1
                for ell=0:delay_spread
                    G_rcp(q+1,mod(q-ell,N*M)+1)=gs(ell+1,q+1);
                end
            end
            % generate received signal after discarding CP
            r=G_rcp*s;
            G =G_rcp; % Assign G_rcp to G for consistency
            % generate delay-Doppler channel matrix (Eq. (6.1))
            H=kron(Im,Fn)*(P'*G*P)*kron(Im,Fn');
        case 'CP-OTFS'
            l_cp=delay_spread;
            % Generate discrete-time baseband channel in TDL form (Eq. (2.22))
            gs=zeros(delay_spread+1,N*(M+l_cp));
            for q=0:N*(M+l_cp)-1
                for i=1:taps
                    gs(l_i(i)+1,q+1)=gs(l_i(i)+1,q+1)+g_i(i)*z^(k_i(i)*(q-l_i(i)));
                end
            end
            % Generate discrete-time baseband channel matrix (Eq. (4.93))
            G_cp=zeros(N*M,N*M);
            for n=0:N-1
                for m=0:M-1
                    for ell=0:delay_spread
                        G_cp(m+n*M+1,n*M+mod(m-ell,M)+1)=gs(ell+1,m+n*M+l_cp+1);
                    end
                end
            end
            % generate received signal after discarding CP per block
            r=G_cp*s;
            G =G_cp; % Assign G_cp to G for consistency
            % generate delay-Doppler channel matrix (Eq. (6.1))
            H=kron(Im,Fn)*(P'*G*P)*kron(Im,Fn');
        case 'ZP-OTFS'
            l_zp=delay_spread;
            % Generate discrete-time baseband channel in TDL form (Eq. (2.22))
            gs=zeros(delay_spread+1,N*(M+l_zp));
            for q=0:N*(M+l_zp)-1
                for i=1:taps
                    gs(l_i(i)+1,q+1)=gs(l_i(i)+1,q+1)+g_i(i)*z^(k_i(i)*(q-l_i(i)));
                end
            end
            % Generate discrete-time baseband channel matrix (Eq. (4.109))
            G_zp=zeros(N*M,N*M);
            for n=0:N-1
                for m=0:M-1
                    for ell=0:delay_spread
                        if(m>=ell)
                            G_zp(m+n*M+1,m+n*M-ell+1)=gs(ell+1,m+n*M+l_zp+1);
                        end
                    end
                end
            end
            % generate received signal after discarding ZP per block
            r=G_zp*s;
            G =G_zp; % Assign G_zp to G for consistency
            % generate delay-Doppler channel matrix (Eq. (6.1))
            H=kron(Im,Fn)*(P'*G*P)*kron(Im,Fn');
        otherwise
            error('Unsupported OTFS type selected.');
    end
end