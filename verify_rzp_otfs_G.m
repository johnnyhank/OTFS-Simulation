% filepath: /Users/zhenglihan/Documents/MATLAB/OTFS/verify_RZP_OTFS.m

% Test parameters
M = 4;          % Block size
N = 2;          % Number of blocks
delay_spread = 2;% Maximum delay
taps = 3;       % Number of channel taps

% Channel parameters
l_i = [0; 1; 2];           % Delays
k_i = [1; 2; 3];          % Doppler shifts
g_i = [0.8; 0.5; 0.3];    % Channel gains
z = exp(1j*2*pi/(N*M));   % Phase factor

% Generate test input signal
s = randn(N*M,1) + 1j*randn(N*M,1);

% Generate channel response in TDL form
gs = zeros(delay_spread+1, N*M);
for q = 0:N*M-1
    for i = 1:taps
        gs(l_i(i)+1,q+1) = gs(l_i(i)+1,q+1) + g_i(i)*z^(k_i(i)*(q-l_i(i)));
    end
end

% Method 1: Matrix multiplication approach
% Generate RZP-OTFS channel matrix
G_rzp=zeros(N*M,N*M);
for n=0:N-1  % Block index
    for m=0:M-1  % In-block index
        q = m + n*M;  
        for ell=0:delay_spread
            if n == 0  % First block with zero padding
                if m >= ell  
                    G_rzp(q+1,q-ell+1)=gs(ell+1,q+1);
                end
            else  % Other blocks
                if q >= ell
                    G_rzp(q+1,q-ell+1)=gs(ell+1,q+1);
                end
            end
        end
    end
end

% Calculate output using matrix multiplication
r1 = G_rzp * s;

% Method 2: Direct implementation of equation (4.68)
r2 = zeros(N*M,1);
for n = 0:N-1
    for m = 0:M-1
        q = m + n*M;
        r2(q+1) = 0;
        
        if n == 0  % First block
            for ell = 0:delay_spread
                if m >= ell
                    r2(q+1) = r2(q+1) + gs(ell+1,q+1)*s(q-ell+1);
                end
            end
        else  % Other blocks
            for ell = 0:delay_spread
                if m >= ell  % In-block samples
                    r2(q+1) = r2(q+1) + gs(ell+1,q+1)*s(q-ell+1);
                else  % Inter-block interference
                    mod_idx = mod(m-ell,M);
                    prev_block_idx = mod_idx+(n-1)*M;
                    if prev_block_idx >= 0
                        r2(q+1) = r2(q+1) + gs(ell+1,q+1)*s(prev_block_idx+1);
                    end
                end
            end
        end
    end
end

% Compare results
error = norm(r1 - r2)/norm(r1);
fprintf('Relative error: %e\n', error);

if error < 1e-10
    fprintf('Verification passed: Matrix implementation matches direct calculation!\n');
else
    fprintf('Warning: Significant difference between methods!\n');
    disp('First few elements of matrix method result:');
    disp(r1(1:min(5,end)));
    disp('First few elements of direct calculation result:');
    disp(r2(1:min(5,end)));
end