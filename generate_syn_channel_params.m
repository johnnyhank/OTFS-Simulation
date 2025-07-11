function [taps, l_i, k_i, g_i] = generate_syn_channel_params(taps, l_max, k_max)
% generate synthetic channel parameters
% taps: the number of propagation paths
% l_max: maximum normalized delay
% k_max: maximum Doppler spread

% generate channel coefficients (Rayleigh fading) with uniform pdp(uniform power distribution)
g_i = sqrt(1/taps).*(sqrt(1/2) * (randn(1,taps)+1i*randn(1,taps)));
% generate delay taps uniformly from [0,l_max]
l_i = [randi([0,l_max],1,taps)];
l_i= l_i-min(l_i);
% generate Doppler taps (assuming uniform spectrum [-k_max,k_max])
k_i = k_max-2*k_max*rand(1,taps);
end