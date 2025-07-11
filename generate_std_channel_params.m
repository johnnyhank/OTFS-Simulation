function [taps, l_i, k_i, g_i] = generate_std_channel_params(param)
% maximum user equipment speed
max_UE_speed_kmh=param.max_UE_speed_kmh;
max_UE_speed = max_UE_speed_kmh*(1000/3600);
% maximum Doppler spread (one-sided)
nu_max = (max_UE_speed*param.fc)/(param.c);
% maximum normalized Doppler spread (one-sided)
k_max = nu_max/param.Doppler_resolution;
% choose delays and pdp for specific 3GPP model
switch upper(param.model)
    case 'EPA'
        delays = [0, 30, 70, 90, 110, 190, 410]*1e-9;
        pdp = [0.0, -1.0, -2.0, -3.0, -8.0, -17.2, -20.8];
    case 'EVA'
        delays= [0, 30, 150, 310, 370, 710, 1090, 1730, 2510]*1e-9;
        pdp= [0.0, -1.5, -1.4, -3.6, -0.6, -9.1, -7.0, -12.0, -16.9];
    case 'ETU'
        delays = [0, 50, 120, 200, 230, 500, 1600, 2300, 5000]*1e-9;
        pdp = [-1.0, -1.0, -1.0, 0.0, 0.0, 0.0, -3.0, -5.0, -7.0];
    otherwise
        error('Unsupported channel model: %s', model);
end
% dB to linear scale
pdp_linear = 10.^(pdp/10);
% normalization
pdp_linear = pdp_linear/sum(pdp_linear);
% number of propagation paths (taps)
taps=length(pdp);
% generate channel coefficients (Rayleigh fading)
g_i = sqrt(pdp_linear).*(sqrt(1/2) * (randn(1,taps)+1i*randn(1,taps)));
% generate delay taps (assuming integer delay taps)
l_i=round(delays./param.delay_resolution);
% Generate Doppler taps (assuming Jakes spectrum)
k_i = (k_max*cos(2*pi*rand(1,taps)));
end
