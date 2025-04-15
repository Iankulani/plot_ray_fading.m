% plot_ray_fading.m
clear;
clf;

% Channel parameters
fc = 9e8; % Carrier frequency (Hz)
fs = 5e4; % Sampling frequency (Hz)
speed_kmh = 120; % Speed in km/h
Ts = 1/fs; % Sampling period

% Convert speed from km/h to m/s
v_ms = speed_kmh / 3.6; 

% Wavelength (m)
wl_m = 3e8 / fc;

% Power Delay Profile (PDP) and delay times (ns)
PDP_dB = [0, -1, -9, -10, -15, -20]; % PDP in dB
t_ns = [0, 310, 710, 1090, 1730, 2510]; % Delay times in ns

% Base Station (BS) and Mobile Station (MS) parameters
BS_theta_LOS_deg = 0; % BS Line of Sight (LOS) angle in degrees
MS_theta_LOS_deg = 0; % MS LOS angle in degrees
BS_AS_deg = 2; % BS Azimuth Spread (AS) in degrees (Laplacian PAS)
BS_AoD_deg = 50 * ones(size(PDP_dB)); % BS Angle of Departure (AoD) in degrees
MS_AS_deg = 35; % MS Azimuth Spread (AS) in degrees (Laplacian PAS)
DoT_deg = 22.5; % Doppler shift angle (degrees)
MS_AoA_deg = 67.5 * ones(size(PDP_dB)); % MS Angle of Arrival (AoA) in degrees

% Generate phase information for subray
[BS_theta_deg, MS_theta_deg, BS_PHI_rad] = gen_phase(BS_theta_LOS_deg, ...
    BS_AS_deg, BS_AoD_deg, MS_theta_LOS_deg, MS_AS_deg, MS_AoA_deg);

% Convert PDP from dB to linear scale
PDP = dB2w(PDP_dB);

% Time vector for simulation
t = [0:1e4-1] * Ts; % Time vector in seconds

% Generate Rayleigh fading channel response
h = ray_fading(20, PDP, BS_PHI_rad, MS_theta_deg, v_ms, DoT_deg, wl_m, t);

% Plot the magnitude of the Rayleigh fading channel response
figure;
plot(t, 10 * log10(abs(h(1,:))));
title(['Ray Channel Model, f_c = ', num2str(fc), ' Hz, T_s = ', num2str(Ts), ' s']);
xlabel('Time [s]');
ylabel('Magnitude [dB]');
grid on;

% --------- Placeholder functions to be defined ---------
% Function to generate the phase information for subray
function [BS_theta_deg, MS_theta_deg, BS_PHI_rad] = gen_phase(BS_theta_LOS_deg, ...
    BS_AS_deg, BS_AoD_deg, MS_theta_LOS_deg, MS_AS_deg, MS_AoA_deg)
    % This function should compute the phase information of the subray
    % based on given parameters. This is a placeholder function.
    BS_theta_deg = BS_theta_LOS_deg; % Placeholder value
    MS_theta_deg = MS_theta_LOS_deg; % Placeholder value
    BS_PHI_rad = rand(size(BS_AoD_deg)); % Random phase generation (for placeholder)
end

% Function to convert dB to linear scale
function linear = dB2w(dB)
    linear = 10 .^ (dB / 10);
end

% Function to simulate Rayleigh fading (Placeholder)
function h = ray_fading(N, PDP, BS_PHI_rad, MS_theta_deg, v_ms, DoT_deg, wl_m, t)
    % This function simulates Rayleigh fading based on given parameters.
    % It's a placeholder for a more detailed fading model.
    h = randn(N, length(t)); % Placeholder: random Rayleigh fading channel
end

