%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Surface Wave-Aerodynamic Roughness Length Model for Air-Sea Interactions %%%
%%%                                                                          %%%
%%%     - This script loads the surface geometric information, calculates    %%%
%%%       the spatial and temporal gradients using first-order finite        %%%
%%%       difference, calculates the zo for sgs surfaces, calculates         %%%
%%%       Lambda through Newton-Raphson iteration                            %%%
%%%       and finally calculates the surface roughness.                      %%%
%%%                                                                          %%%
%%%     - If using multiscale wave case, the following is loaded:            %%%
%%%             * - 2 realizations of the wavy surface (eta_t2, eta_t1)      %%%
%%%             * - Turbulent Reynolds number (Ret)                          %%%
%%%             * - Friction velocity (u_star)                               %%%
%%%             * - Phillips constant (alpha_p)                              %%%
%%%             * - Peak wavenumber of wave (kp) in meters                   %%%
%%%             * - Domain size (Lx, Ly) and resolution (nx, ny)             %%%
%%%             * - Significant wave height (Hs) in meters                   %%%
%%%                                                                          %%%
%%%     - If using monochromatic wave case, the following is loaded:         %%%
%%%             * - 2 realizations of the wavy surface (eta)                 %%%
%%%             * - Turbulent Reynolds number (Ret)                          %%%
%%%             * - Friction velocity (u_star)                               %%%
%%%             * - Height of boundary layer, used for Ret (h)               %%%
%%%             * - Amplitude of wave, normalized by h (amp)                 %%%
%%%             * - Domain size (Lx, Ly) and resolution (nx, ny)             %%%
%%%                                                                          %%%
%%%                                                                          %%%
%%% Script developed by: Manuel Ayala                                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc

%%% Select the wave type and the surface case:
% - Monochromatic --> "mono"
% - Multiscale    --> "multi"
wave_type = 'multi';
load("Surfaces/Y2.mat");

%%% Calculating the spatial gradients, temporal gradients and the wave velocity
%  - If multiscale case selected, a clipping of cx and cy will be applied
dt = 0.001; 
detadt = (eta_t2 - eta_t1)/dt;          
detadx = diff(eta_t2,1,1)/(Lx/nx);              
detady = diff(eta_t2,1,2)/(Ly/ny);              
detadx(nx,:) = detadx(nx-1,:);
detady(:,ny) = detady(:,ny-1);
grad = sqrt(detadx.^2 +  detady.^2);        
nhatx = detadx./grad;
nhaty = detady./grad;
cx = - detadt.*detadx./grad.^2;
cy = - detadt.*detady./grad.^2;
if strcmp(wave_type, 'multi')
k_min = 0.25*(kp);
threshold_wave_speed = sqrt(9.81/k_min);
cx = min(max(cx, -threshold_wave_speed), threshold_wave_speed);
cy = min(max(cy, -threshold_wave_speed), threshold_wave_speed);
end
cx = cx/u_star;
cy = cy/u_star;


%%% Calculating delta (reference height) 
% - delta is normalized by lambdap
% - The surface realization will be normalized by peak wavelength (multiscale)
%   or by wave amplitude (monochromatic)
if strcmp(wave_type, 'multi')
lambdap = (2*pi)/kp;
eta = eta_t2/lambdap;
else
eta = eta_t2/h;
end

eta_mean = mean(eta(:));
eta_p = eta - eta_mean;
ramp_eta_p = eta_p .* (eta_p >= 0);
ramp__eta_p8 = ramp_eta_p.^8;
H_p = (mean(ramp__eta_p8(:)))^(1/8);
delta = 3*H_p;

%%% Calculating the rms of sgs surfaces, surface roughness for 
%%% sgs surfaces and normalizing by lambdap
kappa = 0.4;
if strcmp(wave_type, 'multi')
k_filter = sqrt((pi/(Lx/nx))^2 + (pi/(Ly/ny))^2);
eta_sgs = sqrt((0.2*alpha_p)/kp^2*(1 - exp(-5/4*(kp/k_filter)^2)));
zo_u = eta_sgs*exp(-kappa*8.5)/lambdap;
else
zo_u = 0.0;
end

%%% Setting the Newton-Rahpson iteration 
% - flag_sgs = 1 is to activate effects of 
%   unresolved small-scale sea surface features
% - flag_sgs = 0 is to have a smooth wave.
max_iterations = 10^3;
tolerance = 1e-6;
initial_guess = 20;
if strcmp(wave_type, 'multi')
flag_sgs = 1;
else
flag_sgs = 0;
end

%%% Calculating U through iteration and calculating Lambda
[U] = solve_U(detadx, grad, nhatx, nhaty, cx, cy,  delta, Ret, zo_u, ...
    flag_sgs, max_iterations, tolerance, initial_guess);
lambda = 1/U^2;

%%% Calculating surface roughness for wavy surface (zo) normalized 
%%% by the wave significant wave height (multiscale)
%%% or by wave amplitude (monochromatic)
zo_wave_model = delta*exp(-kappa*(1/sqrt(lambda)));
if strcmp(wave_type, 'multi')
zo_wave_model = zo_wave_model*lambdap/Hs_measured;
else
zo_wave_model = zo_wave_model*(1/amp);
end
fprintf('zo_wave_model: %.6f\n', zo_wave_model);
