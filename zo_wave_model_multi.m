%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Surface Wave-Aerodynamic Roughness Length Model for Air-Sea Interactions %%%
%%%                                                                          %%%
%%%     - This script loads the surface geometric information, calculates    %%%
%%%       the spatial and temporal gradients using first-order finite        %%%
%%%       difference, calculates the zo for sgs surfaces, calculates         %%%
%%%       Lambda through Newton-Raphson iteration                            %%%
%%%       and finally calculates the surface roughness.                      %%%
%%%                                                                          %%%
%%%     - This script is specifically for multiscale waves cases             %%%
%%%                                                                          %%%
%%% Script developed by: Manuel Ayala                                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc

%%% Loading the surface files containing the following:
% - 2 realizations of the wavy surface (eta_t2, eta_t1)
% - Turbulent Reynolds number (Ret)
% - Friction velocity (u_star)
% - Phillips constant (alpha_p)
% - Peak wavenumber of wave (kp) in meters
% - Domain size (Lx, Ly) and resolution (nx, ny)
% - Significant wave height (Hs) in meters
load("Surfaces/S4.mat");
lambdap = (2*pi)/kp;

%%% Calculating the spatial gradients, temporal gradients and the wave velocity
dt = 0.001; 
detadt = (eta_t2 - eta_t1)/dt;          
detadx = diff(eta_t2,1,2)/(Lx/nx);              
detady = diff(eta_t2,1,1)/(Ly/ny);              
detadx(:,nx) = -eta_t2(:,nx)/(Lx/nx);
detady(ny,:) = -eta_t2(ny,:)/(Ly/ny);
grad = sqrt(detadx.^2 +  detady.^2);        
nhatx = detadx./grad;
nhaty = detady./grad;
cx = - detadt.*detadx./grad.^2;
cy = - detadt.*detady./grad.^2;
k_min = 0.25*(kp);
threshold_wave_speed = sqrt(9.81/k_min);
cx = min(max(cx, -threshold_wave_speed), threshold_wave_speed);
cy = min(max(cy, -threshold_wave_speed), threshold_wave_speed);
cx = cx/u_star;
cy = cy/u_star;

%%% Calculating delta (reference height) 
% - delta is normalized by lambdap
% - Normalizing a surface realization by peak wavelength (lambdap)
eta = eta_t2/lambdap;
eta_mean = mean(eta(:));
eta_p = eta - eta_mean;
ramp_eta_p = eta_p .* (eta_p >= 0);
ramp__eta_p8 = ramp_eta_p.^8;
H_p = (mean(ramp__eta_p8(:)))^(1/8);
delta = 3*H_p;

%%% Calculating the rms of sgs surfaces, surface roughness for 
%%% sgs surfaces and normalizing by lambdap
k_filter = sqrt((pi/(Lx/nx))^2 + (pi/(Ly/ny))^2);
eta_sgs = sqrt((0.2*alpha)/kp^2*(1 - exp(-5/4*(kp/k_filter)^2)));
kappa = 0.4;
zo_u = eta_sgs*exp(-kappa*8.5)/lambdap;

%%% Setting the Newton-Rahpson iteration 
% - flag_sgs = 1 is to activate effects of 
%   unresolved small-scale sea surface features
max_iterations = 10^3;
tolerance = 1e-6;
initial_guess = 20;
flag_sgs = 1;

%%% Calculating U through iteration and calculating Lambda
[U] = solve_U(detadx, grad, nhatx, nhaty, cx, cy,  delta, Ret, zo_u, ...
    flag_sgs, max_iterations, tolerance, initial_guess);
lambda = 1/U^2;

%%% Calculating surface roughness for wavy surface (zo) normalized 
%%% by the wave significant wave height 
zo_wave_model = delta*exp(-kappa*(1/sqrt(lambda)));
zo_wave_model = zo_wave_model*lambdap/Hs_measured;
fprintf('zo_wave_model: %.6f\n', zo_wave_model);
