%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Surface Wave-Aerodynamic Roughness Length Model for Air-Sea Interactions %%%
%%%                                                                          %%%
%%%     - This script loads the surface geometric information, calculates    %%%
%%%       the spatial and temporal gradients using first-order finite        %%%
%%%       difference, calculates Lambda through Newton-Raphson iteration     %%%
%%%       and finally calculates the surface roughness.                      %%%
%%%                                                                          %%%
%%%     - This script is specifically for monochromatic waves  cases         %%%
%%%                                                                          %%%
%%% Script developed by: Manuel Ayala                                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc

%%% Loading the surface files containing the following:
% - 2 realizations of the wavy surface (eta)
% - Turbulent Reynolds number (Ret)
% - Friction velocity (u_star)
% - Height of boundary layer, used for Ret (h)
% - Amplitude of wave, normalized by h (amp)
% - Domain size (Lx, Ly) and resolution (nx, ny)
load("W1.mat");

%%% Calculating the spatial gradients, temporal gradients and the wave velocity
dt = 1e-3; 
detadt = (eta_t2 - eta_t1)/dt;         
detadx = diff(eta_t2,1,1)/(Lx/nx);            
detady = diff(eta_t2,1,2)/(Ly/ny);           
detadx(nx,:) = -eta_t2(nx,:)/(Lx/nx);
detady(:,ny) = -eta_t2(:,ny)/(Ly/ny);
grad = sqrt(detadx.^2 + detady.^2);
nhatx = detadx./grad;
nhaty = detady./grad;
cx = - detadt.*detadx./grad.^2;
cy = - detadt.*detady./grad.^2;
cx = cx/u_star;
cy = cy/u_star;

%%% Calculating delta (reference height) 
% - delta is normalized by h
delta = 2.55*amp;

%%% Setting the Newton-Rahpson iteration 
% - flag_sgs = 1 is to activate effects of 
%   unresolved small-scale sea surface features
max_iterations = 10^3;
tolerance = 1e-6;
initial_guess = 20;
flag_sgs = 0;
zo_u = 0;

%%% Calculating U through iteration and calculating Lambda
[U] = solve_U(detadx, grad, nhatx, nhaty, cx, cy,  delta, Ret, zo_u, ...
    flag_sgs, max_iterations, tolerance, initial_guess);
lambda = 1/U^2;
kappa = 0.4;

%%% Calculating surface roughness for wavy surface (zo) normalized 
%%% by wave amplitude
zo_wave_model = delta*exp(-kappa*(1/sqrt(lambda)));
zo_wave_model = zo_wave_model*(1/amp);
fprintf('zo_wave_model: %.6f\n', zo_wave_model);
