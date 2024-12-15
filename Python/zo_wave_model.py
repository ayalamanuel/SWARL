import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from typing import Tuple, Optional

def solve_U(detadx: np.ndarray, grad: np.ndarray, nx: np.ndarray, ny: np.ndarray, 
            c_x: np.ndarray, c_y: np.ndarray, delta: float, Ret: float, zo: float, 
            flag: int, max_iterations: int, tolerance: float, initial_guess: float) -> float:
    """
    Newton-Raphson iteration to estimate U which will be used to estimate Lambda.
    
    Parameters:
        detadx: Spatial gradient in x direction
        grad: Total gradient magnitude
        nx, ny: Surface normal components
        c_x, c_y: Phase velocities (normalized by friction velocity)
        delta: Reference height (normalized)
        Ret: Turbulent Reynolds number
        zo: Surface roughness for sub-grid scale surfaces
        flag: Toggle for unresolved small-scale effects
        max_iterations: Maximum number of Newton-Raphson iterations
        tolerance: Convergence tolerance
        initial_guess: Initial value for U
    
    Returns:
        U: Converged velocity value
    """
    U = initial_guess
    
    def function_and_derivative(U_current: float) -> Tuple[float, float]:
        # Calculate Reynolds number
        Re = U_current * delta * Ret
        
        # Calculate c_fs
        c_fs = 0.0288 * Re**(-1/5) * (1 + 577*Re**(-6/5))**(2/3)
        
        # Calculate c_f with simplified formula
        c_f = 2 * ((1/2*c_fs)**3 + flag*((1/0.4)*np.log(delta/zo))**(-6))**(1/3)
        
        # Derivative of c_fs with respect to Re
        dc_fs_dRe = 0.0288 * (
            (-1/5) * Re**(-6/5) * (1 + 577*Re**(-6/5))**(2/3) +
            Re**(-1/5) * (2/3) * (1 + 577*Re**(-6/5))**(-1/3) *
            577 * (-6/5) * Re**(-11/5)
        )
        
        # Complete derivative of c_f with respect to Re
        dc_f_dRe = 2 * (1/3) * ((1/2*c_fs)**3 + flag*((1/0.4)*np.log(delta/zo))**(-6))**(-2/3) * \
                   (3*(1/2*c_fs)**2 * (1/2)*dc_fs_dRe)
        
        # Convert to derivative with respect to U
        dc_f_dU = dc_f_dRe * (delta * Ret)
        
        # Calculate W
        alpha = np.arctan(grad)
        cx_term = 1 - c_x/U_current
        term1 = (nx**2) * (cx_term**2)
        term2 = (ny**2) * ((c_y/U_current)**2)
        
        numer = nx * cx_term
        dir_term = (numer + np.abs(numer))/numer
        
        W_avg = np.mean((alpha/(np.pi + alpha)) * (term1 + term2) * detadx * 0.5 * dir_term)
        W = W_avg + 0.5*c_f
        
        f = U_current - 1/np.sqrt(W)
        
        # Calculate dW/dU
        dcx_term_dU = c_x/U_current**2
        dterm1_dU = 2 * (nx**2) * cx_term * dcx_term_dU
        dterm2_dU = -2 * (ny**2) * (c_y**2/U_current**3)
        
        dnumer_dU = nx * dcx_term_dU
        ddir_dU = (dnumer_dU * numer - (numer + np.abs(numer)) * dnumer_dU) / numer**2
        
        dW_avg_dU = np.mean((alpha/(np.pi + alpha)) * (
            (dterm1_dU + dterm2_dU) * detadx * 0.5 * dir_term +
            (term1 + term2) * detadx * 0.5 * ddir_dU))
        
        dW_dU = dW_avg_dU + 0.5*dc_f_dU
        df = 1 + (1/2) * W**(-3/2) * dW_dU
        
        return f, df
    
    # Newton-Raphson iteration
    for iter in range(max_iterations):
        U_prev = U
        f, df = function_and_derivative(U)
        U = U - f/df
        error = abs(U - U_prev) / abs(U)
        
        if error < tolerance:
            break
        if iter == max_iterations - 1:
            print('Warning: Maximum iterations reached without convergence')
    
    return U

def calculate_wave_model_mono(data_file: str) -> float:
    """
    Calculate surface roughness for monochromatic waves.
    
    Parameters:
        data_file: Path to .npz file containing wave data
        
    Returns:
        zo_wave_model: Normalized surface roughness
    """
    # Load data
    data = np.load(data_file)
    eta_t1 = data['eta_t1']
    eta_t2 = data['eta_t2']
    Ret = data['Ret']
    u_star = data['u_star']
    h = data['h']
    amp = data['amp']
    Lx = data['Lx']
    Ly = data['Ly']
    nx = data['nx']
    ny = data['ny']
    
    # Calculate gradients
    dt = 1e-3
    detadt = (eta_t2 - eta_t1)/dt
    detadx = np.diff(eta_t2, axis=0)/(Lx/nx)
    detady = np.diff(eta_t2, axis=1)/(Ly/ny)
    detadx = np.vstack([detadx, -eta_t2[-1,:]])/(Lx/nx)
    detady = np.hstack([detady, -eta_t2[:,-1].reshape(-1,1)])/(Ly/ny)
    
    grad = np.sqrt(detadx**2 + detady**2)
    nhatx = detadx/grad
    nhaty = detady/grad
    
    cx = -detadt*detadx/grad**2
    cy = -detadt*detady/grad**2
    cx = cx/u_star
    cy = cy/u_star
    
    # Calculate delta
    delta = 2.55*amp
    
    # Setup iteration parameters
    max_iterations = 10**3
    tolerance = 1e-6
    initial_guess = 20
    flag_sgs = 0
    zo_u = 10**(-20)
    
    # Calculate U and Lambda
    U = solve_U(detadx, grad, nhatx, nhaty, cx, cy, delta, Ret, zo_u,
                flag_sgs, max_iterations, tolerance, initial_guess)
    lambda_val = 1/U**2
    kappa = 0.4
    
    # Calculate surface roughness
    zo_wave_model = delta*np.exp(-kappa*(1/np.sqrt(lambda_val)))
    zo_wave_model = zo_wave_model*(1/amp)
    
    return zo_wave_model

def calculate_wave_model_multi(data_file: str) -> float:
    """
    Calculate surface roughness for multiscale waves.
    
    Parameters:
        data_file: Path to .npz file containing wave data
        
    Returns:
        zo_wave_model: Normalized surface roughness
    """
    # Load data
    data = np.load(data_file)
    eta_t1 = data['eta_t1']
    eta_t2 = data['eta_t2']
    Ret = data['Ret']
    u_star = data['u_star']
    alpha_p = data['alpha_p']
    kp = data['kp']
    Lx = data['Lx']
    Ly = data['Ly']
    nx = data['nx']
    ny = data['ny']
    Hs_measured = data['Hs_measured']
    
    lambdap = (2*np.pi)/kp
    
    # Calculate gradients
    dt = 0.001
    detadt = (eta_t2 - eta_t1)/dt
    detadx = np.diff(eta_t2, axis=1)/(Lx/nx)
    detady = np.diff(eta_t2, axis=0)/(Ly/ny)
    detadx = np.hstack([detadx, -eta_t2[:,-1].reshape(-1,1)])/(Lx/nx)
    detady = np.vstack([detady, -eta_t2[-1,:]])/(Ly/ny)
    
    grad = np.sqrt(detadx**2 + detady**2)
    nhatx = detadx/grad
    nhaty = detady/grad
    
    cx = -detadt*detadx/grad**2
    cy = -detadt*detady/grad**2
    
    k_min = 0.25*kp
    threshold_wave_speed = np.sqrt(9.81/k_min)
    cx = np.clip(cx, -threshold_wave_speed, threshold_wave_speed)
    cy = np.clip(cy, -threshold_wave_speed, threshold_wave_speed)
    cx = cx/u_star
    cy = cy/u_star
    
    # Calculate delta
    eta = eta_t2/lambdap
    eta_mean = np.mean(eta)
    eta_p = eta - eta_mean
    ramp_eta_p = eta_p * (eta_p >= 0)
    ramp__eta_p8 = ramp_eta_p**8
    H_p = (np.mean(ramp__eta_p8))**(1/8)
    delta = 3*H_p
    
    # Calculate surface roughness for sgs surfaces
    k_filter = np.sqrt((np.pi/(Lx/nx))**2 + (np.pi/(Ly/ny))**2)
    eta_sgs = np.sqrt((0.2*alpha_p)/kp**2*(1 - np.exp(-5/4*(kp/k_filter)**2)))
    kappa = 0.4
    zo_u = eta_sgs*np.exp(-kappa*8.5)/lambdap
    
    # Setup iteration parameters
    max_iterations = 10**3
    tolerance = 1e-6
    initial_guess = 20
    flag_sgs = 1
    
    # Calculate U and Lambda
    U = solve_U(detadx, grad, nhatx, nhaty, cx, cy, delta, Ret, zo_u,
                flag_sgs, max_iterations, tolerance, initial_guess)
    lambda_val = 1/U**2
    
    # Calculate surface roughness
    zo_wave_model = delta*np.exp(-kappa*(1/np.sqrt(lambda_val)))
    zo_wave_model = zo_wave_model*lambdap/Hs_measured
    
    return zo_wave_model