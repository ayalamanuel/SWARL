%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Surface Wave-Aerodynamic Roughness Length Model for Air-Sea Interactions %%%
%%%                                                                          %%%
%%%     - This script does the Newton-Raphson iteration to estimate U        %%%
%%%       which will be used to estimate \Lambda                             %%%
%%%                                                                          %%%
%%%     - Make sure that phase velocity (cx, cy) are normalized by friction  %%%
%%%       velocity                                                           %%%
%%%                                                                          %%%
%%%     - Make sure that delta is normalized by the same length scale used   %%%
%%%       in the turbulent Reynolds number (Ret)                             %%%
%%%                                                                          %%%
%%% Script developed by: Manuel Ayala (mayala5@jhu.edu)                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [U] = solve_U(detadx, grad, nx, ny, c_x, c_y, delta, Ret, zo, flag, max_iterations, tolerance, initial_guess)
    U = initial_guess;
    
    function [f, df] = function_and_derivative(U_current)
        % Calculate Reynolds number
        Re = U_current * delta * Ret;
        
        % Calculate c_fs
        c_fs = 0.0288 * Re^(-1/5) * (1 + 577*Re^(-6/5))^(2/3);
        
        % Calculate c_f with simplified formula
        c_f = 2 * ((1/2*c_fs)^3 + flag*((1/0.4)*log(delta/zo))^(-6))^(1/3);
        
        % Derivative of c_fs with respect to Re
        dc_fs_dRe = 0.0288 * (...
            (-1/5) * Re^(-6/5) * (1 + 577*Re^(-6/5))^(2/3) + ...
            Re^(-1/5) * (2/3) * (1 + 577*Re^(-6/5))^(-1/3) * ...
            577 * (-6/5) * Re^(-11/5)...
        );
        
        % Complete derivative of c_f with respect to Re
        dc_f_dRe = 2 * (1/3) * ((1/2*c_fs)^3 + flag*((1/0.4)*log(delta/zo))^(-6))^(-2/3) * ...
                   (3*(1/2*c_fs)^2 * (1/2)*dc_fs_dRe);
        
        % Convert to derivative with respect to U
        dc_f_dU = dc_f_dRe * (delta * Ret);
        
        % Calculate L
        alpha = atan(grad);
        cx_term = 1 - c_x./U_current;
        term1 = (nx.^2) .* (cx_term.^2);
        term2 = (ny.^2) .* ((c_y./U_current).^2);
        numer = nx .* cx_term;
        dir_term = (numer + abs(numer))./numer;
        L = mean((alpha./(pi + alpha)) .* (term1 + term2) .* detadx .* 0.5 .* dir_term, 'all');
        L = L + 0.5*c_f;

        % Calculate f
        f = U_current - 1/sqrt(L);
        
        % Calculate dL/dU
        dcx_term_dU = c_x./U_current^2;
        dterm1_dU = 2 * (nx.^2) .* cx_term .* dcx_term_dU;
        dterm2_dU = -2 * (ny.^2) .* (c_y.^2./U_current^3);
        dnumer_dU = nx .* dcx_term_dU;
        ddir_dU = (dnumer_dU .* numer - (numer + abs(numer)) .* dnumer_dU) ./ numer.^2;
        dL_avg_dU = mean((alpha./(pi + alpha)) .* (...
            (dterm1_dU + dterm2_dU) .* detadx .* 0.5 .* dir_term + ...
            (term1 + term2) .* detadx .* 0.5 .* ddir_dU ), 'all');
        dL_dU = dL_avg_dU + 0.5*dc_f_dU;

        % Calculate df
        df = 1 + (1/2) * L^(-3/2) * dL_dU;
    end
    
    % Newton-Raphson iteration
    for iter = 1:max_iterations
        U_prev = U;
        [f, df] = function_and_derivative(U);
        U = U - f/df;
        error = abs(U - U_prev) / abs(U);
        
        if error < tolerance
            break;
        end
        if iter == max_iterations
            warning('Maximum iterations reached without convergence');
        end
    end
end