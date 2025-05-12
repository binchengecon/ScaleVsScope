% Parameters
g = 0.02; % Growth rate
N = 100; % Number of grid points
q_min = 0.1;
q_max = 10.0;
q_grid = linspace(q_min, q_max, N)';

% Initial guess for Omega
Omega = ones(N,1) / (q_max - q_min); % uniform initial guess

% Tolerance
tol = 1e-6;
max_iter = 1000;

for iter = 1:max_iter
    % Interpolated Omega at shifted q
    Omega_interp = griddedInterpolant(q_grid, Omega, 'linear', 'nearest');
    
    Omega_shifted = (1/(1+g)) * Omega_interp(q_grid/(1+g));
    
    % Compute error
    err = max(abs(Omega_shifted - Omega));
    
    % Update
    Omega = Omega_shifted;
    
    if err < tol
        fprintf('Converged in %d iterations.\n', iter);
        break;
    end
end

% Normalize Omega to integrate to 1
Omega = Omega / trapz(q_grid, Omega);

% Plot
figure;
plot(q_grid, Omega, 'LineWidth', 2);
xlabel('q');
ylabel('\Omega(q)');
title('Stationary Density \Omega(q) on BGP');
grid on;
