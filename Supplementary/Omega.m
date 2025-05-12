% Parameters
g = 0.02; % growth rate
N = 100;
q_min = 0.1;
q_max = 10.0;
q_grid = linspace(q_min, q_max, N)'; % Column vector

% Initialize Omega
Omega = ones(N,1); % guess uniform density
Omega_new = zeros(N,1);

% Tolerance
tol = 1e-6;
max_iter = 1000;

for iter = 1:max_iter
    % Interpolation
    Omega_interp = griddedInterpolant(q_grid, Omega, 'linear', 'linear');
    
    for i = 1:N
        q = q_grid(i);
        q_shifted = q / (1+g);
        Omega_shifted = Omega_interp(q_shifted);
        
        % Example RHS: suppose simple decay (placeholder)
        RHS = Omega(i) * exp(-q); % you can change this
        
        % Update according to equation
        Omega_new(i) = (1+g) * RHS;
    end
    
    if max(abs(Omega_new - Omega)) < tol
        fprintf('Converged in %d iterations.\n', iter);
        break;
    end
    
    Omega = Omega_new;
end

% Normalize if needed (integral to 1)
Omega = Omega / trapz(q_grid, Omega);

% Plot
figure;
plot(q_grid, Omega, 'LineWidth', 2);
xlabel('q');
ylabel('\Omega(q)');
title('Stationary Density on BGP');
grid on;
