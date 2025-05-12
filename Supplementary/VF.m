% Parameters
beta = 0.95;
g = 0.02; % Growth rate
N = 100; % Number of grid points
q_min = 0.1;
q_max = 10.0;
q_grid = linspace(q_min, q_max, N)'; % Column vector

% Utility function
utility = @(a) log(a);

% Initialize value function
V = zeros(N,1);
V_new = zeros(N,1);

% Tolerance and maximum iterations
tol = 1e-6;
max_iter = 1000;

for iter = 1:max_iter
    % Interpolate V using linear interpolation
    V_interp = griddedInterpolant(q_grid, V, 'linear', 'linear');
    
    for i = 1:N
        q = q_grid(i);
        a = q; % Optimal a is q
        u = utility(a);
        q_next = q / (g + 1);
        
        % Evaluate interpolated value at q_next
        V_future = V_interp(q_next);
        
        % Update V_new
        V_new(i) = u + beta * V_future;
    end
    
    % Convergence check
    if max(abs(V_new - V)) < tol
        fprintf('Converged in %d iterations.\n', iter);
        break;
    end
    
    % Update
    V = V_new;
end

% Plot the value function
figure;
plot(q_grid, V, 'LineWidth', 2);
xlabel('q');
ylabel('Value Function V(q)');
title('Value Function on Balanced Growth Path');
grid on;


% Explicit solution
B_exact = 1/(1 - beta);
A_exact = (-beta / (1 - beta)^2) * log(g + 1);
V_exact = A_exact + B_exact * log(q_grid);

% Plot comparison
figure;
plot(q_grid, V, 'b-', 'LineWidth', 2); hold on;
plot(q_grid, V_exact, 'r--', 'LineWidth', 2);
legend('Numerical V(q)', 'Explicit V(q)');
xlabel('q');
ylabel('Value Function V(q)');
title('Comparison of Numerical vs Explicit Value Function');
grid on;