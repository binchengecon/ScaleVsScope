% EGM for Neoclassical Growth Model
% ----------------------------------

% Parameters
alpha = 0.3;         % Capital share
beta = 0.95;         % Discount factor
gamma = 2;           % CRRA coefficient

% Grids
nk = 500;
k_min = 0.01; k_max = 5;
k_prime_grid = linspace(k_min, k_max, nk)'; % Grid for k' (endogenous grid)

% Initial guess for marginal value function V'(k')
V_prime = ones(nk, 1);  % Constant guess to start (can be refined later)

% Tolerance and iteration setup
tol = 1e-6;
max_iter = 1000;

% Fixed k grid for interpolation of policy function
k_grid = linspace(k_min, k_max, nk)';
k_policy = zeros(nk, 1);  % policy function: k'(k)
c_policy = zeros(nk, 1);  % consumption policy: c(k)

for iter = 1:max_iter
    % Step 1: Invert Euler Equation to get c(k')
    c_tilde = (beta * V_prime).^- (1 / gamma);  % u'(c) = beta * V'(k') -> c = (...)

    % Step 2: Recover implied current k using budget constraint
    k_now = (c_tilde + k_prime_grid).^(1 / alpha);  % since c + k' = y = k^alpha

    % Step 3: Interpolate onto fixed k grid
    valid = k_now >= k_min & k_now <= k_max;
    k_now = k_now(valid);
    k_prime_valid = k_prime_grid(valid);
    c_valid = c_tilde(valid);

    % Interpolate policy functions onto fixed k grid
    k_policy = interp1(k_now, k_prime_valid, k_grid, 'linear', 'extrap');
    c_policy = interp1(k_now, c_valid, k_grid, 'linear', 'extrap');

    % Step 4: Update V'(k) using envelope condition
    V_prime_new = (c_policy).^(-gamma) .* alpha .* (k_grid).^(alpha - 1);

    % Step 5: Check convergence
    if max(abs(V_prime_new - V_prime)) < tol
        break;
    end
    V_prime = V_prime_new;
end

fprintf('EGM converged in %d iterations.\n', iter);

% Plot Policy Function
figure;
plot(k_grid, k_policy, 'LineWidth', 1.5);
xlabel('k today'); ylabel('k tomorrow');
title('Capital Policy Function via EGM'); grid on;


