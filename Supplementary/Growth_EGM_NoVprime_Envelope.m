% EGM for Neoclassical Growth Model (based purely on consumption policy)
% -----------------------------------------------------------------------

% Parameters
alpha = 0.3;         % Capital share
beta = 0.95;         % Discount factor
gamma = 2;           % CRRA coefficient

% Grids
nk = 500;
k_min = 0.01; k_max = 5;
k_prime_grid = linspace(k_min, k_max, nk)'; % Grid for next-period capital k'
k_grid = linspace(k_min, k_max, nk)';       % Fixed grid for current capital k

% Initial guess for consumption policy (consume all output)
c_policy_old = k_grid.^alpha;

% Tolerance and iteration setup
tol = 1e-6;
max_iter = 1000;

% Storage for policies
c_policy = zeros(nk, 1);
k_policy = zeros(nk, 1);

for iter = 1:max_iter
    % Step 1: Interpolate c'(k') from previous iteration
    c_prime = interp1(k_grid, c_policy_old, k_prime_grid, 'linear', 'extrap');

    % Step 2: Apply Euler equation with envelope condition
    uc_prime = c_prime.^(-gamma);
    uc_today = beta * uc_prime .* alpha .* k_prime_grid.^(alpha-1);
    c_today = uc_today.^(-1/gamma);

    % Step 3: Recover current k using budget constraint
    k_now = (c_today + k_prime_grid).^(1/alpha);

    % Step 4: Interpolate policies onto fixed k grid
    valid = (k_now >= k_min) & (k_now <= k_max);
    k_now = k_now(valid);
    k_prime_valid = k_prime_grid(valid);
    c_valid = c_today(valid);

    k_policy = interp1(k_now, k_prime_valid, k_grid, 'linear', 'extrap');
    c_policy = interp1(k_now, c_valid, k_grid, 'linear', 'extrap');

    % Step 5: Check convergence
    if max(abs(c_policy - c_policy_old)) < tol
        break;
    end

    c_policy_old = c_policy;
end

fprintf('EGM converged in %d iterations.\n', iter);

% Plot Policy Functions
figure;
% subplot(2,1,1);
plot(k_grid, k_policy, 'LineWidth', 1.5);
xlabel('k today'); ylabel('k tomorrow');
title('Capital Policy Function k''(k)'); grid on;

% subplot(2,1,2);
% plot(k_grid, c_policy, 'LineWidth', 1.5);
% xlabel('k today'); ylabel('Consumption');
% title('Consumption Policy Function c(k)'); grid on;


