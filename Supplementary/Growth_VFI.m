% Parameters
alpha = 0.3; beta = 0.95; gamma = 2;

% Grids
nk = 1000;
k_min = 0.01; k_max = 5;
k_grid = linspace(k_min, k_max, nk)';

% Initialization
V = zeros(nk,1);
policy_k = zeros(nk,1);
policy_c = zeros(nk,1);
tol = 1e-6;
max_iter = 1000;

for iter = 1:max_iter
    V_new = zeros(nk,1);
    for i = 1:nk
        k = k_grid(i);
        y = k^alpha;
        c = y - k_grid;
        c(c <= 0) = NaN;
        u = (c.^(1-gamma)) / (1-gamma);
        RHS = u + beta * V;
        [V_new(i), idx] = max(RHS);
        policy_k(i) = k_grid(idx);
        policy_c(i) = c(idx);
    end
    if max(abs(V_new - V)) < tol
        break;
    end
    V = V_new;
end

% Plot
figure;
plot(k_grid, policy_k);
xlabel('k today'); ylabel('k tomorrow');
title('Policy Function for Capital');

