% Parameters
beta = 0.96;
gamma = 2; % CRRA coefficient
R = 1.03; % gross interest rate
b = 0.2;    % borrowing limit: a >= -b
ygrid = [.1, 2];          % income states
P = [0.9, 0.1; 0.1, 0.9];% income transition matrix
ny = length(ygrid);

% Grids
na = 100; 
a_min = -b; 
a_max = 20;
a_prime_grid = linspace(a_min, a_max, na)'; % grid for a'

% Initial guess for consumption policy c0(a', y)
c_guess = zeros(na, ny);
for iy = 1:ny
    c_guess(:, iy) = R * a_prime_grid + ygrid(iy); % initial guess: consume all resources
end

% Tolerances
tol = 1e-6;
max_iter = 1000;

for iter = 1:max_iter
    c_new = zeros(na, ny);
    a_grid_star = zeros(na, ny); % endogenous current asset grid

    for iy = 1:ny
        % Step 3: Construct RHS of Euler equation
        B = zeros(na, 1);
        for jy = 1:ny
            % uc(c) = c^(-gamma)
            mu_c = c_guess(:, jy).^(-gamma);
            B = B + P(iy, jy) * mu_c;
        end
        B = beta * R * B; % RHS of Euler equation

        % Step 4: Invert Euler equation to get consumption today
        c_tilde = B.^(-1/gamma);

        % Step 5: From budget constraint, solve for a (today's asset)
        a_today = (c_tilde + a_prime_grid - ygrid(iy)) / R;

        % Save
        a_grid_star(:, iy) = a_today;
        c_new(:, iy) = c_tilde;
    end

    % Step 6: Interpolate onto original a grid
    c_interp = zeros(na, ny);
    for iy = 1:ny
        a_star = a_grid_star(:, iy);
        c_tilde = c_new(:, iy);

        % Find the smallest a_today where borrowing constraint binds (a' = -b)
        idx_borrow = find(a_star <= a_min, 1, 'first');
        if isempty(idx_borrow)
            idx_borrow = 1;
        end

        for ia = 1:na
            if a_grid(ia) < a_star(idx_borrow)
                % Borrowing constraint binds: use budget constraint
                c_interp(ia, iy) = R * a_grid(ia) + ygrid(iy) - a_min;
            else
                % Interpolate c at a_grid
                c_interp(ia, iy) = interp1(a_star, c_tilde, a_grid(ia), 'linear', 'extrap');
            end
        end
    end

    % Step 7: Check convergence
    diff = max(abs(c_interp - c_guess), [], 'all');
    if diff < tol
        fprintf('Converged in %d iterations.\n', iter);
        break
    end

    % Update guess
    c_guess = c_interp;
end

% Plot policy function
figure;
plot(a_grid, c_interp)
legend('Low income', 'High income')
xlabel('Assets')
ylabel('Consumption')
title('EGM Consumption Policy Function')
