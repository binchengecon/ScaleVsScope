function eqm = Model_Initialization(p, np, func)

    eqm.g = 0.15;
    eqm.gp1_type = [0.1, 0.1]; % g_l = 0.1, g_h = 0.1


    eqm.bellman_constant = (1+eqm.g) * (1-p.gamma) / (1+p.r);

    eqm.calL = zeros(2, 1); % calL function
    eqm.calE = zeros(2, 1); % calE function

    eqm.Q_type_sigma1 = zeros(2, 1); % Q function
    eqm.Q_sigma1 = 1; % Q_prime function

    eqm.nu_hat = zeros(2, 1); % nu function
    eqm.mu_hat = zeros(2, 2); % mu function

    eqm.v = zeros(2, np.q_num); % value function
    eqm.v(1,:) = p.pi(1) * np.q.^(p.sigma-1); % value function
    eqm.v(2,:) = p.pi(2) * np.q.^(p.sigma-1); % value function
    eqm.v_external_expected = 10 * ones(2,2);  % expected value function   [E_{Omega_1} E_X (v_1 (x q/ (1+g))), E_{Omega_2} E_X (v_1 (x q/ (1+g)))
                                  % expected value function    E_{Omega_1} E_X (v_2 (x q/ (1+g))), E_{Omega_2} E_X (v_2 (x q/ (1+g)))     ]
    eqm.v_internal_expected = 10 * ones(2,np.q_num); % internal expected value function
    eqm.v_idle              = 10 * ones(2, np.q_num); % idle value function                 
                                        
    eqm.lambda = .1 * ones(2, np.q_num); % lambda function
    eqm.lambda_hat = .1 * ones(2,1); % lambda function
    eqm.delta = .1 * ones(2, np.q_num); % delta function
    eqm.delta_hat = .1 * ones(2, 1); % delta function

    eqm.s = [.5; .5]; % s function

    eqm.omega = zeros(2, np.q_num); % omega function


    % eqm.omega(1,:) = exp(-np.q); % omega function
    % eqm.omega(2,:) = exp(-np.q); % omega function
    eqm.omega(1,:) =  1./(np.q * p.tau * sqrt(2 * pi)) .* exp(- (log(np.q) - p.iota).^2 / (2 * p.tau^2)); % omega function
    eqm.omega(2,:) = 1./(np.q * p.tau * sqrt(2 * pi)) .* exp(- (log(np.q) - p.iota).^2 / (2 * p.tau^2)); % omega function

    % eqm.omega = eqm.omega ./ sum(eqm.omega * np.dq, 2); % Vectorized normalization
    eqm.omega = eqm.omega ./ trapz(np.q, eqm.omega, 2); % Vectorized normalization
    eqm.omega_tilde = eqm.omega ./ eqm.s; % omega_tild function
    





    file_name = sprintf('epsilon=%.2f, gamma=%.2f, sigma=%d, phi_h=%.2f, m_bar=%.2f, r=%.2f, alpha=%.2f, theta=%.2f, beta=%.2f', p.epsilon, p.gamma, p.sigma, p.phi_h, p.m_bar, p.r, p.alpha, p.theta, p.beta);
    % data = load(['./data_beta/', file_name, '.mat'] , 'p', 'np' , 'eqm_save', 'iter_history');

    try
        data = load(['./data_beta/', file_name, '.mat'] , 'p', 'np' , 'eqm_save', 'iter_history');
    catch ME
        fprintf('[LOAD ERROR] Failed to load solution: %s.mat\nError: %s\n\n', file_name, ME.message);
    end



    eqm = data.eqm_save;




end