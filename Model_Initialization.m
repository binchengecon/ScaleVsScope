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
    

    % eqm.Omega = zeros(2, np.q_num); % Omega function
    % eqm.Omega(1,:) = cumtrapz(np.q, eqm.omega(1,:));
    % eqm.Omega(2,:) = cumtrapz(np.q, eqm.omega(2,:));
    % eqm.Omega = eqm.Omega ./ max(eqm.Omega); % Vectorized normalization
    
    % eqm.Omega(1,:) = linspace(0,1, np.q_num); % Omega function
    % eqm.Omega(2,:) = linspace(0,1, np.q_num); % Omega function
    % eqm.Omega_tilde = zeros(2, np.q_num); % Omega_tild function
    % eqm.Omega_tilde = eqm.Omega ./ eqm.s; % Vectorized Omega_tilde calculation
    % eqm.omega(:,1:end-1) =  (eqm.Omega(:,2:end) - eqm.Omega(:,1:end-1)) ; % omega function
    % eqm.omega(:,end) = (eqm.Omega(:,end) - eqm.Omega(:,end-1)) ; % omega function
    % eqm.omega = eqm.omega / np.dq; % omega function
    % eqm.omega(1, :) = eqm.omega(1,:) / sum(eqm.omega(1,:) * np.dq); % omega function
    % eqm.omega(2, :) = eqm.omega(2,:) / sum(eqm.omega(2,:) * np.dq); % omega function
    % eqm.omega = eqm.omega ./ sum(eqm.omega * np.dq, 2); % Vectorized normalization

    % eqm.omega_tilde = zeros(2, np.q_num); % omega_tild function
    

end