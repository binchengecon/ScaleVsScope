function eqm = Model_LoM_PDF_Discount(p, np, func, eqm)

    eqm.lambda_omega_tilde_fY_gp1_interp_l = @(q) integral(   @(y) 1./y .* eqm.lambda_interp_l(q .*(1+eqm.g) ./y) .*  eqm.omega_tilde_interp_l(q .*(1+eqm.g) ./y) .* p.theta./ y.^(p.theta+1), 1, Inf);
    eqm.lambda_omega_tilde_fY_gp1_interp_h = @(q) integral(   @(y) 1./y .* eqm.lambda_interp_h(q .*(1+eqm.g) ./y) .*  eqm.omega_tilde_interp_h(q .*(1+eqm.g) ./y) .* p.theta./ y.^(p.theta+1), 1, Inf);
    
    eqm.omega_tilde_fX_gp1_interp_l = @(q) integral(   @(x) 1./x .* eqm.omega_tilde_interp_l(q .*(1+eqm.g) ./x) .* p.theta * p.alpha^(p.theta) ./ x.^(p.theta+1), p.alpha, Inf);
    eqm.omega_tilde_fX_gp1_interp_h = @(q) integral(   @(x) 1./x .* eqm.omega_tilde_interp_h(q .*(1+eqm.g) ./x) .* p.theta * p.alpha^(p.theta) ./ x.^(p.theta+1), p.alpha, Inf);
    
    eqm.omega_tilde_gp1_interp_l = @(q) eqm.omega_tilde_interp_l(q .*(1+eqm.g));
    eqm.omega_tilde_gp1_interp_h = @(q) eqm.omega_tilde_interp_h(q .*(1+eqm.g));



    coef_l = eqm.omega(1,end) / (np.q(end)^(-p.theta - 1));
    coef_h = eqm.omega(2,end) / (np.q(end)^(-p.theta - 1));
    eqm.omega_tilde_gp1_interp_l = @(q) (q .*(1+eqm.g) >  max(np.q)) .* coef_l .* (q .*(1+eqm.g)).^(-p.theta-1) + ...
        (q .*(1+eqm.g) >= min(np.q) & q .*(1+eqm.g) <=  max(np.q)) .* max(0, interp1(np.q, eqm.omega_tilde(1,:), q .*(1+eqm.g), 'spline', 'extrap')); % outside is set  for convergence of calE
    eqm.omega_tilde_gp1_interp_h = @(q) (q .*(1+eqm.g) > max(np.q)) .* coef_h .* (q .*(1+eqm.g)).^(-p.theta-1) + ...
        (q .*(1+eqm.g) >= min(np.q) & q .*(1+eqm.g) <=  max(np.q)) .* max(0, interp1(np.q, eqm.omega_tilde(2,:), q .*(1+eqm.g), 'spline', 'extrap'));
    eqm.lambda_gp1_interp_l = @(q) max(0, min(1, interp1(np.q, eqm.lambda(1,:), q .*(1+eqm.g), 'linear', 'extrap')));
    eqm.lambda_gp1_interp_h = @(q) max(0, min(1, interp1(np.q, eqm.lambda(2,:), q .*(1+eqm.g), 'linear', 'extrap')));

    % eqm.omega_tilde_exo_exponential_l  = 1/p.exp *  exp( -1/p.exp * np.q * (1+eqm.g));
    % eqm.omega_tilde_exo_exponential_h  =  1/p.exp *  exp( -1/p.exp * np.q * (1+eqm.g));

    % omega_entry =    1./(np.q * (1+eqm.g) * p.tau * sqrt(2 * pi)) .* exp(- (log(np.q * (1+eqm.g) ) - p.iota).^2 / (2 * p.tau^2));


    % eqm.omega_tilde_exo_exponential_l  =  omega_entry;
    % eqm.omega_tilde_exo_exponential_h  =  omega_entry;


    eqm.omega_tilde_gp1_zeta_interp_l = @(q) eqm.omega_tilde_interp_l(q .*(1+eqm.g)/(1 - p.zeta));
    eqm.omega_tilde_gp1_zeta_interp_h = @(q) eqm.omega_tilde_interp_h(q .*(1+eqm.g)/(1 - p.zeta));




    LoM_RHS = zeros(2, np.q_num);
    LoM_RHS(1,:) = (1 - p.gamma) * eqm.nu_hat(1) * arrayfun(eqm.lambda_omega_tilde_fY_gp1_interp_l, np.q) ...
                + (1 - p.gamma) * eqm.nu_hat(1) * ( 1 - arrayfun(eqm.lambda_gp1_interp_l, np.q)) .* arrayfun(eqm.omega_tilde_gp1_interp_l, np.q) ...
                + (1 - p.gamma) * eqm.mu_hat(1,1) * arrayfun(eqm.omega_tilde_fX_gp1_interp_l, np.q) ...
                + (1 - p.gamma) * eqm.mu_hat(1,2) * arrayfun(eqm.omega_tilde_fX_gp1_interp_h, np.q) ...
                + p.gamma * p.m(1) * ( arrayfun(eqm.omega_tilde_gp1_zeta_interp_l, np.q) + arrayfun(eqm.omega_tilde_gp1_zeta_interp_h, np.q) )/(1-p.zeta);

    LoM_RHS(2,:) = (1 - p.gamma) * eqm.nu_hat(2) * arrayfun(eqm.lambda_omega_tilde_fY_gp1_interp_h, np.q) ...
                + (1 - p.gamma) * eqm.nu_hat(2) * ( 1 - arrayfun(eqm.lambda_gp1_interp_h, np.q) ) .* arrayfun(eqm.omega_tilde_gp1_interp_h, np.q) ...
                + (1 - p.gamma) * eqm.mu_hat(2,2) * arrayfun(eqm.omega_tilde_fX_gp1_interp_h, np.q) ...
                + (1 - p.gamma) * eqm.mu_hat(2,1) * arrayfun(eqm.omega_tilde_fX_gp1_interp_l, np.q) ...
                + p.gamma * p.m(2) * (arrayfun(eqm.omega_tilde_gp1_zeta_interp_l, np.q) + arrayfun(eqm.omega_tilde_gp1_zeta_interp_h, np.q) )/(1-p.zeta);

    % LoM_RHS = zeros(2, np.q_num);
    % LoM_RHS(1,:) = (1 - p.gamma) * eqm.nu_hat(1) * arrayfun(eqm.lambda_omega_tilde_fY_gp1_interp_l, np.q) ...
    %             + (1 - p.gamma) * eqm.nu_hat(1) * ( 1 - arrayfun(eqm.lambda_gp1_interp_l, np.q)) .* arrayfun(eqm.omega_tilde_gp1_interp_l, np.q) ...
    %             + (1 - p.gamma) * eqm.mu_hat(1,1) * arrayfun(eqm.omega_tilde_fX_gp1_interp_l, np.q) ...
    %             + (1 - p.gamma) * eqm.mu_hat(1,2) * arrayfun(eqm.omega_tilde_fX_gp1_interp_h, np.q) ...
    %             + p.gamma * p.m(1) * eqm.omega_tilde_exo_exponential_l;

    % LoM_RHS(2,:) = (1 - p.gamma) * eqm.nu_hat(2) * arrayfun(eqm.lambda_omega_tilde_fY_gp1_interp_h, np.q) ...
    %             + (1 - p.gamma) * eqm.nu_hat(2) * ( 1 - arrayfun(eqm.lambda_gp1_interp_h, np.q) ) .* arrayfun(eqm.omega_tilde_gp1_interp_h, np.q) ...
    %             + (1 - p.gamma) * eqm.mu_hat(2,2) * arrayfun(eqm.omega_tilde_fX_gp1_interp_h, np.q) ...
    %             + (1 - p.gamma) * eqm.mu_hat(2,1) * arrayfun(eqm.omega_tilde_fX_gp1_interp_l, np.q) ...
    %             + p.gamma * p.m(2) * eqm.omega_tilde_exo_exponential_h;
    
    LoM_RHS      = (1 +eqm.g) * LoM_RHS;

    % coef_l = LoM_RHS(1,end) / (np.q(end)^(-p.theta - 1));
    % coef_h = LoM_RHS(2,end) / (np.q(end)^(-p.theta - 1));
    % eqm.omega_tilde_interp_l = @(q) (q > max(np.q)) .* coef_l .* q.^(-p.theta-1) + ...
    %     (q >= min(np.q) & q <= max(np.q)) .* max(0, interp1(np.q, LoM_RHS(1,:), q * (1 +eqm.g), 'spline')); % outside is set  for convergence of calE
    % eqm.omega_tilde_interp_h = @(q) (q > max(np.q)) .* coef_h .* q.^(-p.theta-1) + ...
    %     (q >= min(np.q) & q <= max(np.q)) .* max(0, interp1(np.q, LoM_RHS(2,:), q * (1 +eqm.g), 'spline'));

    % eqm.omega_tilde(1,:) =arrayfun(eqm.omega_tilde_interp_l, np.q);
    % eqm.omega_tilde(2,:) =arrayfun(eqm.omega_tilde_interp_h, np.q);

    
    eqm.omega_tilde = LoM_RHS;

    % rescale omega to ensure it is a PDF
    eqm.omega            = eqm.omega_tilde ./ eqm.s;
    eqm.omega            = eqm.omega ./ trapz(np.q, eqm.omega, 2); % Vectorized normalization
    eqm.omega_tilde      = eqm.omega .* eqm.s; % Vectorized normalization

    % eqm.omega_tilde(1,:) = interp1(np.q, LoM_RHS(1,:), np.q * (1 +eqm.g), 'spline', 'extrap');
    % eqm.omega_tilde(2,:) = interp1(np.q, LoM_RHS(2,:), np.q * (1 +eqm.g), 'spline', 'extrap');

    % eqm.omega_tilde      = max(0, eqm.omega_tilde);
    % eqm.omega            = eqm.omega_tilde ./ eqm.s;
    % % eqm.omega            = eqm.omega ./ sum(eqm.omega * np.dq, 2); % Vectorized normalization
    % eqm.omega            = eqm.omega ./ trapz(np.q, eqm.omega, 2); % Vectorized normalization
    % eqm.omega_tilde      = eqm.omega .* eqm.s; % Vectorized normalization
    % eqm.Omega(1,:) = cumtrapz(np.q, eqm.omega(1,:));
    % eqm.Omega(2,:) = cumtrapz(np.q, eqm.omega(2,:));
    


end