function eqm = Model_AggregateGrowth_Discount(p, np, func, eqm)

    % eqm.omega_interp_l = @(q) max(0,  interp1(np.q, eqm.omega(1,:), q, 'spline')); % 0 to cap negative  numbers
    % eqm.omega_interp_h = @(q) max(0,  interp1(np.q, eqm.omega(2,:), q, 'spline')); % 0 to cap negative  numbers

    % coef_l = eqm.omega(1,end)/eqm.s(1) / (np.q(end)^(-p.sigma - 1));
    % coef_h = eqm.omega(2,end)/eqm.s(2) / (np.q(end)^(-p.sigma - 1));
    % eqm.omega_interp_l = @(q) (q > max(np.q)) .* coef_l .* q.^(-p.sigma-1) + ...
    %     (q >= min(np.q) & q <= max(np.q)) .* max(0, interp1(np.q, eqm.omega(1,:), q, 'spline')); % outside is set  for convergence of calE
    % eqm.omega_interp_h = @(q) (q > max(np.q)) .* coef_h .* q.^(-p.sigma-1) + ...
    %     (q >= min(np.q) & q <= max(np.q)) .* max(0, interp1(np.q, eqm.omega(2,:), q, 'spline'));


    coef_l = eqm.omega(1,end) / (np.q(end)^(-p.theta - 1));
    coef_h = eqm.omega(2,end) / (np.q(end)^(-p.theta - 1));
    eqm.omega_tilde_interp_l = @(q) (q > max(np.q)) .* coef_l .* q.^(-p.theta-1) + ...
        (q >= min(np.q) & q <= max(np.q)) .* max(0, interp1(np.q, eqm.omega_tilde(1,:), q, 'spline', 'extrap')); % outside is set  for convergence of calE
    eqm.omega_tilde_interp_h = @(q) (q > max(np.q)) .* coef_h .* q.^(-p.theta-1) + ...
        (q >= min(np.q) & q <= max(np.q)) .* max(0, interp1(np.q, eqm.omega_tilde(2,:), q, 'spline', 'extrap'));
    eqm.lambda_interp_l = @(q) max(0, min(1, interp1(np.q, eqm.lambda(1,:), q, 'spline', 'extrap')));
    eqm.lambda_interp_h = @(q) max(0, min(1, interp1(np.q, eqm.lambda(2,:), q, 'spline', 'extrap')));

    % eqm.calE(1)        = integral(  @(q)  q.^(p.sigma - 1) .* eqm.omega_tilde_interp_l(q) , np.q_min, Inf);  % need more than 100 like 200 poitns to achieve accuracy
    % eqm.calE(2)        = integral(  @(q)  q.^(p.sigma - 1) .* eqm.omega_tilde_interp_h(q) , np.q_min, Inf);

    % eqm.calL(1) = integral(  @(q)  q.^(p.sigma - 1) .* eqm.lambda_interp_l(q) .* eqm.omega_tilde_interp_l(q) , np.q_min, Inf);
    % eqm.calL(2) = integral(  @(q)  q.^(p.sigma - 1) .* eqm.lambda_interp_h(q) .* eqm.omega_tilde_interp_h(q) , np.q_min, Inf);



    eqm.calE(1)        = integral(  @(q)  q.^(p.sigma - 1) .* eqm.omega_tilde_interp_l(q) , np.q_min, 3*  np.q_max);  % need more than 100 like 200 poitns to achieve accuracy
    eqm.calE(2)        = integral(  @(q)  q.^(p.sigma - 1) .* eqm.omega_tilde_interp_h(q) , np.q_min, 3*  np.q_max);  % need

    eqm.calL(1) = integral(  @(q)  q.^(p.sigma - 1) .* eqm.lambda_interp_l(q) .* eqm.omega_tilde_interp_l(q) , np.q_min, np.q_max);  % need
    eqm.calL(2) = integral(  @(q)  q.^(p.sigma - 1) .* eqm.lambda_interp_h(q) .* eqm.omega_tilde_interp_h(q) , np.q_min, np.q_max);  % need


    % eqm.calE(1) = trapz(np.q, np.q.^(p.sigma - 1) .* eqm.omega(1,:));
    % eqm.calE(2) = trapz(np.q, np.q.^(p.sigma - 1) .* eqm.omega(2,:));

    % eqm.calL(1) = trapz(np.q, np.q.^(p.sigma - 1) .* eqm.lambda(1,:) .*  eqm.omega_tilde(1,:));
    % eqm.calL(2) = trapz(np.q, np.q.^(p.sigma - 1) .* eqm.lambda(2,:) .*  eqm.omega_tilde(2,:));

    eqm.Q_type_sigma1(1) = eqm.calE(1) / ( eqm.s(1) * p.phi(1)^(p.sigma - 1)  );
    eqm.Q_type_sigma1(2) = eqm.calE(2) / ( eqm.s(2) * p.phi(2)^(p.sigma - 1)  );
    eqm.Q_sigma1 = eqm.Q_type_sigma1(1) + eqm.Q_type_sigma1(2);


    % eqm.gp1_type(1) = (1 - p.gamma) * eqm.nu_hat(1) * eqm.calL(1) / eqm.calE(1) * (p.theta / (p.theta+1-p.sigma)) ...
    %                 + (1 - p.gamma) * eqm.nu_hat(1) * ( 1 - eqm.calL(1) / eqm.calE(1) ) ...
    %                 + (1 - p.gamma) * eqm.mu_hat(1,1) * ( p.theta * p.alpha^(p.sigma - 1) / (p.theta + 1 - p.sigma)   )...
    %                 + (1 - p.gamma) * eqm.mu_hat(1,2) *  eqm.calE(2) / eqm.calE(1) *  (  p.theta * p.alpha^(p.sigma - 1) / (p.theta + 1 - p.sigma)   )...
    %                 + p.gamma * p.m(1) * ( 1 + eqm.calE(2) / eqm.calE(1) ) ;
    % eqm.gp1_type(2) = (1 - p.gamma) * eqm.nu_hat(2) * eqm.calL(2) / eqm.calE(2) * (p.theta / (p.theta+1-p.sigma)) ...
    %                 + (1 - p.gamma) * eqm.nu_hat(2) * ( 1 - eqm.calL(2) / eqm.calE(2) ) ...
    %                 + (1 - p.gamma) * eqm.mu_hat(2,2) * ( p.theta * p.alpha^(p.sigma - 1) / (p.theta + 1 - p.sigma)   )...
    %                 + (1 - p.gamma) * eqm.mu_hat(2,1) *  eqm.calE(1) / eqm.calE(2) *  (  p.theta * p.alpha^(p.sigma - 1) / (p.theta + 1 - p.sigma)   )...
    %                 + p.gamma * p.m(2) * ( 1 + eqm.calE(1) / eqm.calE(2) ) ;

    eqm.gp1_type(1) = (1 - p.gamma) * eqm.nu_hat(1) * eqm.calL(1) / eqm.calE(1) * (p.theta / (p.theta+1-p.sigma)) ...
                    + (1 - p.gamma) * eqm.nu_hat(1) * ( 1 - eqm.calL(1) / eqm.calE(1) ) ...
                    + (1 - p.gamma) * eqm.mu_hat(1,1) * ( p.theta * p.alpha^(p.sigma - 1) / (p.theta + 1 - p.sigma)   )...
                    + (1 - p.gamma) * eqm.mu_hat(1,2) *  eqm.calE(2) / eqm.calE(1) *  (  p.theta * p.alpha^(p.sigma - 1) / (p.theta + 1 - p.sigma)   )...
                    + p.gamma * p.m(1) * (1-p.zeta)^(p.sigma-1) * ( 1 + eqm.calE(2) / eqm.calE(1) ) ;
    eqm.gp1_type(2) = (1 - p.gamma) * eqm.nu_hat(2) * eqm.calL(2) / eqm.calE(2) * (p.theta / (p.theta+1-p.sigma)) ...
                    + (1 - p.gamma) * eqm.nu_hat(2) * ( 1 - eqm.calL(2) / eqm.calE(2) ) ...
                    + (1 - p.gamma) * eqm.mu_hat(2,2) * ( p.theta * p.alpha^(p.sigma - 1) / (p.theta + 1 - p.sigma)   )...
                    + (1 - p.gamma) * eqm.mu_hat(2,1) *  eqm.calE(1) / eqm.calE(2) *  (  p.theta * p.alpha^(p.sigma - 1) / (p.theta + 1 - p.sigma)   )...
                    + p.gamma * p.m(2) * (1-p.zeta)^(p.sigma-1) * ( 1 + eqm.calE(1) / eqm.calE(2) ) ;

    eqm.g  = eqm.Q_type_sigma1(1) * eqm.gp1_type(1) + eqm.Q_type_sigma1(2) * eqm.gp1_type(2);
    eqm.g  = eqm.g / eqm.Q_sigma1; % Normalize g to be between 0 and 1\
    % eqm.g  = max(1, eqm.g); % Ensure g is not negative
    eqm.g  = eqm.g^(1/ (p.sigma - 1)) - 1; % RHS^(1/(sigma-1))-1 to get the growth rate

    eqm.g_share = eqm.Q_type_sigma1 /eqm.Q_sigma1;


    
end