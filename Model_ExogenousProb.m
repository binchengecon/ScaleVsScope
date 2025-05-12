function eqm = Model_ExogenousProb(p, np, func, eqm)



    % eqm.delta_hat(1,:) = trapz(np.q, eqm.delta(1,:) .* eqm.omega(1,:) );
    % eqm.delta_hat(2,:) = trapz(np.q, eqm.delta(2,:) .* eqm.omega(2,:) );

    coef_l = eqm.omega(1,end) /eqm.s(1) / (np.q(end)^(-p.theta - 1));
    coef_h = eqm.omega(2,end) /eqm.s(2) / (np.q(end)^(-p.theta - 1));
    eqm.omega_interp_l = @(q) (q > max(np.q)) .* coef_l .* q.^(-p.theta-1) + ...
        (q >= min(np.q) & q <= max(np.q)) .* max(0, interp1(np.q, eqm.omega(1,:), q, 'linear')); % outside is set  for convergence of calE
    eqm.omega_interp_h = @(q) (q > max(np.q)) .* coef_h .* q.^(-p.theta-1) + ...
        (q >= min(np.q) & q <= max(np.q)) .* max(0, interp1(np.q, eqm.omega(2,:), q, 'linear'));
    eqm.delta_interp_l = @(q) max(0, min(1, interp1(np.q, eqm.delta(1,:), q, 'linear', 'extrap')));
    eqm.delta_interp_h = @(q) max(0, min(1, interp1(np.q, eqm.delta(2,:), q, 'linear', 'extrap')));


    % for i = 1:2
    %     eqm.delta_hat(i) = trapz(np.q, eqm.delta(i,:) .* eqm.omega(i,:) );
    % end
    % integral( @(q) eqm.omega_interp_l(q), np.q_min, np.q_max)

    % eqm.delta_hat(1) = integral( @(q) eqm.delta_interp_l(q) .* eqm.omega_interp_l(q), np.q_min, Inf);
    % eqm.delta_hat(2) = integral( @(q) eqm.delta_interp_h(q) .* eqm.omega_interp_h(q), np.q_min, Inf);


    eqm.delta_hat(1) = integral( @(q) eqm.delta_interp_l(q) .* eqm.omega_interp_l(q), np.q_min, np.q_max);
    eqm.delta_hat(2) = integral( @(q) eqm.delta_interp_h(q) .* eqm.omega_interp_h(q), np.q_min, np.q_max);

    eqm.delta_hat    = max(0, min(1, eqm.delta_hat)); % Ensure delta_hat is not greater than 1

    for i = 1:2
        for iprime = 1:2
            eqm.mu_hat(i,iprime) = eqm.s(i) * eqm.delta_hat(i) * func.P_X(p, p.phi(i)/p.phi(iprime));
        end
    end

    
    flag = eqm.mu_hat(1,1) < 0 || eqm.mu_hat(1,1) > 1 || eqm.mu_hat(2,1) < 0 || eqm.mu_hat(2,1) > 1;
    if flag
        disp('Market share cannot be negative or greater than 1.')

        eqm.mu_hat(1,1) = max(0.01, eqm.mu_hat(1,1));
        eqm.mu_hat(1,1) = min(0.99, eqm.mu_hat(1,1));
        eqm.mu_hat(2,1) = 1 - eqm.mu_hat(1,1);
    end



    
    flag = eqm.mu_hat(1,2) < 0 || eqm.mu_hat(1,2) > 1 || eqm.mu_hat(2,2) < 0 || eqm.mu_hat(2,2) > 1;
    if flag
        disp('Market share cannot be negative or greater than 1.')

        eqm.mu_hat(1,2) = max(0.01, eqm.mu_hat(1,2));
        eqm.mu_hat(1,2) = min(0.99, eqm.mu_hat(1,2));
        eqm.mu_hat(2,2) = 1 - eqm.mu_hat(1,2);
    end



    for i = 1:2
        eqm.nu_hat(i) = 1 - eqm.mu_hat(1,i) - eqm.mu_hat(2,i);
    end



end