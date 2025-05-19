function flag = selector(gamma, m, phi_h, r, theta, alpha, epsilon, beta)


    p_scenario1.gamma = 0.2;
    p_scenario1.m_bar = 0.35;
    p_scenario1.phi_h = 1.05;
    p_scenario1.r = 0.15;
    p_scenario1.theta = 5;
    p_scenario1.alpha = 0.9;
    p_scenario1.epsilon = 0.0;
    p_scenario1.beta = 8;


    p_scenario2.gamma = 0.2;
    p_scenario2.m_bar = 0.35;
    p_scenario2.phi_h = 1.05;
    p_scenario2.r = 0.15;
    p_scenario2.theta = 5;
    p_scenario2.alpha = 0.9;
    p_scenario2.epsilon = 1.0;
    p_scenario2.beta = 8;

    is_scenario1 = ...
        gamma == p_scenario1.gamma && ...
        m == p_scenario1.m_bar && ...
        phi_h == p_scenario1.phi_h && ...
        r == p_scenario1.r && ...
        theta == p_scenario1.theta && ...
        alpha == p_scenario1.alpha && ...
        epsilon == p_scenario1.epsilon && ...
        beta == p_scenario1.beta;
    is_scenario2 = ...
        gamma == p_scenario2.gamma && ...
        m == p_scenario2.m_bar && ...
        phi_h == p_scenario2.phi_h && ...
        r == p_scenario2.r && ...
        theta == p_scenario2.theta && ...
        alpha == p_scenario2.alpha && ...
        epsilon == p_scenario2.epsilon && ...
        beta == p_scenario2.beta;

    if is_scenario1
        flag = 1;
    end

    if is_scenario2
        flag = 2;
    end

    if ~is_scenario1 && ~is_scenario2
        flag = 0;
    end

end