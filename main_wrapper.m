function main_wrapper()
    % Get input arguments
    args = getenv('MATLAB_ARGS');
    tokens = str2double(strsplit(args));

    gamma   = tokens(1);
    m_bar   = tokens(2);
    r       = tokens(3);
    theta   = tokens(4);
    alpha   = tokens(5);
    phi_h   = tokens(6);
    beta   = tokens(7);
    epsilon = tokens(8);

    p_exogenous.r = r;
    p_exogenous.theta = theta;
    p_exogenous.alpha = alpha;
    p_exogenous.m_bar = m_bar;
    p_exogenous.gamma = gamma;
    p_exogenous.phi_h = phi_h;
    p_exogenous.beta = beta;
    p_exogenous.epsilon = epsilon;


    main_func(p_exogenous);  % Your main model
end
