

% gamma = [0.02 0.03 0.07];
% gamma = [0.07];
% gamma = [0.10];
gamma = [0.20];
% gamma = [0.10 0.20];

m = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9];
% m = [0.7, 0.75, 0.8, 0.85, 0.9];
% m = [0.75, 0.8, 0.85, 0.9];
% m = [0.85, 0.9];
% m = [0.9];

% phi_h = [1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3];

% r = [0.05 0.1 0.15 0.2 0.25 0.3];
r = [0.15];
% theta = [4 5];
% theta = [3 4];
theta = [5];
% alpha = [0.5 0.8 0.88];
alpha = [0.6 0.7];


for i = 1:length(gamma)
    for j = 1:length(m)
        % for k = 1:length(phi_h)
            for l = 1:length(r)
                for f = 1:length(theta)
                     for g = 1:length(alpha)
                    % p_exogenous.phi_l = 1;
                    % p_exogenous.phi_h = phi_h(k);
                    p_exogenous.r = r(l);
                    % p_exogenous.theta = theta(f);
                    % p_exogenous.phi_l = 1;
                    % p_exogenous.phi_h = 1.05;
                    p_exogenous.alpha = alpha(g);
                    p_exogenous.theta = theta(f);
                    p_exogenous.m_bar = m(j);
                    p_exogenous.gamma = gamma(i);
                    main_func(p_exogenous);
                    end
                 end
            end
        % end
    end
end

