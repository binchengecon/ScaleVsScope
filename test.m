clear; close all; clc


%%

p = EconomicsParameters();


f = Utility;


np = NumericalParameters(f,p);

eqm.delta_h = 0.1*ones(size(np.q));
eqm.delta_l = 0.1*ones(size(np.q));

eqm.lambda_h = 0.1*linspace(0,1, np.q_num);
eqm.lambda_l = 0.1*linspace(0,1, np.q_num);

eqm.Omega_h = linspace(0,1, np.q_num)';
eqm.Omega_l = linspace(0,1, np.q_num)';
eqm.Omega_h = eqm.Omega_h.^.5;
eqm.Omega_l = eqm.Omega_l.^.5;



eqm.v_h_init = p.Pi_Vec(2) * np.q.^(p.sigma - 1)  + eqm.delta_h.^(1+p.beta)/(1+p.beta);
eqm.v_l_init = p.Pi_Vec(1) * np.q.^(p.sigma - 1)  + eqm.delta_l.^(1+p.beta)/(1+p.beta);


%%

for iter = 1:500

    eqm.v_h = eqm.v_h_init;
    eqm.v_l = eqm.v_l_init;

    
    eqm.delta_h_hat = sum(eqm.delta_h.*f.forward_diff_Omega(eqm.Omega_h));
    eqm.delta_l_hat = sum(eqm.delta_l.*f.forward_diff_Omega(eqm.Omega_l));
    eqm.lambda_h_hat = sum(eqm.lambda_h.*f.forward_diff_Omega(eqm.Omega_h));
    eqm.lambda_l_hat = sum(eqm.lambda_l.*f.forward_diff_Omega(eqm.Omega_l));


    eqm.mu0_hh = eqm.delta_h_hat * f.P_X(p, p.Phi_Matrix(2,2));
    eqm.mu0_hl = eqm.delta_h_hat * f.P_X(p, p.Phi_Matrix(2,1));
    eqm.mu0_lh = eqm.delta_l_hat * f.P_X(p, p.Phi_Matrix(1,2));
    eqm.mu0_ll = eqm.delta_l_hat * f.P_X(p, p.Phi_Matrix(1,1));

    eqm.nu0_h = 1 - eqm.mu0_hh - eqm.mu0_lh;
    eqm.nu0_l = 1 - eqm.mu0_hl - eqm.mu0_ll;


    eqm.s_l = (p.gamma * p.m_bar - (1-p.gamma) * eqm.mu0_lh) / (p.gamma + (1-p.gamma) * eqm.mu0_hl - (1-p.gamma) * eqm.mu0_lh);
    eqm.s_l = real(eqm.s_l);
    eqm.s_l = min(max(eqm.s_l, 0), 0.95);
    eqm.s_h = 1 - eqm.s_l;
    eqm.s_h = real(eqm.s_h);
    eqm.s_h = min(max(eqm.s_h, 0), 1);

    eqm.mu_hh = eqm.s_h * eqm.mu0_hh;
    eqm.mu_hl = eqm.s_l * eqm.mu0_hl;
    eqm.mu_lh = eqm.s_h * eqm.mu0_lh;
    eqm.mu_ll = eqm.s_l * eqm.mu0_ll;

    eqm.nu_h = eqm.s_h * eqm.nu0_h;
    eqm.nu_l = eqm.s_l * eqm.nu0_l;


    %%



    temp_h1 =  (1-p.gamma) * (1-eqm.lambda_h) * eqm.nu0_h .* eqm.Omega_h;

    % Compute np.q_shift as a matrix where q(i)/y(j)
    q_shift = np.q ./ np.y';
    Omega_h_interp = interp1(np.q, eqm.Omega_h, q_shift, 'spline', 'extrap');
    out_of_range = np.q*(1+p.g) < min(np.q) | np.q*(1+p.g) > max(np.q);
    Omega_h_interp(out_of_range) = 1; % Use the last value as a constant
    temp_h2 = (1-p.gamma) * eqm.lambda_h * eqm.nu0_h .* (Omega_h_interp * (p.theta ./ (np.y).^(p.theta+1) .* np.dy));

    q_shift = np.q ./ np.x';
    Omega_h_interp = interp1(np.q, eqm.Omega_h, q_shift, 'spline', 'extrap');
    out_of_range = np.q*(1+p.g) < min(np.q) | np.q*(1+p.g) > max(np.q);
    Omega_h_interp(out_of_range) = 1; % Use the last value as a constant
    temp_h3 = (1-p.gamma) * eqm.mu0_hh * (Omega_h_interp * (p.theta * p.alpha^(p.theta) ./ (np.x).^(p.theta+1) .* np.dx));

    q_shift = np.q ./ np.x';
    Omega_l_interp = interp1(np.q, eqm.Omega_l, q_shift, 'spline', 'extrap');
    out_of_range = np.q*(1+p.g) < min(np.q) | np.q*(1+p.g) > max(np.q);
    Omega_l_interp(out_of_range) = 1; % Use the last value as a constant
    temp_h4 = (1-p.gamma) * eqm.mu0_lh * (Omega_l_interp * (p.theta * p.alpha^(p.theta) ./ (np.x).^(p.theta+1) .* np.dx));

    temp_h5 = p.gamma * p.m_h * (eqm.Omega_h+ eqm.Omega_l);


    temp_h = temp_h1 + temp_h2 + temp_h3 + temp_h4 + temp_h5 ;

    temp_l1 = (1-p.gamma) * (1-eqm.lambda_l) * eqm.nu0_l .* eqm.Omega_l;

    % Compute np.q_shift as a matrix where q(i)/y(j)
    q_shift = np.q ./ np.y';
    Omega_l_interp = interp1(np.q, eqm.Omega_l, q_shift, 'spline', 'extrap');
    out_of_range = np.q*(1+p.g) < min(np.q) | np.q*(1+p.g) > max(np.q);
    Omega_l_interp(out_of_range) =1; % Use the last value as a constant
    temp_l2 = (1-p.gamma) * eqm.lambda_l * eqm.nu0_l .* (Omega_l_interp * (p.theta ./ (np.y).^(p.theta+1) .* np.dy));


    q_shift = np.q ./ np.x';
    Omega_l_interp = interp1(np.q, eqm.Omega_l, q_shift, 'spline', 'extrap');
    out_of_range = np.q*(1+p.g) < min(np.q) | np.q*(1+p.g) > max(np.q);
    Omega_l_interp(out_of_range) = 1; % Use the last value as a constant
    temp_l3 = (1-p.gamma) * eqm.mu0_ll * (Omega_l_interp * (p.theta * p.alpha^(p.theta) ./ (np.x).^(p.theta+1) .* np.dx));

    q_shift = np.q ./ np.x';
    Omega_h_interp = interp1(np.q, eqm.Omega_h, q_shift, 'spline', 'extrap');
    out_of_range = np.q*(1+p.g) < min(np.q) | np.q*(1+p.g) > max(np.q);
    Omega_h_interp(out_of_range) = 1; % Use the last value as a constant
    temp_l4 = (1-p.gamma) * eqm.mu0_hl * (Omega_h_interp * (p.theta * p.alpha^(p.theta) ./ (np.x).^(p.theta+1) .* np.dx));

    temp_l5 = p.gamma * p.m_l * (eqm.Omega_l + eqm.Omega_h);

    temp_l = temp_l1 + temp_l2 + temp_l3 + temp_l4 + temp_l5 ;

    temp_h = temp_h/max(temp_h);
    temp_l = temp_l/max(temp_l);


    eqm.Omega_h_new = interp1(np.q, temp_h, np.q*(1+p.g), 'spline', 'extrap');
    eqm.Omega_l_new = interp1(np.q, temp_l, np.q*(1+p.g), 'spline', 'extrap');
    % Make monotonic increasing
    % eqm.Omega_h_new = cummax(eqm.Omega_h_new);
    % eqm.Omega_l_new = cummax(eqm.Omega_l_new);

    % % Rescale to end at 1
    eqm.Omega_h_new = eqm.Omega_h_new ./ eqm.Omega_h_new(end);
    eqm.Omega_l_new = eqm.Omega_l_new ./ eqm.Omega_l_new(end);

    % % Update
    eqm.Omega_h = eqm.Omega_h_new;
    eqm.Omega_l = eqm.Omega_l_new;

    % plot(np.q, eqm.Omega_h, 'r', 'LineWidth', 2); hold on;
    % plot(np.q, eqm.Omega_h_new, 'b', 'LineWidth', 2); hold on;

    %%

    q_shift = np.q .* np.y'/(1+p.g);
    eqm.v_h_interp = interp1(np.q, eqm.v_h, q_shift, 'spline', 'extrap');
    term_h1 =  eqm.nu0_h * eqm.v_h_interp * (p.theta ./ (np.y).^(p.theta+1) .* np.dy);

    q_shift = np.q /(1+p.g);
    eqm.v_h_interp = interp1(np.q, eqm.v_h, q_shift, 'spline', 'extrap');
    term_h2 = eqm.nu0_h * eqm.nu0_h * eqm.v_h_interp;


    q_shift = np.q .* np.x'/(1+p.g);
    eqm.v_h_interp = interp1(np.q, eqm.v_h, q_shift, 'spline', 'extrap');
    int_x = eqm.v_h_interp * (p.theta * p.alpha^(p.theta) ./ (np.x).^(p.theta+1) .* np.dx);
    int_omega = int_x' * f.forward_diff_Omega(eqm.Omega_h);
    term_h3 = eqm.s_h * f.P_X(p, p.Phi_Matrix(2,2)) * int_omega;

    q_shift = np.q .* np.x'/(1+p.g);
    eqm.v_h_interp = interp1(np.q, eqm.v_h, q_shift, 'spline', 'extrap');
    int_x = eqm.v_h_interp * (p.theta * p.alpha^(p.theta) ./ (np.x).^(p.theta+1) .* np.dx);
    int_omega = int_x' * f.forward_diff_Omega(eqm.Omega_l);
    term_h4 = eqm.s_l * f.P_X(p, p.Phi_Matrix(2,1)) * int_omega;


    eqm.delta_h = (   p.bellman * ( term_h3 + term_h4 )      ).^(1/(1+p.beta)) * ones(size(np.q));
    eqm.lambda_h = (  p.bellman * ( term_h1 - term_h2 )./eqm.v_h   ).^(1/(1+p.beta));


    q_shift = np.q .* np.y'/(1+p.g);
    eqm.v_l_interp = interp1(np.q, eqm.v_l, q_shift, 'spline', 'extrap');
    term_l1 =  eqm.nu0_l * eqm.v_l_interp * (p.theta ./ (np.y).^(p.theta+1) .* np.dy);

    q_shift = np.q /(1+p.g);
    eqm.v_l_interp = interp1(np.q, eqm.v_l, q_shift, 'spline', 'extrap');
    term_l2 = eqm.nu0_l * eqm.nu0_l * eqm.v_l_interp;


    q_shift = np.q .* np.x'/(1+p.g);
    eqm.v_l_interp = interp1(np.q, eqm.v_l, q_shift, 'spline', 'extrap');
    int_x = eqm.v_l_interp * (p.theta * p.alpha^(p.theta) ./ (np.x).^(p.theta+1) .* np.dx);
    int_omega = int_x' * f.forward_diff_Omega(eqm.Omega_h);
    term_l3 = eqm.s_l * f.P_X(p, p.Phi_Matrix(1,1)) * int_omega;

    q_shift = np.q .* np.x'/(1+p.g);
    eqm.v_l_interp = interp1(np.q, eqm.v_l, q_shift, 'spline', 'extrap');
    int_x = eqm.v_l_interp * (p.theta * p.alpha^(p.theta) ./ (np.x).^(p.theta+1) .* np.dx);
    int_omega = int_x' * f.forward_diff_Omega(eqm.Omega_l);
    term_l4 = eqm.s_h * f.P_X(p, p.Phi_Matrix(1,2)) * int_omega;

    eqm.delta_l = (   p.bellman * ( term_l3 + term_l4 )      ).^(1/(1+p.beta)) * ones(size(np.q));
    eqm.lambda_l = (  p.bellman * ( term_l1 + term_l2 )./eqm.v_l   ).^(1/(1+p.beta));

    eqm.delta_l = min( max(eqm.delta_l, 0), 1);
    eqm.delta_h = min( max(eqm.delta_h, 0), 1);

    eqm.lambda_l = min( max(eqm.lambda_l, 0), 1);
    eqm.lambda_h = min( max(eqm.lambda_h, 0), 1);

    eqm.v_l_init = p.Pi_Vec(1) * np.q.^(p.sigma - 1) ...
            - eqm.delta_l.^(1+p.beta)/(1+p.beta) ...
            - eqm.lambda_l.^(1+p.beta)/(1+p.beta) .* eqm.v_l ...
            + p.bellman * ( eqm.lambda_l .* term_l1 + (1 - eqm.lambda_l) .* term_l2 ...
                            + eqm.delta_l * term_l3   + eqm.delta_l * term_l4  );

    eqm.v_h_init = p.Pi_Vec(2) * np.q.^(p.sigma - 1) ...
            - eqm.delta_h.^(1+p.beta)/(1+p.beta) ...
            - eqm.lambda_h.^(1+p.beta)/(1+p.beta) .* eqm.v_h ...
            + p.bellman * ( eqm.lambda_h .* term_h1 + (1 - eqm.lambda_h) .* term_h2 ...
                            + eqm.delta_h * term_h3   + eqm.delta_h * term_h4  );
                            
    eqm.v_h_init = max(eqm.v_h_init, 0);
    eqm.v_l_init = max(eqm.v_l_init, 0);

    % eqm.v_h_init = min(eqm.v_h_init, 1e8);
    % eqm.v_l_init = min(eqm.v_l_init, 1e8);

    %%


    % disp(['Delta_l: ', num2str(max(abs(eqm.v_l_init - eqm.v_l)))])

    if mod(iter, 50) == 0 
        figure;
        subplot(4,2,1);
        plot(np.q, eqm.v_h_init, 'r', 'LineWidth', 2); hold on;
        title("$v_h$" , "Interpreter", "latex")
        xlabel('q')

        subplot(4,2,2);
        plot(np.q, eqm.v_l_init, 'g', 'LineWidth', 2); hold on;
        title("$v_l$" , "Interpreter", "latex")
        xlabel('q')

        subplot(4,2,3);

        plot(np.q, eqm.lambda_h, 'r', 'LineWidth', 2); hold on;
        title("$\lambda_h$" , "Interpreter", "latex")
        xlabel('q')

        subplot(4,2,4);
        plot(np.q, eqm.lambda_l, 'g', 'LineWidth', 2); hold on;
        title("$\lambda_l$" , "Interpreter", "latex")
        xlabel('q')

        subplot(4,2,5);

        plot(np.q, eqm.delta_h, 'r', 'LineWidth', 2); hold on;
        title("$\delta_h$" , "Interpreter", "latex")
        xlabel('q')

        subplot(4,2,6);
        plot(np.q, eqm.delta_l, 'g', 'LineWidth', 2); hold on;
        title("$\delta_l$" , "Interpreter", "latex")
        xlabel('q')

        subplot(4,2,7);
        plot(np.q, eqm.Omega_h, 'r', 'LineWidth', 2); hold on;
        title("$\Omega_h$" , "Interpreter", "latex")
        xlabel('q')

        subplot(4,2,8);
        plot(np.q, eqm.Omega_l, 'g', 'LineWidth', 2); hold on;
        title("$\Omega_l$" , "Interpreter", "latex")

        sgtitle(['Iteration: ', num2str(iter)])
        drawnow
        % Check convergence
        disp(['Iteration: ', num2str(iter)])
        disp(['v_h: ', num2str(max(abs(eqm.v_h_init - eqm.v_h)))])
        disp(['v_h: ', num2str(max(abs(eqm.v_l_init - eqm.v_l)))])

    end


end


% plot(np.q, eqm.v_h_init)
% plot(np.q, eqm.Omega_h)

x1 = [1,2,3];
x2 = [6,7,8];


x12 = x2 .* x1';

x12


x12_loop = zeros(length(x1), length(x2));

for i = 1:length(x1)
    for j = 1:length(x2)
        x12_loop(i,j) = x1(i) * x2(j);
    end
end


x12_loop



x12_x2_loop = zeros(length(x1), length(x2));

for i = 1:length(x1)
    for j = 1:length(x2)
        x12_x2_loop(i,j) = x12_loop(i,j) * x2(j);
    end
end

x12_x2 = x12_loop .* x2;




x12_x2_x1_loop = zeros(length(x1), length(x2));

for i = 1:length(x1)
    for j = 1:length(x2)
        x12_x2_x1_loop(i,j) = x12_loop(i,j) * x2(j) * x1(i);
    end
end

x12_x2_x1 = x12_loop .* x2 .* x1';


%%


clear; close all; clc


%%

p = EconomicsParameters();


func = Utility;


np = NumericalParameters(func,p);


eqm = Model_Initialization(p, np, func);

%%
np.iter_step_max = 80;
s = zeros(2, np.iter_step_max);
nu_hat = zeros(2, np.iter_step_max);
delta_hat = zeros(2, np.iter_step_max);
mu_hat = zeros(2, 2, np.iter_step_max);
g = zeros(1, np.iter_step_max);
gp1_type = zeros(2, np.iter_step_max);
% for iter = 1:np.iter_step_max






% plot(np.q, eqm.omega_tilde(1,:) )
% plot(np.q, eqm.omega(1,:) )




omega_old = eqm.omega;
eqm = Model_ExogenousProb(p, np, func, eqm);

eqm = Model_MarketShare(p, np, func, eqm);

eqm = Model_AggregateGrowth(p, np, func, eqm);

eqm = Model_LoM_PDF_integral_gp1(p, np, func, eqm);

eqm = Model_Bellman_integral_tamelambda(p, np, func, eqm);


% s(:,iter) = eqm.s;
% nu_hat(:,iter) = eqm.nu_hat;
% delta_hat(:,iter) = eqm.delta_hat;
% mu_hat(:,:,iter) = eqm.mu_hat;
% g(iter) = eqm.g;
% gp1_type(:,iter) = eqm.gp1_type;

diff = abs(eqm.omega - omega_old);
diff = max(diff(:));
% fprintf('Iteration %d, max diff = %.3e\n', iter, diff);
fprintf('max diff = %.3e\n', diff);



% figure;
subplot(1,3,1);
hold on;
plot(np.q, squeeze(eqm.omega(1,:)),'r','LineWidth',2);
plot(np.q, squeeze(eqm.omega(2,:)),'b','LineWidth',2);
plot(np.q, gampdf(np.q, 3, 1/2),'k','LineWidth',2);
hold off;
xlabel('q')
ylabel('$\omega$', 'Interpreter', 'latex')
title('$\omega(q)$', 'Interpreter', 'latex')
legend({'$\omega_l$', '$\omega_h$', 'orig'}, 'Interpreter', 'latex', 'Location', 'best')


subplot(1,3,2);
hold on;
plot(np.q, squeeze(eqm.lambda(1,:)),'r','LineWidth',2);
plot(np.q, squeeze(eqm.lambda(2,:)),'b','LineWidth',2);
hold off;
xlabel('q')
ylabel('$\lambda$', 'Interpreter', 'latex')
title('$\lambda(q)$', 'Interpreter', 'latex')
legend({'$\lambda_l$', '$\lambda_h$'}, 'Interpreter', 'latex', 'Location', 'best')

subplot(1,3,3);
hold on;
plot(np.q, squeeze(eqm.delta(1,:)),'r','LineWidth',2);
plot(np.q, squeeze(eqm.delta(2,:)),'b','LineWidth',2);
hold off;
xlabel('q')
ylabel('$\delta$', 'Interpreter', 'latex')
title('$\delta(q)$', 'Interpreter', 'latex')
legend({'$\delta_l$', '$\delta_h$'}, 'Interpreter', 'latex', 'Location', 'best')


sgtitle(sprintf ('$\\epsilon=%.1f, \\gamma=%.1f, \\sigma=%d, \\phi_h=%.2f, \\bar{m}=%.1f$', p.epsilon, p.gamma, p.sigma, p.phi(2), p.m(1)), 'Interpreter', 'latex')
set(gcf, 'Position', [100, 100, 1200, 900])

eqm.mu_hat









saveas(gcf, sprintf ('./figure/epsilon=%.1f, gamma=%.1f, sigma=%d, phi_h=%.2f, m_bar=%.1f.png', p.epsilon, p.gamma, p.sigma, p.phi(2), p.m(1)))




% figure(2);
% subplot(1,3,1);
% hold on;
% plot(np.q, squeeze(eqm.omega(1,:)),'r','LineWidth',2);
% plot(np.q, squeeze(eqm.omega(2,:)),'b','LineWidth',2);
% plot(np.q, gampdf(np.q, 3, 1/2),'k','LineWidth',2);
hold on;
plot(np.q, eqm.v(1,:),'r','LineWidth',2);
plot(np.q, eqm.v(2,:),'b','LineWidth',2); 
hold off;


%%


figure(2);
subplot(1,3,1);
hold on;
plot(np.q, squeeze(eqm.omega(1,:)),'r','LineWidth',2);
plot(np.q, squeeze(eqm.omega(2,:)),'b','LineWidth',2);
plot(np.q, gampdf(np.q, 3, 1/2),'k','LineWidth',2);
hold off;
xlabel('q')
ylabel('$\omega$', 'Interpreter', 'latex')
title('$\omega(q)$', 'Interpreter', 'latex')
legend({'$\omega_l$', '$\omega_h$', 'orig', 'exp(-q)'}, 'Interpreter', 'latex', 'Location', 'best')


subplot(1,3,2);
hold on;
plot(np.q, squeeze(eqm.lambda(1,:)),'r','LineWidth',2);
plot(np.q, squeeze(eqm.lambda(2,:)),'b','LineWidth',2);
hold off;
xlabel('q')
ylabel('$\lambda$', 'Interpreter', 'latex')
title('$\lambda(q)$', 'Interpreter', 'latex')
legend({'$\lambda_l$', '$\lambda_h$'}, 'Interpreter', 'latex', 'Location', 'best')

subplot(1,3,3);
hold on;
plot(np.q, squeeze(eqm.delta(1,:)),'r','LineWidth',2);
plot(np.q, squeeze(eqm.delta(2,:)),'b','LineWidth',2);
hold off;
xlabel('q')
ylabel('$\delta$', 'Interpreter', 'latex')
title('$\delta(q)$', 'Interpreter', 'latex')
legend({'$\delta_l$', '$\delta_h$'}, 'Interpreter', 'latex', 'Location', 'best')
sgtitle(sprintf ('$\\epsilon=%.1f, \\gamma=%.1f, \\sigma=%d, \\phi_h=%.2f, \\bar{m}=%.1f$', p.epsilon, p.gamma, p.sigma, p.phi(2), p.m(1)), 'Interpreter', 'latex')


%%



clear; close all; clc


%%

p = EconomicsParameters();


func = Utility;


np = NumericalParameters(func,p);


eqm = Model_Initialization(p, np, func);

%%
np.iter_step_max = 80;
s = zeros(2, np.iter_step_max);
nu_hat = zeros(2, np.iter_step_max);
delta_hat = zeros(2, np.iter_step_max);
mu_hat = zeros(2, 2, np.iter_step_max);
g = zeros(1, np.iter_step_max);
gp1_type = zeros(2, np.iter_step_max);
% for iter = 1:np.iter_step_max






% plot(np.q, eqm.omega_tilde(1,:) )
% plot(np.q, eqm.omega(1,:) )



iter =0;

iter= iter+1;
omega_old = eqm.omega;
eqm = Model_ExogenousProb(p, np, func, eqm);

eqm = Model_MarketShare(p, np, func, eqm);

eqm = Model_AggregateGrowth(p, np, func, eqm);

% eqm = Model_LoM_PDF_integral_gp1(p, np, func, eqm);
eqm = Model_LoM_PDF_integral_gp1_exogenous_lognormal(p, np, func, eqm);

% eqm = Model_Bellman_integral(p, np, func, eqm);
% eqm = Model_Bellman_integral_tamelambda(p, np, func, eqm);
eqm = Model_Bellman_integral_tamelambda_simple(p, np, func, eqm);




% s(:,iter) = eqm.s;
% nu_hat(:,iter) = eqm.nu_hat;
% delta_hat(:,iter) = eqm.delta_hat;
% mu_hat(:,:,iter) = eqm.mu_hat;
% g(iter) = eqm.g;
% gp1_type(:,iter) = eqm.gp1_type;

diff = abs(eqm.omega - omega_old);
diff = max(diff(:));
% fprintf('Iteration %d, max diff = %.3e\n', iter, diff);
fprintf('iter =%d,  max diff = %.3e\n', iter, diff);



% figure;
subplot(2,2,1);
hold on;
plot(np.q, squeeze(eqm.omega(1,:)),'r','LineWidth',2);
plot(np.q, squeeze(eqm.omega(2,:)),'b','LineWidth',2);
plot(np.q,   1./(np.q * (1+eqm.g) * p.tau * sqrt(2 * pi)) .* exp(- (log(np.q * (1+eqm.g)) - p.iota).^2 / (2 * p.tau^2)),'k','LineWidth',2);
hold off;
xlabel('q')
ylabel('$\omega$', 'Interpreter', 'latex')
title('$\omega(q)$', 'Interpreter', 'latex')
legend({'$\omega_l$', '$\omega_h$', 'orig'}, 'Interpreter', 'latex', 'Location', 'best')


subplot(2,2,2);
hold on;
plot(np.q, squeeze(eqm.lambda(1,:)),'r','LineWidth',2);
plot(np.q, squeeze(eqm.lambda(2,:)),'b','LineWidth',2);
hold off;
xlabel('q')
ylabel('$\lambda$', 'Interpreter', 'latex')
title('$\lambda(q)$', 'Interpreter', 'latex')
legend({'$\lambda_l$', '$\lambda_h$'}, 'Interpreter', 'latex', 'Location', 'best')

subplot(2,2,3);
hold on;
plot(np.q, squeeze(eqm.delta(1,:)),'r','LineWidth',2);
plot(np.q, squeeze(eqm.delta(2,:)),'b','LineWidth',2);
hold off;
xlabel('q')
ylabel('$\delta$', 'Interpreter', 'latex')
title('$\delta(q)$', 'Interpreter', 'latex')
legend({'$\delta_l$', '$\delta_h$'}, 'Interpreter', 'latex', 'Location', 'best')


% figure(1);
subplot(2,2,4);
hold on;
plot(np.q, squeeze(eqm.v(1,:)),'r','LineWidth',2);
plot(np.q, squeeze(eqm.v(2,:)),'b','LineWidth',2);
plot(np.q, p.pi(1) .* np.q.^(p.sigma-1),'k','LineWidth',2);
hold off;
xlabel('q')
ylabel('$v$', 'Interpreter', 'latex')
title('$v(q)$', 'Interpreter', 'latex')
legend({'$v_l$', '$v_h$', 'orig'}, 'Interpreter', 'latex', 'Location', 'best')

% subplot(2,3,5);
% hold on;
% plot(g, 'g','LineWidth',2);
% plot(gp1_type(1,:),'r','LineWidth',2);
% plot(gp1_type(2,:),'b','LineWidth',2);
% hold off;
% xlabel('Iteration','Interpreter','latex','FontSize',14);
% ylabel('Growth Rate','Interpreter','latex','FontSize',14);
% title('Growth Rate Convergence','Interpreter','latex','FontSize',14);
% legend({'$g$','$g_l$', '$g_h$'}, 'Interpreter', 'latex', 'Location', 'best')


% subplot(2,3,6);
% hold on;
% plot(s(1,:),'r','LineWidth',2);
% plot(s(2,:),'b','LineWidth',2);
% hold off;
% xlabel('Iteration','Interpreter','latex','FontSize',14);
% ylabel('Market Share','Interpreter','latex','FontSize',14);
% legend('$s_l$','$s_h$','Interpreter','latex','Location','best','FontSize',14);
% title('Market Share Convergence','Interpreter','latex','FontSize',14);

sgtitle(sprintf ('$\\epsilon=%.1f, \\gamma=%.1f, \\sigma=%d, \\phi_h=%.2f, \\bar{m}=%.1f$', p.epsilon, p.gamma, p.sigma, p.phi(2), p.m(1)), 'Interpreter', 'latex')
% set(gcf, 'Position', [100, 100, 1200, 900])

eqm.omega(:,end)







%%
% Parameters
qhat = linspace(0.01, 3, 1000);
q0 = 0.5;  % target mode of point distribution

% Set global defaults for LaTeX
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

% -------------------------------
% Figure setup with two subplots
% -------------------------------
figure;

% ---------- Subplot 1: Lognormal ----------
subplot(1, 2, 1);
hold on;

% Try multiple sigma values (fix mu = log(q0))
sigmas = [0.3, 0.1, 0.05];
colors = {'r', 'g', 'b'};
for i = 1:length(sigmas)
    sigma_logn = sigmas(i);
    mu_logn = log(q0);
    lognormal_pdf = lognpdf(qhat, mu_logn, sigma_logn);
    lognormal_pdf = lognormal_pdf / trapz(qhat, lognormal_pdf);
    plot(qhat, lognormal_pdf, 'Color', colors{i}, 'LineWidth', 2, ...
        'DisplayName', sprintf('$\\sigma = %.2f$', sigma_logn));
end

xline(q0, 'k--', 'LineWidth', 1.5, 'DisplayName', '$\phi$', 'Interpreter', 'latex');
title('Lognormal Approximations', 'Interpreter', 'latex', 'FontSize', fontSize);
xlabel('$\hat{q}$', 'Interpreter', 'latex', 'FontSize', fontSize);
ylabel('Density', 'Interpreter', 'latex', 'FontSize', fontSize);
legend('Location', 'northeast', 'FontSize', fontSize);
grid on;

% ---------- Subplot 2: Gamma ----------
subplot(1, 2, 2);
hold on;

% Try multiple theta values (fix mode at q0: (k - 1)*theta = q0)
thetas = [0.1, 0.03, 0.01];
for i = 1:length(thetas)
    theta_gamma = thetas(i);
    k_gamma = q0 / theta_gamma + 1;
    gamma_pdf = gampdf(qhat, k_gamma, theta_gamma);
    gamma_pdf = gamma_pdf / trapz(qhat, gamma_pdf);
    plot(qhat, gamma_pdf, 'Color', colors{i}, 'LineWidth', 2, ...
        'DisplayName', sprintf('$\\theta = %.2f$', theta_gamma));
end

xline(q0, 'k--', 'LineWidth', 1.5, 'DisplayName', '$\phi$', 'Interpreter', 'latex');
title('Lognormal Approximations', 'Interpreter', 'latex', 'FontSize', fontSize);
xlabel('$\hat{q}$', 'Interpreter', 'latex', 'FontSize', fontSize);
ylabel('Density', 'Interpreter', 'latex', 'FontSize', fontSize);
legend('Location', 'northeast', 'FontSize', fontSize);
grid on;

sgtitle('Approximating Point Distribution at $\phi = 0.5$','Interpreter', 'latex', 'FontSize', fontSize);


% Set figure size and position
set(gcf, 'Position', [100, 100, 1200, 800]);
% Set consistent font size
fontSize = 25;
set(groot, 'defaultAxesFontSize', fontSize);
set(groot, 'defaultLegendFontSize', fontSize);
set(groot, 'defaultTextFontSize', fontSize);

% Save figure
saveas(gcf, sprintf ('./figure/point_distribution_approximation.png'))



%%


% Parameters
alpha = 0.88;
theta = 5;

% Domains
x = linspace(alpha, 5, 500);
y = linspace(1, 5, 500);
z = linspace(0.01, 5, 500);

% CDF of X
F_X = 1 - (alpha ./ x).^theta;

% CDF of Y
F_Y = 1 - (1 ./ y).^theta;

% CDF of Z (piecewise)
F_Z = zeros(size(z));
idx1 = z < alpha;
idx2 = z >= alpha;
F_Z(idx1) = (z(idx1).^theta) / (2 * alpha^theta);
F_Z(idx2) = 1 - (alpha^theta) ./ (2 * z(idx2).^theta);

% Plot
figure; hold on; grid on;
plot(x, F_X, 'r-', 'LineWidth', 2);
plot(y, F_Y, 'b--', 'LineWidth', 2);
plot(z, F_Z, 'k-.', 'LineWidth', 2);
legend('$F_X(x)$', '$F_Y(y)$', '$F_Z(z)$', 'Interpreter', 'latex');
xlabel('$x, y, z$', 'Interpreter', 'latex');
ylabel('CDF Value', 'Interpreter', 'latex');
title('CDFs of X, Y, and Z (Pareto-related)', 'Interpreter', 'latex');






% Parameters
alpha = 0.88;
theta = 5;

% Domains
x = linspace(alpha, 5, 500);
y = linspace(1, 5, 500);
z = linspace(0.01, 5, 500);

% CDF of X
F_X = 1 - (alpha ./ x).^theta;

% CDF of Y
F_Y = 1 - (1 ./ y).^theta;

% CDF of Z (piecewise)
F_Z = zeros(size(z));
idx1 = z < alpha;
idx2 = z >= alpha;
F_Z(idx1) = (z(idx1).^theta) / (2 * alpha^theta);
F_Z(idx2) = 1 - (alpha^theta) ./ (2 * z(idx2).^theta);

% Plot
figure; hold on; grid on;
plot(x, 1-F_X, 'r-', 'LineWidth', 2);
plot(y, 1-F_Y, 'b--', 'LineWidth', 2);
plot(z, 1-F_Z, 'k-.', 'LineWidth', 2);
legend('$1-F_X(x)$', '$1-F_Y(y)$', '$1-F_Z(z)$', 'Interpreter', 'latex', 'FontSize', 20);
xlabel('$x, y, z$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('CDF Value', 'Interpreter', 'latex', 'FontSize', 20);
title('CDFs of X, Y, and Z (Pareto-related)', 'Interpreter', 'latex', 'FontSize', 20);