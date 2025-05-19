

clear; close all; clc


%%

p = EconomicsParameters();


% p.m(1) = p_exogenous.m_bar; % m_l = 0.9, m_h = 0.1
% p.m(2) = 1 - p.m(1); % m_l = 0.9, m_h = 0.1
% p.gamma = p_exogenous.gamma;
% % p.phi(2) = p_exogenous.phi_h;
% p.theta = p_exogenous.theta;
% p.alpha = p_exogenous.alpha;
% p.r = p_exogenous.r;

func = Utility;


np = NumericalParameters(func,p);


eqm = Model_Initialization(p, np, func);

%%
np.iter_step_max = 100;
s = zeros(2, np.iter_step_max);
nu_hat = zeros(2, np.iter_step_max);
delta_hat = zeros(2, np.iter_step_max);
mu_hat = zeros(2, 2, np.iter_step_max);
g = zeros(1, np.iter_step_max);
gp1_type = zeros(2, np.iter_step_max);

for iter = 1:np.iter_step_max

    omega_old = eqm.omega;
    
    eqm = Model_ExogenousProb(p, np, func, eqm);

    eqm = Model_MarketShare(p, np, func, eqm);

    % eqm = Model_AggregateGrowth(p, np, func, eqm);
    eqm = Model_AggregateGrowth_Discount(p, np, func, eqm);

    % eqm = Model_LoM_PDF(p, np, func, eqm);
    eqm = Model_LoM_PDF_Discount(p, np, func, eqm);

    eqm = Model_Bellman(p, np, func, eqm);




    s(:,iter) = eqm.s;
    nu_hat(:,iter) = eqm.nu_hat;
    delta_hat(:,iter) = eqm.delta_hat;
    mu_hat(:,:,iter) = eqm.mu_hat;
    g(iter) = eqm.g;
    gp1_type(:,iter) = eqm.gp1_type-1;


    diff = abs(eqm.omega - omega_old);
    diff = max(diff(:));

    if diff < 1e-8
        fprintf('Converged at iteration %d\n', iter);
        break;
    end

    if mod(iter, 10) == 0
        
        figure;
        subplot(2,3,1);
        hold on;
        plot(np.q, squeeze(eqm.omega(1,:)),'r','LineWidth',2);
        plot(np.q, squeeze(eqm.omega(2,:)),'b','LineWidth',2);
        plot(np.q,   1./(np.q * (1+eqm.g) * p.tau * sqrt(2 * pi)) .* exp(- (log(np.q * (1+eqm.g)) - p.iota).^2 / (2 * p.tau^2)),'k','LineWidth',2);
        hold off;
        xlabel('q')
        ylabel('$\omega$', 'Interpreter', 'latex')
        title('$\omega(q)$', 'Interpreter', 'latex')
        legend({'$\omega_l$', '$\omega_h$', 'orig'}, 'Interpreter', 'latex', 'Location', 'best')
        
        
        subplot(2,3,2);
        hold on;
        plot(np.q, squeeze(eqm.lambda(1,:)),'r','LineWidth',2);
        plot(np.q, squeeze(eqm.lambda(2,:)),'b','LineWidth',2);
        hold off;
        xlabel('q')
        ylabel('$\lambda$', 'Interpreter', 'latex')
        title('$\lambda(q)$', 'Interpreter', 'latex')
        legend({'$\lambda_l$', '$\lambda_h$'}, 'Interpreter', 'latex', 'Location', 'best')
        
        subplot(2,3,3);
        hold on;
        plot(np.q, squeeze(eqm.delta(1,:)),'r','LineWidth',2);
        plot(np.q, squeeze(eqm.delta(2,:)),'b','LineWidth',2);
        hold off;
        xlabel('q')
        ylabel('$\delta$', 'Interpreter', 'latex')
        title('$\delta(q)$', 'Interpreter', 'latex')
        legend({'$\delta_l$', '$\delta_h$'}, 'Interpreter', 'latex', 'Location', 'best')
        

        % figure(1);
        subplot(2,3,4);
        hold on;
        plot(np.q, squeeze(eqm.v(1,:)),'r','LineWidth',2);
        plot(np.q, squeeze(eqm.v(2,:)),'b','LineWidth',2);
        plot(np.q, p.pi(1) .* np.q.^(p.sigma-1),'k','LineWidth',2);
        hold off;
        xlabel('q')
        ylabel('$v$', 'Interpreter', 'latex')
        title('$v(q)$', 'Interpreter', 'latex')
        legend({'$v_l$', '$v_h$', 'orig'}, 'Interpreter', 'latex', 'Location', 'best')
        
        subplot(2,3,5);
        hold on;
        plot(g, 'g','LineWidth',2);
        plot(gp1_type(1,:),'r','LineWidth',2);
        plot(gp1_type(2,:),'b','LineWidth',2);
        hold off;
        xlabel('Iteration','Interpreter','latex','FontSize',14);
        ylabel('Growth Rate','Interpreter','latex','FontSize',14);
        title('Growth Rate Convergence','Interpreter','latex','FontSize',14);
        legend({'$g$','$g_l$', '$g_h$'}, 'Interpreter', 'latex', 'Location', 'best')


        subplot(2,3,6);
        hold on;
        plot(s(1,:),'r','LineWidth',2);
        plot(s(2,:),'b','LineWidth',2);
        hold off;
        xlabel('Iteration','Interpreter','latex','FontSize',14);
        ylabel('Market Share','Interpreter','latex','FontSize',14);
        legend('$s_l$','$s_h$','Interpreter','latex','Location','best','FontSize',14);
        title('Market Share Convergence','Interpreter','latex','FontSize',14);

        
        sgtitle(sprintf ('$\\epsilon=%.1f, \\gamma=%.1f, \\sigma=%d, \\phi_h=%.2f, \\bar{m}=%.1f$', p.epsilon, p.gamma, p.sigma, p.phi(2), p.m(1)), 'Interpreter', 'latex');
        set(gcf, 'Position', [100, 100, 1600, 1200]);

        fprintf('Iteration %d, max diff = %.3e\n', iter, diff);

        % saveas(gcf, sprintf ('./figure2/epsilon=%.1f, gamma=%.1f, sigma=%d, phi_h=%.2f, m_bar=%.1f_iter=%d.png', p.epsilon, p.gamma, p.sigma, p.phi(2), p.m(1), iter));

    end

end


























p.error = diff;
p.iter = iter;

eqm_save.v = eqm.v;
eqm_save.omega = eqm.omega;
eqm_save.omega_tilde = eqm.omega_tilde;
eqm_save.lambda = eqm.lambda;
eqm_save.delta = eqm.delta;
eqm_save.s = eqm.s;
eqm_save.g = eqm.g;
eqm_save.gp1_type = eqm.gp1_type;
eqm_save.mu_hat = eqm.mu_hat;
eqm_save.nu_hat = eqm.nu_hat;
eqm_save.delta_hat = eqm.delta_hat;

% filepath: c:\Users\super\Dropbox\Esteban_Projects\Scale Vs Scope\1 code\matlab_bin_05_09_2025 LogNormal Wild_Defense\main.m
title_name = sprintf ('$\\epsilon=%.2f, \\gamma=%.2f, \\sigma=%d, \\phi_h=%.2f, \\bar{m}=%.2f, r = %.2f$', p.epsilon, p.gamma, p.sigma, p.phi(2), p.m(1), p.r);
file_name = sprintf('./data/epsilon=%.2f, gamma=%.2f, sigma=%d, phi_h=%.2f, m_bar=%.2f, r=%.2f', p.epsilon, p.gamma, p.sigma, p.phi(2), p.m(1), p.r);
save([file_name, '.mat'] , 'p', 'np' , 'eqm_save');


figure;
subplot(2,3,1);
hold on;
plot(np.q, squeeze(eqm.omega(1,:)),'r','LineWidth',2);
plot(np.q, squeeze(eqm.omega(2,:)),'b','LineWidth',2);
plot(np.q,   1./(np.q * (1+eqm.g) * p.tau * sqrt(2 * pi)) .* exp(- (log(np.q * (1+eqm.g)) - p.iota).^2 / (2 * p.tau^2)),'k','LineWidth',2);
hold off;
xlabel('q')
ylabel('$\omega$', 'Interpreter', 'latex')
title('$\omega(q)$', 'Interpreter', 'latex')
legend({'$\omega_l$', '$\omega_h$', 'orig'}, 'Interpreter', 'latex', 'Location', 'best')


subplot(2,3,2);
hold on;
plot(np.q, squeeze(eqm.lambda(1,:)),'r','LineWidth',2);
plot(np.q, squeeze(eqm.lambda(2,:)),'b','LineWidth',2);
hold off;
xlabel('q')
ylabel('$\lambda$', 'Interpreter', 'latex')
title('$\lambda(q)$', 'Interpreter', 'latex')
legend({'$\lambda_l$', '$\lambda_h$'}, 'Interpreter', 'latex', 'Location', 'best')

subplot(2,3,3);
hold on;
plot(np.q, squeeze(eqm.delta(1,:)),'r','LineWidth',2);
plot(np.q, squeeze(eqm.delta(2,:)),'b','LineWidth',2);
hold off;
xlabel('q')
ylabel('$\delta$', 'Interpreter', 'latex')
title('$\delta(q)$', 'Interpreter', 'latex')
legend({'$\delta_l$', '$\delta_h$'}, 'Interpreter', 'latex', 'Location', 'best')


% figure(1);
subplot(2,3,4);
hold on;
plot(np.q, squeeze(eqm.v(1,:)),'r','LineWidth',2);
plot(np.q, squeeze(eqm.v(2,:)),'b','LineWidth',2);
plot(np.q, p.pi(1) .* np.q.^(p.sigma-1),'k','LineWidth',2);
hold off;
xlabel('q')
ylabel('$v$', 'Interpreter', 'latex')
title('$v(q)$', 'Interpreter', 'latex')
legend({'$v_l$', '$v_h$', 'orig'}, 'Interpreter', 'latex', 'Location', 'best')

subplot(2,3,5);
hold on;
plot(1:iter, g(1:iter), 'g','LineWidth',2);
plot(1:iter, gp1_type(1,1:iter),'r','LineWidth',2);
plot(1:iter, gp1_type(2,1:iter),'b','LineWidth',2);
hold off;
xlabel('Iteration','Interpreter','latex','FontSize',14);
ylabel('Growth Rate','Interpreter','latex','FontSize',14);
title('Growth Rate Convergence','Interpreter','latex','FontSize',14);
legend({'$g$','$g_l$', '$g_h$'}, 'Interpreter', 'latex', 'Location', 'best')


subplot(2,3,6);
hold on;
plot(1:iter, s(1,1:iter),'r','LineWidth',2);
plot(1:iter, s(2,1:iter),'b','LineWidth',2);
hold off;
xlabel('Iteration','Interpreter','latex','FontSize',14);
ylabel('Market Share','Interpreter','latex','FontSize',14);
legend('$s_l$','$s_h$','Interpreter','latex','Location','best','FontSize',14);
title('Market Share Convergence','Interpreter','latex','FontSize',14);


sgtitle( title_name, 'Interpreter', 'latex');
set(gcf, 'Position', [100, 100, 1600, 1200]);


saveas(gcf, ['./figure_converge/', file_name, sprintf('_iter=%d.png', iter)]);





