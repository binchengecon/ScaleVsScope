function main_func(p_exogenous)

    %%

    p = EconomicsParameters();


    p.m(1) = p_exogenous.m_bar; % m_l = 0.9, m_h = 0.1
    p.m(2) = 1 - p.m(1); % m_l = 0.9, m_h = 0.1
    p.gamma = p_exogenous.gamma;
    p.phi(2) = p_exogenous.phi_h;
    p.theta = p_exogenous.theta;
    p.alpha = p_exogenous.alpha;
    p.r = p_exogenous.r;
    p.beta = p_exogenous.beta;
    p.epsilon = p_exogenous.epsilon;

    func = Utility;


    np = NumericalParameters(func,p);


    eqm = Model_Initialization(p, np, func);

    %%
    np.iter_step_max = 800;
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

        eqm = Model_AggregateGrowth(p, np, func, eqm);

        eqm = Model_LoM_PDF(p, np, func, eqm);

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
            
            fprintf('Iteration %d, max diff = %.3e\n', iter, diff);

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
    eqm_save.g_share = eqm.g_share;

    iter_history.s = s;
    iter_history.nu_hat = nu_hat;
    iter_history.delta_hat = delta_hat;
    iter_history.mu_hat = mu_hat;
    iter_history.g = g;
    iter_history.gp1_type = gp1_type;

    % filepath: c:\Users\super\Dropbox\Esteban_Projects\Scale Vs Scope\1 code\matlab_bin_05_09_2025 LogNormal Wild_Defense\main.m
    title_name = sprintf ('$\\epsilon=%.2f, \\gamma=%.2f, \\sigma=%d, \\phi_h=%.2f, \\bar{m}=%.2f, r = %.2f, \\alpha=%.2f, \\theta=%.2f, \\beta=%.2f$', p.epsilon, p.gamma, p.sigma, p.phi(2), p.m(1), p.r, p.alpha, p.theta, p.beta);
    file_name = sprintf('epsilon=%.2f, gamma=%.2f, sigma=%d, phi_h=%.2f, m_bar=%.2f, r=%.2f, alpha=%.2f, theta=%.2f, beta=%.2f', p.epsilon, p.gamma, p.sigma, p.phi(2), p.m(1), p.r, p.alpha, p.theta, p.beta);

    data_folderpath = './data_beta/';
    figure_folderpath = './figure_beta/';

    if ~exist(figure_folderpath, 'dir')
        mkdir(figure_folderpath);
    end
    if ~exist(data_folderpath, 'dir')
        mkdir(data_folderpath);
    end

    save([data_folderpath, file_name, '.mat'] , 'p', 'np' , 'eqm_save', 'iter_history');
    fprintf('Save Success\n');

    % ...existing code...

    % Set a single variable for font size control
    fontSize = 18;      % Axis labels, ticks
    fontSizeTitle = 22; % Titles
    fontSizeLegend = 18; % Legends

    figure;
    subplot(2,4,1);
    hold on;
    plot(np.q, squeeze(eqm.omega(1,:)),'r','LineWidth',2);
    plot(np.q, squeeze(eqm.omega(2,:)),'b','LineWidth',2);
    hold off;
    xlabel('q','FontSize',fontSize)
    ylabel('$\omega$', 'Interpreter', 'latex','FontSize',fontSize)
    title('$\omega(q)$', 'Interpreter', 'latex','FontSize',fontSizeTitle)
    legend({'$\omega_l$', '$\omega_h$'}, 'Interpreter', 'latex', 'Location', 'best','FontSize',fontSizeLegend)

    subplot(2,4,2);
    hold on;
    plot(np.q, squeeze(eqm.lambda(1,:)),'r','LineWidth',2);
    plot(np.q, squeeze(eqm.lambda(2,:)),'b','LineWidth',2);
    hold off;
    xlabel('q','FontSize',fontSize)
    ylabel('$\lambda$', 'Interpreter', 'latex','FontSize',fontSize)
    title('$\lambda(q)$', 'Interpreter', 'latex','FontSize',fontSizeTitle)
    legend({'$\lambda_l$', '$\lambda_h$'}, 'Interpreter', 'latex', 'Location', 'best','FontSize',fontSizeLegend)

    subplot(2,4,3);
    hold on;
    plot(np.q, squeeze(eqm.delta(1,:)),'r','LineWidth',2);
    plot(np.q, squeeze(eqm.delta(2,:)),'b','LineWidth',2);
    hold off;
    xlabel('q','FontSize',fontSize)
    ylabel('$\delta$', 'Interpreter', 'latex','FontSize',fontSize)
    title('$\delta(q)$', 'Interpreter', 'latex','FontSize',fontSizeTitle)
    legend({'$\delta_l$', '$\delta_h$'}, 'Interpreter', 'latex', 'Location', 'best','FontSize',fontSizeLegend)

    subplot(2,4,4);
    hold on;
    plot(np.q, squeeze(eqm.v(1,:)),'r','LineWidth',2);
    plot(np.q, squeeze(eqm.v(2,:)),'b','LineWidth',2);
    hold off;
    xlabel('q','FontSize',fontSize)
    ylabel('$v$', 'Interpreter', 'latex','FontSize',fontSize)
    title('$v(q)$', 'Interpreter', 'latex','FontSize',fontSizeTitle)
    legend({'$v_l$', '$v_h$'}, 'Interpreter', 'latex', 'Location', 'best','FontSize',fontSizeLegend)

    subplot(2,4,5);
    hold on;
    plot(1:iter, g(1:iter), 'g','LineWidth',2);
    plot(1:iter, gp1_type(1,1:iter),'r','LineWidth',2);
    plot(1:iter, gp1_type(2,1:iter),'b','LineWidth',2);
    yline(0, 'k--', 'LineWidth', 1.5);
    hold off;
    xlabel('Iteration','Interpreter','latex','FontSize',fontSize);
    ylabel('Growth Rate','Interpreter','latex','FontSize',fontSize);
    title('Growth Rate Convergence','Interpreter','latex','FontSize',fontSizeTitle);
    legend({'$g$','$g_l$', '$g_h$'}, 'Interpreter', 'latex', 'Location', 'best','FontSize',fontSizeLegend)

    subplot(2,4,6);
    hold on;
    plot(1:iter, s(1,1:iter),'r','LineWidth',2);
    plot(1:iter, s(2,1:iter),'b','LineWidth',2);
    hold off;
    xlabel('Iteration','Interpreter','latex','FontSize',fontSize);
    ylabel('Market Share','Interpreter','latex','FontSize',fontSize);
    legend('$s_l$','$s_h$','Interpreter','latex','Location','best','FontSize',fontSizeLegend);
    title('Market Share Convergence','Interpreter','latex','FontSize',fontSizeTitle);

    subplot(2,4,7);
    hold on;
    plot(1:iter, nu_hat(1,1:iter),'r','LineWidth',2);
    plot(1:iter, nu_hat(2,1:iter),'b','LineWidth',2);
    hold off;
    xlabel('Iteration','Interpreter','latex','FontSize',fontSize);
    ylabel('Reserve Probability','Interpreter','latex','FontSize',fontSize);
    legend('$\nu_1$','$\nu_2$','Interpreter','latex','Location','best','FontSize',fontSizeLegend);
    title('Reserve Probability Convergence','Interpreter','latex','FontSize',fontSizeTitle);

    subplot(2,4,8);
    hold on;
    plot(1:iter, squeeze(mu_hat(1,1,1:iter)),'r','LineWidth',2);
    plot(1:iter, squeeze(mu_hat(1,2,1:iter)),'g','LineWidth',2);
    plot(1:iter, squeeze(mu_hat(2,1,1:iter)),'b','LineWidth',2);
    plot(1:iter, squeeze(mu_hat(2,2,1:iter)),'k','LineWidth',2);
    hold off;
    xlabel('Iteration','Interpreter','latex','FontSize',fontSize);
    ylabel('Steal Probability','Interpreter','latex','FontSize',fontSize);
    legend('$\mu_{11}$','$\mu_{12}$','$\mu_{21}$','$\mu_{22}$','Interpreter','latex','Location','best','FontSize',fontSizeLegend);
    title('Steal Probability Convergence','Interpreter','latex','FontSize',fontSizeTitle);

    sgtitle( title_name, 'Interpreter', 'latex', 'FontSize', fontSizeTitle);

    % set(gcf, 'Position', [100, 100, 2200, 1200]);
    % saveas(gcf, [figure_folderpath, file_name, sprintf('_iter=%d.png', iter)]);
    % ...existing code...
    % subplot(2,4,3);
    % hold on;
    % plot(delta_hat(1,:),'r','LineWidth',2);
    % plot(delta_hat(2,:),'b','LineWidth',2);
    % hold off;
    % xlabel('Iteration','Interpreter','latex','FontSize',14);
    % ylabel('Exogenous Probability','Interpreter','latex','FontSize',14);
    % legend('$\delta_1$','$\delta_2$','Interpreter','latex','Location','best','FontSize',14);
    % title('Exogenous Probability Convergence','Interpreter','latex','FontSize',14);


    % sgtitle( title_name, 'Interpreter', 'latex');
    set(gcf, 'Position', [100, 100, 2200, 1200]);


    saveas(gcf, [figure_folderpath, file_name, sprintf('_iter=%d.png', iter)]);






end