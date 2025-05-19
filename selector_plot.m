function selector_plot(p)

    file_name = sprintf('epsilon=%.2f, gamma=%.2f, sigma=%d, phi_h=%.2f, m_bar=%.2f, r=%.2f, alpha=%.2f, theta=%.2f, beta=%.2f', p.epsilon, p.gamma, p.sigma, p.phi_h, p.m_bar, p.r, p.alpha, p.theta, p.beta);
    % data = load(['./data_beta/', file_name, '.mat'] , 'p', 'np' , 'eqm_save', 'iter_history');

    data = load(['./data_beta/', file_name, '.mat'] , 'p', 'np' , 'eqm_save', 'iter_history');
    figure_folderpath = './figure_scenario/';


    if ~exist(figure_folderpath, 'dir')
        mkdir(figure_folderpath);
    end

    p_iter = data.p;
    np_iter = data.np;
    eqm_iter = data.eqm_save;
    iter_history = data.iter_history;
    
    iter = p_iter.iter;

    % Set a single variable for font size control
    fontSize = 18;      % Axis labels, ticks
    fontSizeTitle = 22; % Titles
    fontSizeLegend = 18; % Legends

    title_name = sprintf ('$\\epsilon=%.2f, \\gamma=%.2f, \\sigma=%d, \\phi_h=%.2f, \\bar{m}=%.2f, r = %.2f, \\alpha=%.2f, \\theta=%.2f, \\beta=%.2f$', p_iter.epsilon, p_iter.gamma, p_iter.sigma, p_iter.phi(2), p_iter.m(1), p_iter.r, p_iter.alpha, p_iter.theta, p_iter.beta);
    % idx = find(np_iter.q <= 10);

    figure;
    subplot(2,4,1);
    hold on;
    plot(np_iter.q, eqm_iter.omega(1,:),'r','LineWidth',2);
    plot(np_iter.q, eqm_iter.omega(2,:),'b','LineWidth',2);
    xlim([0 8]);
    hold off;
    grid on;
    xlabel('q','FontSize',fontSize)
    ylabel('$\omega$', 'Interpreter', 'latex','FontSize',fontSize)
    title('$\omega(q)$', 'Interpreter', 'latex','FontSize',fontSizeTitle)
    legend({'$\omega_l$', '$\omega_h$'}, 'Interpreter', 'latex', 'Location', 'best','FontSize',fontSizeLegend)

    subplot(2,4,2);
    hold on;
    plot(np_iter.q, eqm_iter.lambda(1,:),'r','LineWidth',2);
    plot(np_iter.q, eqm_iter.lambda(2,:),'b','LineWidth',2);
    ylim([0 1]);
    xlim([0 8]);
    grid on;
    hold off;
    xlabel('q','FontSize',fontSize)
    ylabel('$\lambda$', 'Interpreter', 'latex','FontSize',fontSize)
    title('$\lambda(q)$', 'Interpreter', 'latex','FontSize',fontSizeTitle)
    legend({'$\lambda_l$', '$\lambda_h$'}, 'Interpreter', 'latex', 'Location', 'best','FontSize',fontSizeLegend)

    subplot(2,4,3);
    hold on;
    plot(np_iter.q, eqm_iter.delta(1,:),'r','LineWidth',2);
    plot(np_iter.q, eqm_iter.delta(2,:),'b','LineWidth',2);
    ylim([0.5 1]);
    xlim([0 8]);
    grid on;
    hold off;
    xlabel('q','FontSize',fontSize)
    ylabel('$\delta$', 'Interpreter', 'latex','FontSize',fontSize)
    title('$\delta(q)$', 'Interpreter', 'latex','FontSize',fontSizeTitle)
    legend({'$\delta_l$', '$\delta_h$'}, 'Interpreter', 'latex', 'Location', 'best','FontSize',fontSizeLegend)

    subplot(2,4,4);
    hold on;
    plot(np_iter.q, eqm_iter.v(1,:),'r','LineWidth',2);
    plot(np_iter.q, eqm_iter.v(2,:),'b','LineWidth',2);
    xlim([0 8]);
    grid on;
    hold off;
    xlabel('q','FontSize',fontSize)
    ylabel('$v$', 'Interpreter', 'latex','FontSize',fontSize)
    title('$v(q)$', 'Interpreter', 'latex','FontSize',fontSizeTitle)
    legend({'$v_l$', '$v_h$'}, 'Interpreter', 'latex', 'Location', 'best','FontSize',fontSizeLegend)

    subplot(2,4,5);
    hold on;
    plot(1:iter, iter_history.g(1:iter), 'g','LineWidth',2);
    plot(1:iter, iter_history.gp1_type(1,1:iter),'r','LineWidth',2);
    plot(1:iter, iter_history.gp1_type(2,1:iter),'b','LineWidth',2);
    yline(0, 'k--', 'LineWidth', 1.5);
    ylim([0.0 0.3])
    grid on;
    hold off;
    xlabel('Iteration','Interpreter','latex','FontSize',fontSize);
    ylabel('Growth Rate','Interpreter','latex','FontSize',fontSize);
    title('Growth Rate Convergence','Interpreter','latex','FontSize',fontSizeTitle);
    legend({'$g$','$g_l$', '$g_h$'}, 'Interpreter', 'latex', 'Location', 'best','FontSize',fontSizeLegend)

    subplot(2,4,6);
    hold on;
    plot(1:iter, iter_history.s(1,1:iter),'r','LineWidth',2);
    plot(1:iter, iter_history.s(2,1:iter),'b','LineWidth',2);
    ylim([0.3 0.7])
    grid on;
    hold off;
    xlabel('Iteration','Interpreter','latex','FontSize',fontSize);
    ylabel('Market Share','Interpreter','latex','FontSize',fontSize);
    legend('$s_l$','$s_h$','Interpreter','latex','Location','best','FontSize',fontSizeLegend);
    title('Market Share Convergence','Interpreter','latex','FontSize',fontSizeTitle);

    subplot(2,4,7);
    hold on;
    plot(1:iter, iter_history.nu_hat(1,1:iter),'r','LineWidth',2);
    plot(1:iter, iter_history.nu_hat(2,1:iter),'b','LineWidth',2);
    grid on;
    hold off;
    xlabel('Iteration','Interpreter','latex','FontSize',fontSize);
    ylabel('Reserve Probability','Interpreter','latex','FontSize',fontSize);
    legend('$\nu_l$','$\nu_h$','Interpreter','latex','Location','best','FontSize',fontSizeLegend);
    title('Reserve Probability Convergence','Interpreter','latex','FontSize',fontSizeTitle);

    subplot(2,4,8);
    hold on;
    plot(1:iter, squeeze(iter_history.mu_hat(1,1,1:iter)),'r','LineWidth',2);
    plot(1:iter, squeeze(iter_history.mu_hat(1,2,1:iter)),'g','LineWidth',2);
    plot(1:iter, squeeze(iter_history.mu_hat(2,1,1:iter)),'b','LineWidth',2);
    plot(1:iter, squeeze(iter_history.mu_hat(2,2,1:iter)),'k','LineWidth',2);
    grid on;
    hold off;
    ylim([0 1]);
    xlabel('Iteration','Interpreter','latex','FontSize',fontSize);
    ylabel('Steal Probability','Interpreter','latex','FontSize',fontSize);
    legend('$\mu_{ll}$','$\mu_{lh}$','$\mu_{hl}$','$\mu_{hh}$','Interpreter','latex','Location','northeast','FontSize',fontSizeLegend);
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


    saveas(gcf, [figure_folderpath, file_name, 'Scenario.png']);



    
end