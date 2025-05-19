function p  = EconomicsParameters()


    % p.alpha = 0.88;
    p.alpha = .90;
    p.beta = 8; 
    % p.gamma = 0.1; % Exogenous Death Rate
    p.gamma = 0.2; % Exogenous Death Rate
    p.theta = 5;
    p.phi = [1, 1.05]; % phi_l = 1, phi_h = 1.05
    % p.phi = [1, 1.05]; % phi_l = 1, phi_h = 1.05
    % p.phi = [1, 1.05]; % phi_l = 1, phi_h = 1.05
    p.m = zeros(1,2); % m_l = 0.9, m_h = 0.1
    p.m(1) = 0.1; % m_l = 0.9, m_h = 0.1
    p.m(2) = 1 - p.m(1); % m_l = 0.9, m_h = 0.1
    % p.m = [0.6, 0.4]; % m_l = 0.9, m_h = 0.1
    p.sigma = 3;
    p.r = 0.15;
    % p.phi_ratio = p.phi'./p.phi;
    p.pi = 1/(p.sigma - 1) * ( (p.sigma-1)/ p.sigma)^p.sigma ./ (p.phi).^(p.sigma - 1);
    p.epsilon = 1;
    p.zeta= 0.1;
    % figure(2);


    p.tau = .5;
    p.iota = - (p.sigma - 1) * p.tau^2/2;

    % test= linspace(0.01, 10, 1000);

    % omega_entry = 1./(test * p.tau * sqrt(2 * pi)) .* exp(- (log(test) - p.iota).^2 / (2 * p.tau^2));

    % plot(test, omega_entry, 'r', 'LineWidth', 2);




    % % Range of tau values
    % tau_values = linspace(0.2, 1.5, 7); % 10 values from 0.1 to 2
    % test = linspace(0.01, 10, 1000); % Test range

    % % Prepare figure
    % figure;
    % hold on;

    % % Loop over tau values
    % for tau = tau_values
    %     p.tau = tau;
    %     p.iota = - (p.sigma - 1) * p.tau^2 / 2;

    %     % Compute omega_entry
    %     omega_entry = 1 ./ (test * p.tau * sqrt(2 * pi)) .* exp(- (log(test) - p.iota).^2 / (2 * p.tau^2));

    %     % Plot the curve
    %     plot(test, omega_entry, 'LineWidth', 2, 'DisplayName', sprintf('\\tau = %.2f', tau));
    % end

    % % Customize plot with larger font sizes
    % xlabel('$\hat{q}$', 'FontSize', 20, 'Interpreter', 'latex'); % Larger font size for x-axis label
    % ylabel('$\omega^{entry}(\hat{q})$', 'FontSize', 20, 'Interpreter', 'latex'); % Larger font size for y-axis label
    % title('$\omega^{entry}(\hat{q})$', 'FontSize', 20, 'Interpreter', 'latex'); % Larger font size for y-axis label
    % legend show;
    % set(legend, 'FontSize', 20); % Larger font size for legend
    % grid on;
    % hold off;
    % set(gcf, 'Position', [100, 100, 1600, 1200]); % Set figure size

    % saveas(gcf, './figure/omega_entry.png');


end