

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





%%



% gamma = [0.02 0.03 0.07];
% gamma = [0.07];
% gamma = [0.10];
% gamma = [0.20];
gamma = [0.10 0.20];

m = [0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9];
% m = [0.7, 0.75, 0.8, 0.85, 0.9];
% m = [0.75, 0.8, 0.85, 0.9];
% m = [0.85, 0.9];
% m = [0.9];

% phi_h = [1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3];
phi_h = [1.05];

% r = [0.05 0.1 0.15 0.2 0.25 0.3];
r = [0.10 0.15];
% theta = [4 5];
% theta = [3 4];
theta = [5];
% alpha = [0.5 0.8 0.88];
alpha = [0.8 0.9];
% alpha = [0.5 0.6 0.7 0.8 0.9];
epsilon = [0.0, 1.0]

result_Expected_q = zeros(length(gamma), length(m), length(r), length(theta), length(alpha),3);
result_10quantile_q = zeros(length(gamma), length(m), length(r), length(theta), length(alpha),3);
result_50quantile_q = zeros(length(gamma), length(m), length(r), length(theta), length(alpha),3);
result_90quantile_q = zeros(length(gamma), length(m), length(r), length(theta), length(alpha),3);
result_s = zeros(length(gamma), length(m), length(r), length(theta), length(alpha),2);
result_g_3 = zeros(length(gamma), length(m), length(r), length(theta), length(alpha),3);
result_mu_hat = zeros(length(gamma), length(m), length(r), length(theta), length(alpha),2,2);
result_nu_hat = zeros(length(gamma), length(m), length(r), length(theta), length(alpha),2);
result_delta_hat = zeros(length(gamma), length(m), length(r), length(theta), length(alpha),2);
result_lambda_hat = zeros(length(gamma), length(m), length(r), length(theta), length(alpha),2);
result_flag = zeros(length(gamma), length(m), length(r), length(theta), length(alpha),1);

for i = 1:length(gamma)
    for j = 1:length(m)
        % for k = 1:length(phi_h)
            for l = 1:length(r)
                for f = 1:length(theta)
                    for g = 1:length(alpha)
                        
                        p.r = r(l);
                        p.alpha = alpha(g);
                        p.theta = theta(f);
                        p.m_bar = m(j);
                        p.gamma = gamma(i);
                        p.phi_h = 1.05;
                        p.epsilon = 1.0;
                        p.sigma =3.0;
                        file_name = sprintf('epsilon=%.2f, gamma=%.2f, sigma=%d, phi_h=%.2f, m_bar=%.2f, r=%.2f, alpha=%.2f, theta=%.2f', p.epsilon, p.gamma, p.sigma, p.phi_h, p.m_bar, p.r, p.alpha, p.theta);
                        data = load(['./data3/', file_name, '.mat'] , 'p', 'np' , 'eqm_save');
                        % Regular expression to extract the iteration number

                        % Directory containing figure files
                        folder_path = './figure_converge3';

                        % Get list of all PNG files in the folder
                        file_list = dir(fullfile(folder_path, '*.png'));

                        % Loop through each file and extract iteration number
                        for k = 1:length(file_list)
                            file_name_var = file_list(k).name;
                            
                            % Check if the file name starts with the base file name
                            if startsWith(file_name_var, file_name)
                                % Regular expression to extract the iteration number
                                pattern = '_iter=(\d+)\.png';
                                tokens = regexp(file_name_var, pattern, 'tokens');
                                
                                % Convert the extracted token to a number
                                if ~isempty(tokens)
                                    iter = str2double(tokens{1}{1});
                                    fprintf('File: %s, Extracted iteration number: %d\n', file_name_var, iter);
                                    break; % Exit the loop once the iteration number is found

                                else
                                    fprintf('File: %s, No iteration number found.\n', file_name_var);
                                end
                            end
                        end


                        p_iter = data.p;
                        np_iter = data.np;
                        eqm_iter = data.eqm_save;
                        
                        lambda_interp_l = @(q) interp1(np_iter.q, eqm_iter.lambda(1,:), q, 'spline', 'extrap');
                        lambda_interp_h = @(q) interp1(np_iter.q, eqm_iter.lambda(2,:), q, 'spline', 'extrap');
                        omega_interp_l = @(q) interp1(np_iter.q, eqm_iter.omega(1,:), q, 'spline', 'extrap');
                        omega_interp_h = @(q) interp1(np_iter.q, eqm_iter.omega(2,:), q, 'spline', 'extrap');
                        omega_entry = @(q) 1./(q * p_iter.tau * sqrt(2*pi)) .* exp(-(log(q) - p_iter.iota).^2 / (2 * p_iter.tau^2));

                        result_Expected_q(i,j, l,f,g, 1) = integral(@(q) q .* omega_interp_l(q), np_iter.q_min, np_iter.q_max);  % Choose appropriate limits
                        result_10quantile_q(i,j, l,f,g, 1) = fzero(@(q) integral(omega_interp_l, np_iter.q_min, q) - 0.10, [np_iter.q_min, np_iter.q_max]);
                        result_50quantile_q(i,j, l,f,g, 1) = fzero(@(q) integral(omega_interp_l, np_iter.q_min, q) - 0.50, [np_iter.q_min, np_iter.q_max]);
                        result_90quantile_q(i,j, l,f,g, 1) = fzero(@(q) integral(omega_interp_l, np_iter.q_min, q) - 0.90, [np_iter.q_min, np_iter.q_max]);



                        result_Expected_q(i,j, l,f,g, 2) = integral(@(q) q .* omega_interp_h(q), np_iter.q_min, np_iter.q_max);  % Choose appropriate limits
                        result_10quantile_q(i,j, l,f,g, 2) = fzero(@(q) integral(omega_interp_h, np_iter.q_min, q) - 0.10, [np_iter.q_min, np_iter.q_max]);
                        result_50quantile_q(i,j, l,f,g, 2) = fzero(@(q) integral(omega_interp_h, np_iter.q_min, q) - 0.50, [np_iter.q_min, np_iter.q_max]);
                        result_90quantile_q(i,j, l,f,g, 2) = fzero(@(q) integral(omega_interp_h, np_iter.q_min, q) - 0.90, [np_iter.q_min, np_iter.q_max]);



                        result_Expected_q(i,j, l,f,g, 3) = integral(@(q) q .* omega_entry(q), np_iter.q_min, np_iter.q_max);  % Choose appropriate limits
                        result_10quantile_q(i,j, l,f,g, 3) = fzero(@(q) integral(omega_entry, np_iter.q_min, q) - 0.10, [np_iter.q_min, np_iter.q_max]);
                        result_50quantile_q(i,j, l,f,g, 3) = fzero(@(q) integral(omega_entry, np_iter.q_min, q) - 0.50, [np_iter.q_min, np_iter.q_max]);
                        result_90quantile_q(i,j, l,f,g, 3) = fzero(@(q) integral(omega_entry, np_iter.q_min, q) - 0.90, [np_iter.q_min, np_iter.q_max]);


                        % result_s = zeros(length(gamma), length(m), length(r), length(theta), length(alpha),2);
                        % result_g_3 = zeros(length(gamma), length(m), length(r), length(theta), length(alpha),3);
                        % result_mu_hat = zeros(length(gamma), length(m), length(r), length(theta), length(alpha),2,2);
                        % result_nu_hat = zeros(length(gamma), length(m), length(r), length(theta), length(alpha),2);
                        % result_delta_hat = zeros(length(gamma), length(m), length(r), length(theta), length(alpha),2);
                        % result_flag = zeros(length(gamma), length(m), length(r), length(theta), length(alpha),1);

                        result_s(i,j, l,f,g, :) = eqm_iter.s;
                        % result_g_3(i,j, l,f,g, 1) = (eqm_iter.gp1_type(1))^(1/p_iter.sigma-1)-1;
                        % result_g_3(i,j, l,f,g, 2) = (eqm_iter.gp1_type(2))^(1/p_iter.sigma-1)-1;
                        result_g_3(i,j, l,f,g, 1) = (eqm_iter.gp1_type(1))-1;
                        result_g_3(i,j, l,f,g, 2) = (eqm_iter.gp1_type(2))-1;
                        result_g_3(i,j, l,f,g, 3) = eqm_iter.g;

                        result_mu_hat(i,j, l,f,g, :, :) = eqm_iter.mu_hat;
                        result_nu_hat(i,j, l,f,g, :) = eqm_iter.nu_hat;
                        result_delta_hat(i,j, l,f,g, :) = eqm_iter.delta_hat;
                        result_lambda_hat(i,j, l,f,g, 1) = integral(@(q) omega_interp_l(q) .* lambda_interp_l(q), np_iter.q_min, np_iter.q_max);
                        result_lambda_hat(i,j, l,f,g, 2) = integral(@(q) omega_interp_h(q) .* lambda_interp_h(q), np_iter.q_min, np_iter.q_max);


                        if iter<500
                            result_flag(i,j, l,f,g, 1) = 1;
                        else
                            result_flag(i,j, l,f,g, 1) = 0;
                        end


                    end
                end
            end
        % end
    end
end



% l = 1;
% f = 1;
% g = 1;
% i = 1;
% j = 1;

% p.r = r(l);
% p.alpha = alpha(g);
% p.theta = theta(f);
% p.m_bar = m(j);
% p.gamma = gamma(i);
% p.phi_h = 0.95;
% p.epsilon = 1.0;
% p.sigma =3.0;
% file_name = sprintf('epsilon=%.2f, gamma=%.2f, sigma=%d, phi_h=%.2f, m_bar=%.2f, r=%.2f, alpha=%.2f, theta=%.2f', p.epsilon, p.gamma, p.sigma, p.phi_h, p.m_bar, p.r, p.alpha, p.theta);
% data = load(['./data3/', file_name, '.mat'] , 'p', 'np' , 'eqm_save');
% % Regular expression to extract the iteration number

% % Directory containing figure files
% folder_path = './figure_converge3';

% % Get list of all PNG files in the folder
% file_list = dir(fullfile(folder_path, '*.png'));

% % Loop through each file and extract iteration number
% for k = 1:length(file_list)
%     file_name_var = file_list(k).name;
    
%     % Check if the file name starts with the base file name
%     if startsWith(file_name_var, file_name)
%         % Regular expression to extract the iteration number
%         pattern = '_iter=(\d+)\.png';
%         tokens = regexp(file_name_var, pattern, 'tokens');
        
%         % Convert the extracted token to a number
%         if ~isempty(tokens)
%             iter = str2double(tokens{1}{1});
%             fprintf('File: %s, Extracted iteration number: %d\n', file_name_var, iter);
%             break; % Exit the loop once the iteration number is found

%         else
%             fprintf('File: %s, No iteration number found.\n', file_name_var);
%         end
%     end
% end


% p_iter = data.p;
% np_iter = data.np;
% eqm_iter = data.eqm_save;

% % omega_entry = @(q) 1./(q * p_iter.tau * sqrt(2*pi)) .* exp(-(log(q) - p_iter.iota).^2 / (2 * p_iter.tau^2));

% % E_q = integral(@(q) q .* omega_entry(q), np_iter.q_min, np_iter.q_max);  % Choose appropriate limits
% % q_10 = fzero(@(q) integral(omega_entry, np_iter.q_min, q) - 0.10, [np_iter.q_min, np_iter.q_max]);
% % q_50 = fzero(@(q) integral(omega_entry, np_iter.q_min, q) - 0.50, [np_iter.q_min, np_iter.q_max]);
% % q_90 = fzero(@(q) integral(omega_entry, np_iter.q_min, q) - 0.90, [np_iter.q_min, np_iter.q_max]);

                        
% omega_interp_l = @(q) interp1(np_iter.q, eqm_iter.omega(1,:), q, 'spline', 'extrap');
% omega_interp_h = @(q) interp1(np_iter.q, eqm_iter.omega(2,:), q, 'spline', 'extrap');
% omega_entry = @(q) 1./(q * p_iter.tau * sqrt(2*pi)) .* exp(-(log(q) - p_iter.iota).^2 / (2 * p_iter.tau^2));

% result_Expected_q(i,j, l,f,g, 1) = integral(@(q) q .* omega_interp_l(q), np_iter.q_min, np_iter.q_max);  % Choose appropriate limits
% result_10quantile_q(i,j, l,f,g, 1) = fzero(@(q) integral(omega_interp_l, np_iter.q_min, q) - 0.10, [np_iter.q_min, np_iter.q_max]);
% result_50quantile_q(i,j, l,f,g, 1) = fzero(@(q) integral(omega_interp_l, np_iter.q_min, q) - 0.50, [np_iter.q_min, np_iter.q_max]);
% result_90quantile_q(i,j, l,f,g, 1) = fzero(@(q) integral(omega_interp_l, np_iter.q_min, q) - 0.90, [np_iter.q_min, np_iter.q_max]);


% result_Expected_q(i,j, l,f,g, 2) = integral(@(q) q .* omega_interp_h(q), np_iter.q_min, np_iter.q_max);  % Choose appropriate limits
% result_10quantile_q(i,j, l,f,g, 2) = fzero(@(q) integral(omega_interp_h, np_iter.q_min, q) - 0.10, [np_iter.q_min, np_iter.q_max]);
% result_50quantile_q(i,j, l,f,g, 2) = fzero(@(q) integral(omega_interp_h, np_iter.q_min, q) - 0.50, [np_iter.q_min, np_iter.q_max]);
% result_90quantile_q(i,j, l,f,g, 2) = fzero(@(q) integral(omega_interp_h, np_iter.q_min, q) - 0.90, [np_iter.q_min, np_iter.q_max]);



% result_Expected_q(i,j, l,f,g, 3) = integral(@(q) q .* omega_entry(q), np_iter.q_min, np_iter.q_max);  % Choose appropriate limits
% result_10quantile_q(i,j, l,f,g, 3) = fzero(@(q) integral(omega_entry, np_iter.q_min, q) - 0.10, [np_iter.q_min, np_iter.q_max]);
% result_50quantile_q(i,j, l,f,g, 3) = fzero(@(q) integral(omega_entry, np_iter.q_min, q) - 0.50, [np_iter.q_min, np_iter.q_max]);
% result_90quantile_q(i,j, l,f,g, 3) = fzero(@(q) integral(omega_entry, np_iter.q_min, q) - 0.90, [np_iter.q_min, np_iter.q_max]);





% result_s(i,j, l,f,g, :) = eqm_iter.s;
% result_g_3(i,j, l,f,g, 1) = (eqm_iter.gp1_type(1))^(1/p_iter.sigma-1);
% result_g_3(i,j, l,f,g, 2) = (eqm_iter.gp1_type(2))^(1/p_iter.sigma-1);
% result_g_3(i,j, l,f,g, 3) = eqm_iter.g;

% result_mu_hat(i,j, l,f,g, :, :) = eqm_iter.mu_hat;
% result_nu_hat(i,j, l,f,g, :) = eqm_iter.nu_hat;
% result_delta_hat(i,j, l,f,g, :) = eqm_iter.delta_hat;

% if iter<500
%     result_flag(i,j, l,f,g, 1) = 1;
% else
%     result_flag(i,j, l,f,g, 1) = 0;
% end



%%
color_set = {'r', 'g', 'b', 'k', 'm', 'c'};

for i = 1:length(gamma)

    figure;

    subplot(1, 2, 1);
    hold on;
    legend_handles = []; % Initialize an array to store scatter handles for the legend
    legend_labels = {};  % Initialize legend labels

    for g = 1:length(alpha)

        color = color_set{g};
        legend_labels{g} = sprintf('$\\alpha = %.2f$', alpha(g)); % Add legend label for each alpha

        % Create a dummy scatter plot for the legend
        h = scatter(NaN, NaN, color, 'filled');
        legend_handles = [legend_handles, h];

        for l = 1:length(r)
            for f = 1:length(theta)
                for j = 1:length(m)
                    if result_flag(i,j, l,f,g, 1) == 1
                        scatter(m(j), result_s(i,j, l,f,g, 1), color, 'filled');
                    elseif result_flag(i,j, l,f,g, 1) == 0
                        scatter(m(j), result_s(i,j, l,f,g, 1), 'x', 'MarkerEdgeColor', color);
                    end
                end
            end
        end

    end

    xlabel('$\bar{m}$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('$s_l$', 'Interpreter', 'latex', 'FontSize', 20);
    ylim([0, 1]);
    yline(0.5, 'k--', 'LineWidth', 1.5);
    grid on;
    legend(legend_handles, legend_labels, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best'); % Set legend font size
    title('Market Share For Type L', 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
    hold off;


    subplot(1, 2, 2);
    hold on;
    legend_handles = []; % Initialize an array to store scatter handles for the legend
    legend_labels = {};  % Initialize legend labels

    for g = 1:length(alpha)

        color = color_set{g};
        legend_labels{g} = sprintf('$\\alpha = %.2f$', alpha(g)); % Add legend label for each alpha

        % Create a dummy scatter plot for the legend
        h = scatter(NaN, NaN, color, 'filled');
        legend_handles = [legend_handles, h];

        for l = 1:length(r)
            for f = 1:length(theta)
                for j = 1:length(m)
                    if result_flag(i,j, l,f,g, 1) == 1
                        scatter(m(j), result_s(i,j, l,f,g, 2), color, 'filled');
                    elseif result_flag(i,j, l,f,g, 1) == 0
                        scatter(m(j), result_s(i,j, l,f,g, 2), 'x', 'MarkerEdgeColor', color);
                    end
                end
            end
        end

    end

    xlabel('$\bar{m}$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('$s_h$', 'Interpreter', 'latex', 'FontSize', 20);
    ylim([0, 1]);
    grid on;
    yline(0.5, 'k--', 'LineWidth', 1.5);
    legend(legend_handles, legend_labels, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best'); % Set legend font size
    title('Market Share For Type H', 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
    hold off;

    sgtitle(sprintf('$\\gamma = %.2f$', gamma(i)), 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
    set(gcf, 'Position', [100, 100, 2200, 1200]); % Set figure size
    saveas(gcf, sprintf ('./figure_analysis/MarketShare_gamma=%.2f.png', gamma(i)));
end


%%


%%
color_set = {'r', 'g', 'b', 'k', 'm', 'c'};

for i = 1:length(gamma)

    figure;

    subplot(1, 3, 1);
    
    legend_handles = []; % Initialize an array to store scatter handles for the legend
    legend_labels = {};  % Initialize legend labels

    hold on;
    for g = 1:length(alpha)

        color = color_set{g};
        legend_labels{g} = sprintf('$\\alpha = %.2f$', alpha(g)); % Add legend label for each alpha

        % Create a dummy scatter plot for the legend
        h = scatter(NaN, NaN, color, 'filled');
        legend_handles = [legend_handles, h];

        for l = 1:length(r)
            for f = 1:length(theta)
                for j = 1:length(m)
                    if result_flag(i,j, l,f,g, 1) == 1
                        scatter(m(j), result_g_3(i,j, l,f,g, 1), color, 'filled');
                    elseif result_flag(i,j, l,f,g, 1) == 0
                        scatter(m(j), result_g_3(i,j, l,f,g, 1), 'x', 'MarkerEdgeColor', color);
                    end
                end
            end
        end

    end

    xlabel('$\bar{m}$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('$g_l$', 'Interpreter', 'latex', 'FontSize', 20);
    ylim([-0.35, 0.4]);
    grid on;
    yline(0, 'k--', 'LineWidth', 1.5);
    legend(legend_handles, legend_labels, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best'); % Set legend font size
    title('Growth Rate For Type L', 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
    hold off;


    subplot(1, 3, 2);
    hold on;
    
    legend_handles = []; % Initialize an array to store scatter handles for the legend
    legend_labels = {};  % Initialize legend labels


    for g = 1:length(alpha)

        color = color_set{g};
        legend_labels{g} = sprintf('$\\alpha = %.2f$', alpha(g)); % Add legend label for each alpha

        % Create a dummy scatter plot for the legend
        h = scatter(NaN, NaN, color, 'filled');
        legend_handles = [legend_handles, h];

        for l = 1:length(r)
            for f = 1:length(theta)
                for j = 1:length(m)
                    if result_flag(i,j, l,f,g, 1) == 1
                        scatter(m(j), result_g_3(i,j, l,f,g, 2), color, 'filled');
                    elseif result_flag(i,j, l,f,g, 1) == 0
                        scatter(m(j), result_g_3(i,j, l,f,g, 2), 'x', 'MarkerEdgeColor', color);
                    end
                end
            end
        end

    end

    xlabel('$\bar{m}$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('$g_h$', 'Interpreter', 'latex', 'FontSize', 20);
    ylim([-0.35, 0.4]);
    grid on;
    yline(0, 'k--', 'LineWidth', 1.5);
    legend(legend_handles, legend_labels, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best'); % Set legend font size
    title('Growth Rate For Type H', 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
    hold off;





    subplot(1, 3, 3);
    hold on;
    legend_handles = []; % Initialize an array to store scatter handles for the legend
    legend_labels = {};  % Initialize legend labels


    for g = 1:length(alpha)

        color = color_set{g};
        legend_labels{g} = sprintf('$\\alpha = %.2f$', alpha(g)); % Add legend label for each alpha

        % Create a dummy scatter plot for the legend
        h = scatter(NaN, NaN, color, 'filled');
        legend_handles = [legend_handles, h];

        for l = 1:length(r)
            for f = 1:length(theta)
                for j = 1:length(m)
                    if result_flag(i,j, l,f,g, 1) == 1
                        scatter(m(j), result_g_3(i,j, l,f,g, 3), color, 'filled');
                    elseif result_flag(i,j, l,f,g, 1) == 0
                        scatter(m(j), result_g_3(i,j, l,f,g, 3), 'x', 'MarkerEdgeColor', color);
                    end
                end
            end
        end

    end

    xlabel('$\bar{m}$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('$g$', 'Interpreter', 'latex', 'FontSize', 20);
    ylim([-0.35, 0.4]);
    grid on;
    yline(0, 'k--', 'LineWidth', 1.5);
    legend(legend_handles, legend_labels, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best'); % Set legend font size
    title('Aggregate Growth Rate', 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
    hold off;


    sgtitle(sprintf('$\\gamma = %.2f$', gamma(i)), 'Interpreter', 'latex', 'FontSize', 20); % Set title font size



    set(gcf, 'Position', [100, 100, 2200, 1200]); % Set figure size
    saveas(gcf, sprintf ('./figure_analysis/GrowthRate_gamma=%.2f.png', gamma(i)));

end



%%
color_set = {'r', 'g', 'b', 'k', 'm', 'c'};

for i = 1:length(gamma)

    figure;

    subplot(1, 2, 1);
    
    legend_handles = []; % Initialize an array to store scatter handles for the legend
    legend_labels = {};  % Initialize legend labels

    hold on;
    for g = 1:length(alpha)

        color = color_set{g};
        legend_labels{g} = sprintf('$\\alpha = %.2f$', alpha(g)); % Add legend label for each alpha

        % Create a dummy scatter plot for the legend
        h = scatter(NaN, NaN, color, 'filled');
        legend_handles = [legend_handles, h];

        for l = 1:length(r)
            for f = 1:length(theta)
                for j = 1:length(m)
                    if result_flag(i,j, l,f,g, 1) == 1
                        scatter(m(j), result_mu_hat(i,j, l,f,g, 1,1) + result_mu_hat(i,j, l,f,g, 1,2), color, 'filled');
                    elseif result_flag(i,j, l,f,g, 1) == 0
                        scatter(m(j), result_mu_hat(i,j, l,f,g, 1,1) + result_mu_hat(i,j, l,f,g, 1,2), 'x', 'MarkerEdgeColor', color);
                    end
                end
            end
        end

    end

    xlabel('$\bar{m}$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('$\hat{\mu}_{11} + \hat{\mu}_{12}$', 'Interpreter', 'latex', 'FontSize', 20);
    ylim([0, 1]);
    grid on;
    % yline(0, 'k--', 'LineWidth', 1.5);
    legend(legend_handles, legend_labels, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best'); % Set legend font size
    title('Type L Steal Prob $\hat{\mu}_{11} + \hat{\mu}_{12}$', 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
    hold off;


    subplot(1, 2, 2);
    hold on;
    
    legend_handles = []; % Initialize an array to store scatter handles for the legend
    legend_labels = {};  % Initialize legend labels


    for g = 1:length(alpha)

        color = color_set{g};
        legend_labels{g} = sprintf('$\\alpha = %.2f$', alpha(g)); % Add legend label for each alpha

        % Create a dummy scatter plot for the legend
        h = scatter(NaN, NaN, color, 'filled');
        legend_handles = [legend_handles, h];

        for l = 1:length(r)
            for f = 1:length(theta)
                for j = 1:length(m)
                    if result_flag(i,j, l,f,g, 1) == 1
                        scatter(m(j), result_mu_hat(i,j, l,f,g, 2,1) + result_mu_hat(i,j, l,f,g, 2,2), color, 'filled');
                    elseif result_flag(i,j, l,f,g, 1) == 0
                        scatter(m(j), result_mu_hat(i,j, l,f,g, 2,1) + result_mu_hat(i,j, l,f,g, 2,2), 'x', 'MarkerEdgeColor', color);
                    end
                end
            end
        end

    end

    xlabel('$\bar{m}$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel(' $\hat{\mu}_{21} + \hat{\mu}_{22}$', 'Interpreter', 'latex', 'FontSize', 20);
    ylim([0, 1]);
    grid on;
    % yline(0, 'k--', 'LineWidth', 1.5);
    legend(legend_handles, legend_labels, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best'); % Set legend font size
    title('Type H Steal Prob $\hat{\mu}_{21} + \hat{\mu}_{22}$', 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
    hold off;





    sgtitle(sprintf('$\\gamma = %.2f$', gamma(i)), 'Interpreter', 'latex', 'FontSize', 20); % Set title font size



    set(gcf, 'Position', [100, 100, 2200, 1200]); % Set figure size
    saveas(gcf, sprintf ('./figure_analysis/StealProbability_gamma=%.2f.png', gamma(i)));

end


%%




%%
color_set = {'r', 'g', 'b', 'k', 'm', 'c'};

for i = 1:length(gamma)

    figure;

    subplot(1, 2, 1);
    
    legend_handles = []; % Initialize an array to store scatter handles for the legend
    legend_labels = {};  % Initialize legend labels

    hold on;
    for g = 1:length(alpha)

        color = color_set{g};
        legend_labels{g} = sprintf('$\\alpha = %.2f$', alpha(g)); % Add legend label for each alpha

        % Create a dummy scatter plot for the legend
        h = scatter(NaN, NaN, color, 'filled');
        legend_handles = [legend_handles, h];

        for l = 1:length(r)
            for f = 1:length(theta)
                for j = 1:length(m)
                    if result_flag(i,j, l,f,g, 1) == 1
                        scatter(m(j), result_nu_hat(i,j, l,f,g, 1) , color, 'filled');
                    elseif result_flag(i,j, l,f,g, 1) == 0
                        scatter(m(j), result_nu_hat(i,j, l,f,g, 1) , 'x', 'MarkerEdgeColor', color);
                    end
                end
            end
        end


    end

    xlabel('$\bar{m}$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('  $\hat{\nu}_{1}$', 'Interpreter', 'latex', 'FontSize', 20);
    ylim([0, 1]);
    grid on;
    % yline(0, 'k--', 'LineWidth', 1.5);
    legend(legend_handles, legend_labels, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best'); % Set legend font size
    title('Type L Save Prob $\hat{\nu}_{1}$', 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
    hold off;


    subplot(1, 2, 2);
    hold on;
    
    legend_handles = []; % Initialize an array to store scatter handles for the legend
    legend_labels = {};  % Initialize legend labels


    for g = 1:length(alpha)

        color = color_set{g};
        legend_labels{g} = sprintf('$\\alpha = %.2f$', alpha(g)); % Add legend label for each alpha

        % Create a dummy scatter plot for the legend
        h = scatter(NaN, NaN, color, 'filled');
        legend_handles = [legend_handles, h];

        for l = 1:length(r)
            for f = 1:length(theta)
                for j = 1:length(m)
                    if result_flag(i,j, l,f,g, 1) == 1
                        scatter(m(j), result_nu_hat(i,j, l,f,g, 2) , color, 'filled');
                    elseif result_flag(i,j, l,f,g, 1) == 0
                        scatter(m(j), result_nu_hat(i,j, l,f,g, 2) , 'x', 'MarkerEdgeColor', color);
                    end
                end
            end
        end

    end

    xlabel('$\bar{m}$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('  $\hat{\nu}_{2}$', 'Interpreter', 'latex', 'FontSize', 20);
    ylim([0, 1]);
    grid on;
    % yline(0, 'k--', 'LineWidth', 1.5);
    legend(legend_handles, legend_labels, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best'); % Set legend font size
    title('Type H Save Prob $\hat{\nu}_{2}$', 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
    hold off;





    sgtitle(sprintf('$\\gamma = %.2f$', gamma(i)), 'Interpreter', 'latex', 'FontSize', 20); % Set title font size



    set(gcf, 'Position', [100, 100, 2200, 1200]); % Set figure size
    saveas(gcf, sprintf ('./figure_analysis/SaveProbability_gamma=%.2f.png', gamma(i)));

end





%%
color_set = {'r', 'g', 'b', 'k', 'm', 'c'};

for i = 1:length(gamma)

    figure;

    subplot(1, 2, 1);
    
    legend_handles = []; % Initialize an array to store scatter handles for the legend
    legend_labels = {};  % Initialize legend labels

    hold on;
    for g = 1:length(alpha)

        color = color_set{g};
        legend_labels{g} = sprintf('$\\alpha = %.2f$', alpha(g)); % Add legend label for each alpha

        % Create a dummy scatter plot for the legend
        h = scatter(NaN, NaN, color, 'filled');
        legend_handles = [legend_handles, h];

        for l = 1:length(r)
            for f = 1:length(theta)
                for j = 1:length(m)
                    if result_flag(i,j, l,f,g, 1) == 1
                        scatter(m(j), result_delta_hat(i,j, l,f,g, 1) , color, 'filled');
                    elseif result_flag(i,j, l,f,g, 1) == 0
                        scatter(m(j), result_delta_hat(i,j, l,f,g, 1) , 'x', 'MarkerEdgeColor', color);
                    end
                end
            end
        end


    end

    xlabel('$\bar{m}$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('  $\hat{\delta}_{1}$', 'Interpreter', 'latex', 'FontSize', 20);
    ylim([0, 1]);
    grid on;
    % yline(0, 'k--', 'LineWidth', 1.5);
    legend(legend_handles, legend_labels, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best'); % Set legend font size
    title('Type H Expected External Innovation Probability $\hat{\delta}_{1}$', 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
    hold off;


    subplot(1, 2, 2);
    hold on;
    
    legend_handles = []; % Initialize an array to store scatter handles for the legend
    legend_labels = {};  % Initialize legend labels


    for g = 1:length(alpha)

        color = color_set{g};
        legend_labels{g} = sprintf('$\\alpha = %.2f$', alpha(g)); % Add legend label for each alpha

        % Create a dummy scatter plot for the legend
        h = scatter(NaN, NaN, color, 'filled');
        legend_handles = [legend_handles, h];

        for l = 1:length(r)
            for f = 1:length(theta)
                for j = 1:length(m)
                    if result_flag(i,j, l,f,g, 1) == 1
                        scatter(m(j), result_delta_hat(i,j, l,f,g, 2) , color, 'filled');
                    elseif result_flag(i,j, l,f,g, 1) == 0
                        scatter(m(j), result_delta_hat(i,j, l,f,g, 2) , 'x', 'MarkerEdgeColor', color);
                    end
                end
            end
        end

    end

    xlabel('$\bar{m}$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('  $\hat{\delta}_{2}$', 'Interpreter', 'latex', 'FontSize', 20);
    ylim([0, 1]);
    grid on;
    % yline(0, 'k--', 'LineWidth', 1.5);
    legend(legend_handles, legend_labels, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best'); % Set legend font size
    title('Type H Expected External Innovation Probability $\hat{\delta}_{2}$', 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
    hold off;





    sgtitle(sprintf('$\\gamma = %.2f$', gamma(i)), 'Interpreter', 'latex', 'FontSize', 20); % Set title font size



    set(gcf, 'Position', [100, 100, 2200, 1200]); % Set figure size
    saveas(gcf, sprintf ('./figure_analysis/Expected_External_gamma=%.2f.png', gamma(i)));

end




%%

%%
color_set = {'r', 'g', 'b', 'k', 'm', 'c'};

for i = 1:length(gamma)

    figure;

    subplot(1, 2, 1);
    
    legend_handles = []; % Initialize an array to store scatter handles for the legend
    legend_labels = {};  % Initialize legend labels

    hold on;
    for g = 1:length(alpha)

        color = color_set{g};
        legend_labels{g} = sprintf('$\\alpha = %.2f$', alpha(g)); % Add legend label for each alpha

        % Create a dummy scatter plot for the legend
        h = scatter(NaN, NaN, color, 'filled');
        legend_handles = [legend_handles, h];

        for l = 1:length(r)
            for f = 1:length(theta)
                for j = 1:length(m)
                    if result_flag(i,j, l,f,g, 1) == 1
                        scatter(m(j), result_lambda_hat(i,j, l,f,g, 1) , color, 'filled');
                    elseif result_flag(i,j, l,f,g, 1) == 0
                        scatter(m(j), result_lambda_hat(i,j, l,f,g, 1) , 'x', 'MarkerEdgeColor', color);
                    end
                end
            end
        end


    end

    xlabel('$\bar{m}$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('  $\hat{\lambda}_{1}$', 'Interpreter', 'latex', 'FontSize', 20);
    ylim([0, 1]);
    grid on;
    % yline(0, 'k--', 'LineWidth', 1.5);
    legend(legend_handles, legend_labels, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best'); % Set legend font size
    title('Type H Expected Internal Innovation Probability $\hat{\lambda}_{1}$', 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
    hold off;


    subplot(1, 2, 2);
    hold on;
    
    legend_handles = []; % Initialize an array to store scatter handles for the legend
    legend_labels = {};  % Initialize legend labels


    for g = 1:length(alpha)

        color = color_set{g};
        legend_labels{g} = sprintf('$\\alpha = %.2f$', alpha(g)); % Add legend label for each alpha

        % Create a dummy scatter plot for the legend
        h = scatter(NaN, NaN, color, 'filled');
        legend_handles = [legend_handles, h];

        for l = 1:length(r)
            for f = 1:length(theta)
                for j = 1:length(m)
                    if result_flag(i,j, l,f,g, 1) == 1
                        scatter(m(j), result_lambda_hat(i,j, l,f,g, 2) , color, 'filled');
                    elseif result_flag(i,j, l,f,g, 1) == 0
                        scatter(m(j), result_lambda_hat(i,j, l,f,g, 2) , 'x', 'MarkerEdgeColor', color);
                    end
                end
            end
        end

    end

    xlabel('$\bar{m}$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('  $\hat{\lambda}_{2}$', 'Interpreter', 'latex', 'FontSize', 20);
    ylim([0, 1]);
    grid on;
    % yline(0, 'k--', 'LineWidth', 1.5);
    legend(legend_handles, legend_labels, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best'); % Set legend font size
    title('Type H Expected Internal Innovation Probability $\hat{\lambda}_{2}$', 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
    hold off;





    sgtitle(sprintf('$\\gamma = %.2f$', gamma(i)), 'Interpreter', 'latex', 'FontSize', 20); % Set title font size



    set(gcf, 'Position', [100, 100, 2200, 1200]); % Set figure size
    saveas(gcf, sprintf ('./figure_analysis/Expected_Internal_gamma=%.2f.png', gamma(i)));

end



%%


%%
color_set = {'r', 'g', 'b', 'k', 'm', 'c'};

for i = 1:length(gamma)

    figure;

    for g = 1:length(alpha)

        subplot(1, length(alpha), g);

        legend_handles = []; % Initialize an array to store scatter handles for the legend
        legend_labels = {};  % Initialize legend labels

        hold on;

        color = color_set{1};        
        h1 = scatter(NaN, NaN, color, 'filled'); % Dummy scatter for legend
        legend_handles = [legend_handles, h1];
        legend_labels{end+1} = '$E(q_l)$';

        % legend_labels{1} = sprintf('$E(q_l)$'); % Add legend label for each alpha

        % Create a dummy scatter plot for the legend
        % legend_handles = [legend_handles, h];

        for l = 1:length(r)
            for f = 1:length(theta)
                for j = 1:length(m)
                    if result_flag(i,j, l,f,g, 1) == 1
                        scatter(m(j), result_Expected_q(i,j, l,f,g, 1) , color, 'filled');
                    elseif result_flag(i,j, l,f,g, 1) == 0
                        scatter(m(j), result_Expected_q(i,j, l,f,g, 1) , 'x', 'MarkerEdgeColor', color);
                    end
                end
            end
        end


        color = color_set{2};

        h2 = scatter(NaN, NaN, color, 'filled'); % Dummy scatter for legend
        legend_handles = [legend_handles, h2];
        legend_labels{end+1} = '$E(q_h)$';
        % Create a dummy scatter plot for the legend
        % legend_handles = [legend_handles, h];

        for l = 1:length(r)
            for f = 1:length(theta)
                for j = 1:length(m)
                    if result_flag(i,j, l,f,g, 1) == 1
                        scatter(m(j), result_Expected_q(i,j, l,f,g, 2) , color, 'filled');
                    elseif result_flag(i,j, l,f,g, 1) == 0
                        scatter(m(j), result_Expected_q(i,j, l,f,g, 2) , 'x', 'MarkerEdgeColor', color);
                    end
                end
            end
        end



        color = color_set{3};

        h3 = scatter(NaN, NaN, color, 'filled'); % Dummy scatter for legend
        legend_handles = [legend_handles, h3];
        legend_labels{end+1} = '$E(q_{entry})$';

        % Create a dummy scatter plot for the legend
        % legend_handles = [legend_handles, h];

        for l = 1:length(r)
            for f = 1:length(theta)
                for j = 1:length(m)
                    if result_flag(i,j, l,f,g, 1) == 1
                        scatter(m(j), result_Expected_q(i,j, l,f,g, 3) , color, 'filled');
                    elseif result_flag(i,j, l,f,g, 1) == 0
                        scatter(m(j), result_Expected_q(i,j, l,f,g, 3) , 'x', 'MarkerEdgeColor', color);
                    end
                end
            end
        end

        xlabel('$\bar{m}$', 'Interpreter', 'latex', 'FontSize', 20);
        ylabel('  $E(q)$', 'Interpreter', 'latex', 'FontSize', 20);
        ylim([0.8, 3.6]);
        grid on;
        % yline(0, 'k--', 'LineWidth', 1.5);
        legend(legend_handles, legend_labels, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best'); % Set legend font size
        title(sprintf(' $\\alpha = %.2f$', alpha(g)), 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
        set(gcf, 'Position', [100, 100, 2600, 1200]); % Set figure size
        hold off;



    end

    sgtitle(sprintf('$\\gamma = %.2f$', gamma(i)), 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
    saveas(gcf, sprintf ('./figure_analysis/Expected_q_gamma=%.2f.png', gamma(i)));
end




%%


%%
color_set = {'r', 'g', 'b', 'k', 'm', 'c'};

for i = 1:length(gamma)

    figure;

    for g = 1:length(alpha)

        subplot(1, length(alpha), g);

        legend_handles = []; % Initialize an array to store scatter handles for the legend
        legend_labels = {};  % Initialize legend labels

        hold on;

        color = color_set{1};        
        h1 = scatter(NaN, NaN, color, 'filled'); % Dummy scatter for legend
        legend_handles = [legend_handles, h1];
        legend_labels{end+1} = '$q_{l}(0.1)$';

        % legend_labels{1} = sprintf('$E(q_l)$'); % Add legend label for each alpha

        % Create a dummy scatter plot for the legend
        % legend_handles = [legend_handles, h];

        for l = 1:length(r)
            for f = 1:length(theta)
                for j = 1:length(m)
                    if result_flag(i,j, l,f,g, 1) == 1
                        scatter(m(j), result_10quantile_q(i,j, l,f,g, 1) , color, 'filled');
                    elseif result_flag(i,j, l,f,g, 1) == 0
                        scatter(m(j), result_10quantile_q(i,j, l,f,g, 1) , 'x', 'MarkerEdgeColor', color);
                    end
                end
            end
        end


        color = color_set{2};

        h2 = scatter(NaN, NaN, color, 'filled'); % Dummy scatter for legend
        legend_handles = [legend_handles, h2];
        legend_labels{end+1} = '$q_{h}(0.1)$';
        % Create a dummy scatter plot for the legend
        % legend_handles = [legend_handles, h];

        for l = 1:length(r)
            for f = 1:length(theta)
                for j = 1:length(m)
                    if result_flag(i,j, l,f,g, 1) == 1
                        scatter(m(j), result_10quantile_q(i,j, l,f,g, 2) , color, 'filled');
                    elseif result_flag(i,j, l,f,g, 1) == 0
                        scatter(m(j), result_10quantile_q(i,j, l,f,g, 2) , 'x', 'MarkerEdgeColor', color);
                    end
                end
            end
        end



        color = color_set{3};

        h3 = scatter(NaN, NaN, color, 'filled'); % Dummy scatter for legend
        legend_handles = [legend_handles, h3];
        legend_labels{end+1} = '$q_{entry}(0.1)$';

        % Create a dummy scatter plot for the legend
        % legend_handles = [legend_handles, h];

        for l = 1:length(r)
            for f = 1:length(theta)
                for j = 1:length(m)
                    if result_flag(i,j, l,f,g, 1) == 1
                        scatter(m(j), result_10quantile_q(i,j, l,f,g, 3) , color, 'filled');
                    elseif result_flag(i,j, l,f,g, 1) == 0
                        scatter(m(j), result_10quantile_q(i,j, l,f,g, 3) , 'x', 'MarkerEdgeColor', color);
                    end
                end
            end
        end

        xlabel('$\bar{m}$', 'Interpreter', 'latex', 'FontSize', 20);
        ylabel('  $q_{10\%}$', 'Interpreter', 'latex', 'FontSize', 20);
        ylim([0.3, 0.6]);
        grid on;
        % yline(0, 'k--', 'LineWidth', 1.5);
        legend(legend_handles, legend_labels, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best'); % Set legend font size
        title(sprintf(' $\\alpha = %.2f$', alpha(g)), 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
        set(gcf, 'Position', [100, 100, 2600, 1200]); % Set figure size
        hold off;



    end

    sgtitle(sprintf('$\\gamma = %.2f$', gamma(i)), 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
    saveas(gcf, sprintf ('./figure_analysis/10quantile_gamma=%.2f.png', gamma(i)));
end





%%
color_set = {'r', 'g', 'b', 'k', 'm', 'c'};

for i = 1:length(gamma)

    figure;

    for g = 1:length(alpha)

        subplot(1, length(alpha), g);

        legend_handles = []; % Initialize an array to store scatter handles for the legend
        legend_labels = {};  % Initialize legend labels

        hold on;

        color = color_set{1};        
        h1 = scatter(NaN, NaN, color, 'filled'); % Dummy scatter for legend
        legend_handles = [legend_handles, h1];
        legend_labels{end+1} = '$q_{l}(0.5)$';

        % legend_labels{1} = sprintf('$E(q_l)$'); % Add legend label for each alpha

        % Create a dummy scatter plot for the legend
        % legend_handles = [legend_handles, h];

        for l = 1:length(r)
            for f = 1:length(theta)
                for j = 1:length(m)
                    if result_flag(i,j, l,f,g, 1) == 1
                        scatter(m(j), result_50quantile_q(i,j, l,f,g, 1) , color, 'filled');
                    elseif result_flag(i,j, l,f,g, 1) == 0
                        scatter(m(j), result_50quantile_q(i,j, l,f,g, 1) , 'x', 'MarkerEdgeColor', color);
                    end
                end
            end
        end


        color = color_set{2};

        h2 = scatter(NaN, NaN, color, 'filled'); % Dummy scatter for legend
        legend_handles = [legend_handles, h2];
        legend_labels{end+1} = '$q_{h}(0.5)$';
        % Create a dummy scatter plot for the legend
        % legend_handles = [legend_handles, h];

        for l = 1:length(r)
            for f = 1:length(theta)
                for j = 1:length(m)
                    if result_flag(i,j, l,f,g, 1) == 1
                        scatter(m(j), result_50quantile_q(i,j, l,f,g, 2) , color, 'filled');
                    elseif result_flag(i,j, l,f,g, 1) == 0
                        scatter(m(j), result_50quantile_q(i,j, l,f,g, 2) , 'x', 'MarkerEdgeColor', color);
                    end
                end
            end
        end



        color = color_set{3};

        h3 = scatter(NaN, NaN, color, 'filled'); % Dummy scatter for legend
        legend_handles = [legend_handles, h3];
        legend_labels{end+1} = '$q_{entry}(0.5)$';

        % Create a dummy scatter plot for the legend
        % legend_handles = [legend_handles, h];

        for l = 1:length(r)
            for f = 1:length(theta)
                for j = 1:length(m)
                    if result_flag(i,j, l,f,g, 1) == 1
                        scatter(m(j), result_50quantile_q(i,j, l,f,g, 3) , color, 'filled');
                    elseif result_flag(i,j, l,f,g, 1) == 0
                        scatter(m(j), result_50quantile_q(i,j, l,f,g, 3) , 'x', 'MarkerEdgeColor', color);
                    end
                end
            end
        end

        xlabel('$\bar{m}$', 'Interpreter', 'latex', 'FontSize', 20);
        ylabel('  $q_{50\%}$', 'Interpreter', 'latex', 'FontSize', 20);
        ylim([0.7, 2]);
        grid on;
        % yline(0, 'k--', 'LineWidth', 1.5);
        legend(legend_handles, legend_labels, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best'); % Set legend font size
        title(sprintf(' $\\alpha = %.2f$', alpha(g)), 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
        set(gcf, 'Position', [100, 100, 2600, 1200]); % Set figure size
        hold off;



    end

    sgtitle(sprintf('$\\gamma = %.2f$', gamma(i)), 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
    saveas(gcf, sprintf ('./figure_analysis/50quantile_gamma=%.2f.png', gamma(i)));
end




%%
color_set = {'r', 'g', 'b', 'k', 'm', 'c'};

for i = 1:length(gamma)

    figure;

    for g = 1:length(alpha)

        subplot(1, length(alpha), g);

        legend_handles = []; % Initialize an array to store scatter handles for the legend
        legend_labels = {};  % Initialize legend labels

        hold on;

        color = color_set{1};        
        h1 = scatter(NaN, NaN, color, 'filled'); % Dummy scatter for legend
        legend_handles = [legend_handles, h1];
        legend_labels{end+1} = '$q_{l}(0.9)$';

        % legend_labels{1} = sprintf('$E(q_l)$'); % Add legend label for each alpha

        % Create a dummy scatter plot for the legend
        % legend_handles = [legend_handles, h];

        for l = 1:length(r)
            for f = 1:length(theta)
                for j = 1:length(m)
                    if result_flag(i,j, l,f,g, 1) == 1
                        scatter(m(j), result_90quantile_q(i,j, l,f,g, 1) , color, 'filled');
                    elseif result_flag(i,j, l,f,g, 1) == 0
                        scatter(m(j), result_90quantile_q(i,j, l,f,g, 1) , 'x', 'MarkerEdgeColor', color);
                    end
                end
            end
        end


        color = color_set{2};

        h2 = scatter(NaN, NaN, color, 'filled'); % Dummy scatter for legend
        legend_handles = [legend_handles, h2];
        legend_labels{end+1} = '$q_{h}(0.9)$';
        % Create a dummy scatter plot for the legend
        % legend_handles = [legend_handles, h];

        for l = 1:length(r)
            for f = 1:length(theta)
                for j = 1:length(m)
                    if result_flag(i,j, l,f,g, 1) == 1
                        scatter(m(j), result_90quantile_q(i,j, l,f,g, 2) , color, 'filled');
                    elseif result_flag(i,j, l,f,g, 1) == 0
                        scatter(m(j), result_90quantile_q(i,j, l,f,g, 2) , 'x', 'MarkerEdgeColor', color);
                    end
                end
            end
        end



        color = color_set{3};

        h3 = scatter(NaN, NaN, color, 'filled'); % Dummy scatter for legend
        legend_handles = [legend_handles, h3];
        legend_labels{end+1} = '$q_{entry}(0.9)$';

        % Create a dummy scatter plot for the legend
        % legend_handles = [legend_handles, h];

        for l = 1:length(r)
            for f = 1:length(theta)
                for j = 1:length(m)
                    if result_flag(i,j, l,f,g, 1) == 1
                        scatter(m(j), result_90quantile_q(i,j, l,f,g, 3) , color, 'filled');
                    elseif result_flag(i,j, l,f,g, 1) == 0
                        scatter(m(j), result_90quantile_q(i,j, l,f,g, 3) , 'x', 'MarkerEdgeColor', color);
                    end
                end
            end
        end

        xlabel('$\bar{m}$', 'Interpreter', 'latex', 'FontSize', 20);
        ylabel('  $q_{90\%}$', 'Interpreter', 'latex', 'FontSize', 20);
        ylim([1, 10]);
        grid on;
        % yline(0, 'k--', 'LineWidth', 1.5);
        legend(legend_handles, legend_labels, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best'); % Set legend font size
        title(sprintf(' $\\alpha = %.2f$', alpha(g)), 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
        set(gcf, 'Position', [100, 100, 2600, 1200]); % Set figure size
        hold off;



    end

    sgtitle(sprintf('$\\gamma = %.2f$', gamma(i)), 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
    saveas(gcf, sprintf ('./figure_analysis/90quantile_gamma=%.2f.png', gamma(i)));
end


