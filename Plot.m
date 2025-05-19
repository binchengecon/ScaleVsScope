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
epsilon = [0.0, 1.0];
beta = [2 8];
result_Expected_q = zeros(length(gamma), length(m), legnth(phi_h), length(r), length(theta), length(alpha), length(epsilon), length(beta), 3);
result_10quantile_q = zeros(length(gamma), length(m), legnth(phi_h), length(r), length(theta), length(alpha), length(epsilon), length(beta), 3);
result_50quantile_q = zeros(length(gamma), length(m), legnth(phi_h), length(r), length(theta), length(alpha), length(epsilon), length(beta), 3);
result_90quantile_q = zeros(length(gamma), length(m), legnth(phi_h), length(r), length(theta), length(alpha), length(epsilon), length(beta), 3);
result_s = zeros(length(gamma), length(m), length(r), legnth(phi_h), length(theta), length(alpha), length(epsilon), length(beta), 2);
result_g_3 = zeros(length(gamma), length(m), length(r), legnth(phi_h), length(theta), length(alpha), length(epsilon), length(beta), 3);
result_mu_hat = zeros(length(gamma), length(m), length(r), legnth(phi_h), length(theta), length(alpha), length(epsilon), length(beta), 2,2);
result_nu_hat = zeros(length(gamma), length(m), length(r), legnth(phi_h), length(theta), length(alpha), length(epsilon), length(beta), 2);
result_delta_hat = zeros(length(gamma), length(m), length(r), legnth(phi_h), length(theta), length(alpha), length(epsilon), length(beta), 2);
result_lambda_hat = zeros(length(gamma), length(m), length(r), legnth(phi_h), length(theta), length(alpha), length(epsilon), length(beta), 2);
result_flag = zeros(length(gamma), length(m), length(r), legnth(phi_h), length(theta), length(alpha), length(epsilon), length(beta), 1);

bad_files = {};

for i = 1:length(gamma)
    for j = 1:length(m)
        for k = 1:length(phi_h)
            for l = 1:length(r)
                for f = 1:length(theta)
                    for g = 1:length(alpha)
                        for e = 1:length(epsilon)
                            for o = 1:length(beta)
                                % for o = 1:length(beta)
                        
                                p.r = r(l);
                                p.alpha = alpha(g);
                                p.theta = theta(f);
                                p.m_bar = m(j);
                                p.gamma = gamma(i);
                                p.phi_h = phi_h(k);
                                p.epsilon = epsilon(e);
                                p.beta = beta(o);
                                p.sigma =3.0;
                                file_name = sprintf('epsilon=%.2f, gamma=%.2f, sigma=%d, phi_h=%.2f, m_bar=%.2f, r=%.2f, alpha=%.2f, theta=%.2f, beta=%.2f', p.epsilon, p.gamma, p.sigma, p.phi_h, p.m_bar, p.r, p.alpha, p.theta, p.beta);
                                % data = load(['./data_beta/', file_name, '.mat'] , 'p', 'np' , 'eqm_save', 'iter_history');

                                try
                                    data = load(['./data_beta/', file_name, '.mat'] , 'p', 'np' , 'eqm_save', 'iter_history');
                                catch ME
                                    fprintf('[LOAD ERROR] Failed to load: %s.mat\nError: %s\n\n', file_name, ME.message);
                                    bad_files{end+1} = file_name;
                                    continue;
                                end
                                % Regular expression to extract the iteration number

                                % % Directory containing figure files
                                % folder_path = './figure_beta';

                                % % Get list of all PNG files in the folder
                                % file_list = dir(fullfile(folder_path, '*.png'));

                                % % Loop through each file and extract iteration number
                                % for kk = 1:length(file_list)
                                %     file_name_var = file_list(kk).name;
                                    
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


                                p_iter = data.p;
                                np_iter = data.np;
                                eqm_iter = data.eqm_save;
                                iter_history = data.iter_history;
                                
                                lambda_interp_l = @(q) interp1(np_iter.q, eqm_iter.lambda(1,:), q, 'spline', 'extrap');
                                lambda_interp_h = @(q) interp1(np_iter.q, eqm_iter.lambda(2,:), q, 'spline', 'extrap');
                                omega_interp_l = @(q) interp1(np_iter.q, eqm_iter.omega(1,:), q, 'spline', 'extrap');
                                omega_interp_h = @(q) interp1(np_iter.q, eqm_iter.omega(2,:), q, 'spline', 'extrap');
                                omega_entry = @(q) 1./(q * p_iter.tau * sqrt(2*pi)) .* exp(-(log(q) - p_iter.iota).^2 / (2 * p_iter.tau^2));

                                result_Expected_q(i,j, k, l,f,g, e, o, 1) = integral(@(q) q .* omega_interp_l(q), np_iter.q_min, np_iter.q_max);  % Choose appropriate limits
                                result_10quantile_q(i,j, k, l,f,g, e, o, 1) = fzero(@(q) integral(omega_interp_l, np_iter.q_min, q) - 0.10, [np_iter.q_min, np_iter.q_max]);
                                result_50quantile_q(i,j, k, l,f,g, e, o, 1) = fzero(@(q) integral(omega_interp_l, np_iter.q_min, q) - 0.50, [np_iter.q_min, np_iter.q_max]);
                                result_90quantile_q(i,j, k, l,f,g, e, o, 1) = fzero(@(q) integral(omega_interp_l, np_iter.q_min, q) - 0.90, [np_iter.q_min, np_iter.q_max]);



                                result_Expected_q(i,j, k, l,f,g, e, o, 2) = integral(@(q) q .* omega_interp_h(q), np_iter.q_min, np_iter.q_max);  % Choose appropriate limits
                                result_10quantile_q(i,j, k, l,f,g, e, o, 2) = fzero(@(q) integral(omega_interp_h, np_iter.q_min, q) - 0.10, [np_iter.q_min, np_iter.q_max]);
                                result_50quantile_q(i,j, k, l,f,g, e, o, 2) = fzero(@(q) integral(omega_interp_h, np_iter.q_min, q) - 0.50, [np_iter.q_min, np_iter.q_max]);
                                result_90quantile_q(i,j, k, l,f,g, e, o, 2) = fzero(@(q) integral(omega_interp_h, np_iter.q_min, q) - 0.90, [np_iter.q_min, np_iter.q_max]);



                                result_Expected_q(i,j, k, l,f,g, e, o, 3) = integral(@(q) q .* omega_entry(q), np_iter.q_min, np_iter.q_max);  % Choose appropriate limits
                                result_10quantile_q(i,j, k, l,f,g, e, o, 3) = fzero(@(q) integral(omega_entry, np_iter.q_min, q) - 0.10, [np_iter.q_min, np_iter.q_max]);
                                result_50quantile_q(i,j, k, l,f,g, e, o, 3) = fzero(@(q) integral(omega_entry, np_iter.q_min, q) - 0.50, [np_iter.q_min, np_iter.q_max]);
                                result_90quantile_q(i,j, k, l,f,g, e, o, 3) = fzero(@(q) integral(omega_entry, np_iter.q_min, q) - 0.90, [np_iter.q_min, np_iter.q_max]);


                                % result_s = zeros(length(gamma), length(m), length(r), length(theta), length(alpha),2);
                                % result_g_3 = zeros(length(gamma), length(m), length(r), length(theta), length(alpha),3);
                                % result_mu_hat = zeros(length(gamma), length(m), length(r), length(theta), length(alpha),2,2);
                                % result_nu_hat = zeros(length(gamma), length(m), length(r), length(theta), length(alpha),2);
                                % result_delta_hat = zeros(length(gamma), length(m), length(r), length(theta), length(alpha),2);
                                % result_flag = zeros(length(gamma), length(m), length(r), length(theta), length(alpha),1);

                                result_s(i,j, k, l,f,g, e, o, :) = eqm_iter.s;
                                % result_g_3(i,j, k, l,f,g, e, o, 1) = (eqm_iter.gp1_type(1))^(1/p_iter.sigma-1)-1;
                                % result_g_3(i,j, k, l,f,g, e, o, 2) = (eqm_iter.gp1_type(2))^(1/p_iter.sigma-1)-1;
                                result_g_3(i,j, k, l,f,g, e, o, 1) = (eqm_iter.gp1_type(1))-1;
                                result_g_3(i,j, k, l,f,g, e, o, 2) = (eqm_iter.gp1_type(2))-1;
                                result_g_3(i,j, k, l,f,g, e, o, 3) = eqm_iter.g;

                                result_mu_hat(i,j, k, l,f,g, e, o, :, :) = eqm_iter.mu_hat;
                                result_nu_hat(i,j, k, l,f,g, e, o, :) = eqm_iter.nu_hat;
                                result_delta_hat(i,j, k, l,f,g, e, o, :) = eqm_iter.delta_hat;
                                result_lambda_hat(i,j, k, l,f,g, e, o, 1) = integral(@(q) omega_interp_l(q) .* lambda_interp_l(q), np_iter.q_min, np_iter.q_max);
                                result_lambda_hat(i,j, k, l,f,g, e, o, 2) = integral(@(q) omega_interp_h(q) .* lambda_interp_h(q), np_iter.q_min, np_iter.q_max);


                                if p_iter.iter<800
                                    result_flag(i,j, k, l,f,g, e, o, 1) = 2;
                                elseif p_iter.iter==800 && p_iter.error<1e-4
                                    result_flag(i,j, k, l,f,g, e, o, 1) = 1;
                                elseif p_iter.iter==800 && p_iter.error>=1e-4
                                    result_flag(i,j, k, l,f,g, e, o, 1) = 0;
                                end


                                % if iter<800
                                %     result_flag(i,j, l,f,g, 1) = 1;
                                % else
                                %     result_flag(i,j, l,f,g, 1) = 0;
                                % end



                            end
                        end
                    end
                end
            end
        end
    end
end


figure_save_path = './figure_analysis_more';

if ~exist(figure_save_path, 'dir')
    mkdir(figure_save_path);
end

pointsize= 100;

%%
color_set = {'r', 'g', 'b', 'm', 'c', 'k'};  % one is orange


% gamma = [0.10 0.20];

% epsilon = [0.0, 1.0];

for i = 1:length(gamma)
    for e = 1:length(epsilon)
        for l = 1:length(r)

            figure;

            subplot(1, 2, 1);
            hold on;
            legend_handles = []; % Initialize an array to store scatter handles for the legend
            legend_labels = {};  % Initialize legend labels

            legend_num =0;
            for g = 1:length(alpha)
                for o = 1:length(beta)
                    legend_num = legend_num + 1;
                    color = color_set{legend_num};
                    legend_labels{legend_num} = sprintf('$\\alpha = %.2f, \\beta=%.2f$', alpha(g), beta(o)); % Add legend label for each alpha

                    % Create a dummy scatter plot for the legend
                    h = scatter(NaN, NaN, pointsize, color, 'filled');
                    legend_handles = [legend_handles, h];

                    for f = 1:length(theta)
                        for k = 1:length(phi_h)
                            for j = 1:length(m)

                                if result_flag(i,j, k, l,f,g, e, o, 1) == 2
                                    scatter(m(j), result_s(i,j, k, l,f,g, e, o, 1), pointsize, color, 'filled');
                                elseif result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                    scatter(m(j), result_s(i,j, k, l,f,g, e, o, 1), pointsize, 'o', 'MarkerEdgeColor', color, 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
                                elseif result_flag(i,j, k, l,f,g, e, o, 1) == 0
                                    scatter(m(j), result_s(i,j, k, l,f,g, e, o, 1), pointsize, 'x', 'MarkerEdgeColor', color, 'LineWidth', 1.5);
                                end
                                % if result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 1), pointsize, color, 'filled');
                                % % elseif result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                % %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 1), pointsize, 'o', 'MarkerEdgeColor', color);
                                % elseif result_flag(i,j, k, l,f,g, e, o, 1) == 0
                                %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 1), pointsize, 'x', 'MarkerEdgeColor', color, 'LineWidth', 1.5);
                                % end
                                if selector(gamma(i), m(j), phi_h(k), r(l), theta(f), alpha(g), epsilon(e), beta(o)) == 1
                                    h_txt = text(m(j), result_s(i,j, k, l,f,g, e, o, 1) + 0.2, ...
                                                                                    'Scenario 1', ...
                                                                                    'Color', 'k', ...
                                                                                    'FontSize', 20, ...
                                                                                    'Interpreter', 'latex', ...
                                                                                    'BackgroundColor', 'w', ...      % Boxed background
                                                                                    'EdgeColor', 'k', ...            % Box border
                                                                                    'Margin', 4);                    % Padding inside box
                                    plot(m(j), result_s(i,j, k, l,f,g, e, o, 1), ...
                                    'p', ...                % Star marker
                                    'MarkerSize', 30, ...   % Large size
                                    'MarkerEdgeColor', 'k', ...
                                    'MarkerFaceColor', 'none', ... % Hollow
                                    'LineWidth', 3);
                                end

                                if selector(gamma(i), m(j), phi_h(k), r(l), theta(f), alpha(g), epsilon(e), beta(o)) == 2
                                    h_txt = text(m(j), result_s(i,j, k, l,f,g, e, o, 1) + 0.2, ...
                                                                                    'Scenario 2', ...
                                                                                    'Color', 'k', ...
                                                                                    'FontSize', 20, ...
                                                                                    'Interpreter', 'latex', ...
                                                                                    'BackgroundColor', 'w', ...      % Boxed background
                                                                                    'EdgeColor', 'k', ...            % Box border
                                                                                    'Margin', 4);                    % Padding inside box
                                    plot(m(j), result_s(i,j, k, l,f,g, e, o, 1), ...
                                    'p', ...                % Star marker
                                    'MarkerSize', 30, ...   % Large size
                                    'MarkerEdgeColor', 'k', ...
                                    'MarkerFaceColor', 'none', ... % Hollow
                                    'LineWidth', 3);
                                end
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
            legend_num =0;

            for g = 1:length(alpha)
                for o = 1:length(beta)
                    legend_num = legend_num + 1;

                    color = color_set{legend_num};
                    legend_labels{legend_num} = sprintf('$\\alpha = %.2f, \\beta=%.2f$', alpha(g), beta(o)); % Add legend label for each alpha

                    % Create a dummy scatter plot for the legend
                    h = scatter(NaN, NaN, pointsize, color, 'filled');
                    legend_handles = [legend_handles, h];

                    for f = 1:length(theta)
                        for k = 1:length(phi_h)
                            for j = 1:length(m)

                                
                                if result_flag(i,j, k, l,f,g, e, o, 1) == 2
                                    scatter(m(j), result_s(i,j, k, l,f,g, e, o, 2), pointsize, color, 'filled');
                                elseif result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                    scatter(m(j), result_s(i,j, k, l,f,g, e, o, 2), pointsize, 'o', 'MarkerEdgeColor', color, 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
                                elseif result_flag(i,j, k, l,f,g, e, o, 1) == 0
                                    scatter(m(j), result_s(i,j, k, l,f,g, e, o, 2), pointsize, 'x', 'MarkerEdgeColor', color, 'LineWidth', 1.5);
                                end

                                % if result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 2), pointsize, color, 'filled');
                                % % elseif result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                % %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 2), pointsize, 'o', 'MarkerEdgeColor', color, 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
                                % elseif result_flag(i,j, k, l,f,g, e, o, 1) == 0
                                %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 2), pointsize, 'x', 'MarkerEdgeColor', color, 'LineWidth', 1.5);
                                % end

                                if selector(gamma(i), m(j), phi_h(k), r(l), theta(f), alpha(g), epsilon(e), beta(o)) == 1
                                    h_txt = text(m(j), result_s(i,j, k, l,f,g, e, o, 2) + 0.2, ...
                                                                                    'Scenario 1', ...
                                                                                    'Color', 'k', ...
                                                                                    'FontSize', 20, ...
                                                                                    'Interpreter', 'latex', ...
                                                                                    'BackgroundColor', 'w', ...      % Boxed background
                                                                                    'EdgeColor', 'k', ...            % Box border
                                                                                    'Margin', 4);                    % Padding inside box
                                    plot(m(j), result_s(i,j, k, l,f,g, e, o, 2), ...
                                    'p', ...                % Star marker
                                    'MarkerSize', 30, ...   % Large size
                                    'MarkerEdgeColor', 'k', ...
                                    'MarkerFaceColor', 'none', ... % Hollow
                                    'LineWidth', 3);
                                end

                                if selector(gamma(i), m(j), phi_h(k), r(l), theta(f), alpha(g), epsilon(e), beta(o)) == 2
                                    h_txt = text(m(j), result_s(i,j, k, l,f,g, e, o, 2) + 0.2, ...
                                                                                    'Scenario 2', ...
                                                                                    'Color', 'k', ...
                                                                                    'FontSize', 20, ...
                                                                                    'Interpreter', 'latex', ...
                                                                                    'BackgroundColor', 'w', ...      % Boxed background
                                                                                    'EdgeColor', 'k', ...            % Box border
                                                                                    'Margin', 4);                    % Padding inside box
                                    plot(m(j), result_s(i,j, k, l,f,g, e, o, 2), ...
                                    'p', ...                % Star marker
                                    'MarkerSize', 30, ...   % Large size
                                    'MarkerEdgeColor', 'k', ...
                                    'MarkerFaceColor', 'none', ... % Hollow
                                    'LineWidth', 3);
                                end

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

            sgtitle(sprintf('$\\gamma = %.2f, \\epsilon=%.2f, r=%.2f$', gamma(i), epsilon(e), r(l)), 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
            set(gcf, 'Position', [100, 100, 2200, 1000]); % Set figure size
            saveas(gcf, sprintf ([figure_save_path, '/MarketShare_gamma=%.2f_epsilon=%.2f_r=%.2f.png'], gamma(i), epsilon(e), r(l)));
            
        end
    end
end





%%
color_set = {'r', 'g', 'b', 'm', 'c', 'k'};  % one is orange



for i = 1:length(gamma)
    for e = 1:length(epsilon)
        for l = 1:length(r)

            figure;

            % subplot(1, 2, 1);
            hold on;
            legend_handles = []; % Initialize an array to store scatter handles for the legend
            legend_labels = {};  % Initialize legend labels

            legend_num =0;
            for g = 1:length(alpha)
                for o = 1:length(beta)
                    legend_num = legend_num + 1;
                    color = color_set{legend_num};
                    legend_labels{legend_num} = sprintf('$\\alpha = %.2f, \\beta=%.2f$', alpha(g), beta(o)); % Add legend label for each alpha

                    % Create a dummy scatter plot for the legend
                    h = scatter(NaN, NaN, pointsize, color, 'filled');
                    legend_handles = [legend_handles, h];

                    for f = 1:length(theta)
                        for k = 1:length(phi_h)
                            for j = 1:length(m)

                                if result_flag(i,j, k, l,f,g, e, o, 1) == 2
                                    scatter(m(j), result_g_3(i,j, k, l,f,g, e, o, 3), pointsize, color, 'filled');
                                elseif result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                    scatter(m(j), result_g_3(i,j, k, l,f,g, e, o, 3), pointsize, 'o', 'MarkerEdgeColor', color, 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
                                elseif result_flag(i,j, k, l,f,g, e, o, 1) == 0
                                    scatter(m(j), result_g_3(i,j, k, l,f,g, e, o, 3), pointsize, 'x', 'MarkerEdgeColor', color, 'LineWidth', 1.5);
                                end
                                % if result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 1), pointsize, color, 'filled');
                                % % elseif result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                % %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 1), pointsize, 'o', 'MarkerEdgeColor', color);
                                % elseif result_flag(i,j, k, l,f,g, e, o, 1) == 0
                                %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 1), pointsize, 'x', 'MarkerEdgeColor', color, 'LineWidth', 1.5);
                                % end

                                if selector(gamma(i), m(j), phi_h(k), r(l), theta(f), alpha(g), epsilon(e), beta(o)) == 1
                                    h_txt = text(m(j), result_g_3(i,j, k, l,f,g, e, o, 3) - 0.1, ...
                                                                                    'Scenario 1', ...
                                                                                    'Color', 'k', ...
                                                                                    'FontSize', 20, ...
                                                                                    'Interpreter', 'latex', ...
                                                                                    'BackgroundColor', 'w', ...      % Boxed background
                                                                                    'EdgeColor', 'k', ...            % Box border
                                                                                    'Margin', 4);                    % Padding inside box
                                    plot(m(j), result_g_3(i,j, k, l,f,g, e, o, 3), ...
                                    'p', ...                % Star marker
                                    'MarkerSize', 30, ...   % Large size
                                    'MarkerEdgeColor', 'k', ...
                                    'MarkerFaceColor', 'none', ... % Hollow
                                    'LineWidth', 3);
                                end

                                if selector(gamma(i), m(j), phi_h(k), r(l), theta(f), alpha(g), epsilon(e), beta(o)) == 2
                                    h_txt = text(m(j), result_g_3(i,j, k, l,f,g, e, o, 3) - 0.1, ...
                                                                                    'Scenario 2', ...
                                                                                    'Color', 'k', ...
                                                                                    'FontSize', 20, ...
                                                                                    'Interpreter', 'latex', ...
                                                                                    'BackgroundColor', 'w', ...      % Boxed background
                                                                                    'EdgeColor', 'k', ...            % Box border
                                                                                    'Margin', 4);                    % Padding inside box
                                    plot(m(j), result_g_3(i,j, k, l,f,g, e, o, 3), ...
                                    'p', ...                % Star marker
                                    'MarkerSize', 30, ...   % Large size
                                    'MarkerEdgeColor', 'k', ...
                                    'MarkerFaceColor', 'none', ... % Hollow
                                    'LineWidth', 3);
                                end


                            end
                        end
                    end

                end
            end

            xlabel('$\bar{m}$', 'Interpreter', 'latex', 'FontSize', 20);
            ylabel('$g$', 'Interpreter', 'latex', 'FontSize', 20);
            % ylim([0, 1]);
            yline(0, 'k--', 'LineWidth', 1.5);
            grid on;
            legend(legend_handles, legend_labels, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best'); % Set legend font size
            title(sprintf('Aggregate Growth Rate $\\gamma = %.2f, \\epsilon=%.2f, r=%.2f$', gamma(i), epsilon(e), r(l)), 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
            hold off;

            % sgtitle('Interpreter', 'latex', 'FontSize', 20); % Set title font size

            set(gcf, 'Position', [100, 100, 1600, 1200]); % Set figure size
            saveas(gcf, sprintf ([figure_save_path, '/AggregateGrowthRate_gamma=%.2f_epsilon=%.2f_r=%.2f.png'], gamma(i), epsilon(e), r(l)));
            
        end
    end
end



color_set = {'r', 'g', 'b', 'm', 'c', 'k'};

for i = 1:length(gamma)
    for e = 1:length(epsilon)
        for l = 1:length(r)

            figure;

            legend_handles = [];
            legend_labels = {};

            for g = 1:length(alpha)
                subplot(1, 2, g);
                hold on;
                legend_handles = [];
                legend_labels = {};

                for o = 1:length(beta)
                    % Assign color based on g and o
                    if g == 1
                        color = color_set{o};      % o=1: 'r', o=2: 'g'
                    elseif g == 2
                        color = color_set{o+2};    % o=1: 'b', o=2: 'm'
                    end

                    legend_labels{o} = sprintf('$\\alpha = %.2f, \\beta=%.2f$', alpha(g), beta(o));
                    h = scatter(NaN, NaN, pointsize, color, 'filled');
                    legend_handles = [legend_handles, h];

                    for f = 1:length(theta)
                        for k = 1:length(phi_h)
                            for j = 1:length(m)

                                if result_flag(i,j, k, l,f,g, e, o, 1) == 2
                                    scatter(m(j), result_g_3(i,j, k, l,f,g, e, o, 3), pointsize, color, 'filled');
                                elseif result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                    scatter(m(j), result_g_3(i,j, k, l,f,g, e, o, 3), pointsize, 'o', 'MarkerEdgeColor', color, 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
                                elseif result_flag(i,j, k, l,f,g, e, o, 1) == 0
                                    scatter(m(j), result_g_3(i,j, k, l,f,g, e, o, 3), pointsize, 'x', 'MarkerEdgeColor', color, 'LineWidth', 1.5);
                                end

                                if selector(gamma(i), m(j), phi_h(k), r(l), theta(f), alpha(g), epsilon(e), beta(o)) == 1
                                    h_txt = text(m(j), result_g_3(i,j, k, l,f,g, e, o, 3) - 0.1, ...
                                                                                    'Scenario 1', ...
                                                                                    'Color', 'k', ...
                                                                                    'FontSize', 20, ...
                                                                                    'Interpreter', 'latex', ...
                                                                                    'BackgroundColor', 'w', ...      % Boxed background
                                                                                    'EdgeColor', 'k', ...            % Box border
                                                                                    'Margin', 4);                    % Padding inside box
                                    plot(m(j), result_g_3(i,j, k, l,f,g, e, o, 3), ...
                                    'p', ...                % Star marker
                                    'MarkerSize', 30, ...   % Large size
                                    'MarkerEdgeColor', 'k', ...
                                    'MarkerFaceColor', 'none', ... % Hollow
                                    'LineWidth', 3);
                                end

                                if selector(gamma(i), m(j), phi_h(k), r(l), theta(f), alpha(g), epsilon(e), beta(o)) == 2
                                    h_txt = text(m(j), result_g_3(i,j, k, l,f,g, e, o, 3) - 0.1, ...
                                                                                    'Scenario 2', ...
                                                                                    'Color', 'k', ...
                                                                                    'FontSize', 20, ...
                                                                                    'Interpreter', 'latex', ...
                                                                                    'BackgroundColor', 'w', ...      % Boxed background
                                                                                    'EdgeColor', 'k', ...            % Box border
                                                                                    'Margin', 4);                    % Padding inside box
                                    plot(m(j), result_g_3(i,j, k, l,f,g, e, o, 3), ...
                                    'p', ...                % Star marker
                                    'MarkerSize', 30, ...   % Large size
                                    'MarkerEdgeColor', 'k', ...
                                    'MarkerFaceColor', 'none', ... % Hollow
                                    'LineWidth', 3);
                                end


                            end
                        end
                    end
                end

                xlabel('$\bar{m}$', 'Interpreter', 'latex', 'FontSize', 20);
                ylabel('$g$', 'Interpreter', 'latex', 'FontSize', 20);
                if gamma(i) == 0.1
                    ylim([-0.3, 0.11]);
                elseif gamma(i) == 0.2
                    ylim([-0.3, 0.11]);
                end
                yline(0, 'k--', 'LineWidth', 1.5);
                grid on;
                legend(legend_handles, legend_labels, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best');
            title(sprintf('Aggregate Growth Rate $\\gamma = %.2f, \\epsilon=%.2f, r=%.2f$', gamma(i), epsilon(e), r(l)), 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
                hold off;
            end
            sgtitle(sprintf('$\\gamma = %.2f, \\epsilon=%.2f, r=%.2f$', gamma(i), epsilon(e), r(l)), 'Interpreter', 'latex', 'FontSize', 20); % Set title font size

            set(gcf, 'Position', [100, 100, 2200, 1000]);
            saveas(gcf, sprintf([figure_save_path, '/AggregateGrowthRate_sbs_gamma=%.2f_epsilon=%.2f_r=%.2f.png'], gamma(i), epsilon(e), r(l)));
        end
    end
end









%%
color_set = {'r', 'g', 'b', 'm', 'c', 'k'};  % one is orange


gamma = [0.10 0.20];

epsilon = [0.0, 1.0];

for i = 1:length(gamma)
    for e = 1:length(epsilon)
        for l = 1:length(r)

            figure;

            subplot(1, 2, 1);
            hold on;
            legend_handles = []; % Initialize an array to store scatter handles for the legend
            legend_labels = {};  % Initialize legend labels

            legend_num =0;
            for g = 1:length(alpha)
                for o = 1:length(beta)
                    legend_num = legend_num + 1;
                    color = color_set{legend_num};
                    legend_labels{legend_num} = sprintf('$\\alpha = %.2f, \\beta=%.2f$', alpha(g), beta(o)); % Add legend label for each alpha

                    % Create a dummy scatter plot for the legend
                    h = scatter(NaN, NaN, pointsize, color, 'filled');
                    legend_handles = [legend_handles, h];

                    for f = 1:length(theta)
                        for k = 1:length(phi_h)
                            for j = 1:length(m)

                                if result_flag(i,j, k, l,f,g, e, o, 1) == 2
                                    scatter(m(j), result_g_3(i,j, k, l,f,g, e, o, 1), pointsize, color, 'filled');
                                elseif result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                    scatter(m(j), result_g_3(i,j, k, l,f,g, e, o, 1), pointsize, 'o', 'MarkerEdgeColor', color, 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
                                elseif result_flag(i,j, k, l,f,g, e, o, 1) == 0
                                    scatter(m(j), result_g_3(i,j, k, l,f,g, e, o, 1), pointsize, 'x', 'MarkerEdgeColor', color, 'LineWidth', 1.5);
                                end
                                % if result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 1), pointsize, color, 'filled');
                                % % elseif result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                % %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 1), pointsize, 'o', 'MarkerEdgeColor', color);
                                % elseif result_flag(i,j, k, l,f,g, e, o, 1) == 0
                                %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 1), pointsize, 'x', 'MarkerEdgeColor', color, 'LineWidth', 1.5);
                                % end

                                if selector(gamma(i), m(j), phi_h(k), r(l), theta(f), alpha(g), epsilon(e), beta(o)) == 1
                                    h_txt = text(m(j), result_g_3(i,j, k, l,f,g, e, o, 1) + 0.05, ...
                                                                                    'Scenario 1', ...
                                                                                    'Color', 'k', ...
                                                                                    'FontSize', 20, ...
                                                                                    'Interpreter', 'latex', ...
                                                                                    'BackgroundColor', 'w', ...      % Boxed background
                                                                                    'EdgeColor', 'k', ...            % Box border
                                                                                    'Margin', 4);                    % Padding inside box
                                    plot(m(j), result_g_3(i,j, k, l,f,g, e, o, 1), ...
                                    'p', ...                % Star marker
                                    'MarkerSize', 30, ...   % Large size
                                    'MarkerEdgeColor', 'k', ...
                                    'MarkerFaceColor', 'none', ... % Hollow
                                    'LineWidth', 3);
                                end

                                if selector(gamma(i), m(j), phi_h(k), r(l), theta(f), alpha(g), epsilon(e), beta(o)) == 2
                                    h_txt = text(m(j), result_g_3(i,j, k, l,f,g, e, o, 1) + 0.05, ...
                                                                                    'Scenario 2', ...
                                                                                    'Color', 'k', ...
                                                                                    'FontSize', 20, ...
                                                                                    'Interpreter', 'latex', ...
                                                                                    'BackgroundColor', 'w', ...      % Boxed background
                                                                                    'EdgeColor', 'k', ...            % Box border
                                                                                    'Margin', 4);                    % Padding inside box
                                    plot(m(j), result_g_3(i,j, k, l,f,g, e, o, 1), ...
                                    'p', ...                % Star marker
                                    'MarkerSize', 30, ...   % Large size
                                    'MarkerEdgeColor', 'k', ...
                                    'MarkerFaceColor', 'none', ... % Hollow
                                    'LineWidth', 3);
                                end
                            end
                        end
                    end

                end
            end

            xlabel('$\bar{m}$', 'Interpreter', 'latex', 'FontSize', 20);
            ylabel('$s_l$', 'Interpreter', 'latex', 'FontSize', 20);
            % ylim([0, 1]);
            yline(0.0, 'k--', 'LineWidth', 1.5);
            grid on;
            legend(legend_handles, legend_labels, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best'); % Set legend font size
            title('Growth Rate For Type L', 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
            hold off;


            subplot(1, 2, 2);
            hold on;
            legend_handles = []; % Initialize an array to store scatter handles for the legend
            legend_labels = {};  % Initialize legend labels
            legend_num =0;

            for g = 1:length(alpha)
                for o = 1:length(beta)
                    legend_num = legend_num + 1;

                    color = color_set{legend_num};
                    legend_labels{legend_num} = sprintf('$\\alpha = %.2f, \\beta=%.2f$', alpha(g), beta(o)); % Add legend label for each alpha

                    % Create a dummy scatter plot for the legend
                    h = scatter(NaN, NaN, pointsize, color, 'filled');
                    legend_handles = [legend_handles, h];

                    for f = 1:length(theta)
                        for k = 1:length(phi_h)
                            for j = 1:length(m)

                                
                                if result_flag(i,j, k, l,f,g, e, o, 1) == 2
                                    scatter(m(j), result_g_3(i,j, k, l,f,g, e, o, 2), pointsize, color, 'filled');
                                elseif result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                    scatter(m(j), result_g_3(i,j, k, l,f,g, e, o, 2), pointsize, 'o', 'MarkerEdgeColor', color, 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
                                elseif result_flag(i,j, k, l,f,g, e, o, 1) == 0
                                    scatter(m(j), result_g_3(i,j, k, l,f,g, e, o, 2), pointsize, 'x', 'MarkerEdgeColor', color, 'LineWidth', 1.5);
                                end

                                % if result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 2), pointsize, color, 'filled');
                                % % elseif result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                % %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 2), pointsize, 'o', 'MarkerEdgeColor', color, 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
                                % elseif result_flag(i,j, k, l,f,g, e, o, 1) == 0
                                %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 2), pointsize, 'x', 'MarkerEdgeColor', color, 'LineWidth', 1.5);
                                % end


                                if selector(gamma(i), m(j), phi_h(k), r(l), theta(f), alpha(g), epsilon(e), beta(o)) == 1
                                    h_txt = text(m(j), result_g_3(i,j, k, l,f,g, e, o, 2) + 0.05, ...
                                                                                    'Scenario 1', ...
                                                                                    'Color', 'k', ...
                                                                                    'FontSize', 20, ...
                                                                                    'Interpreter', 'latex', ...
                                                                                    'BackgroundColor', 'w', ...      % Boxed background
                                                                                    'EdgeColor', 'k', ...            % Box border
                                                                                    'Margin', 4);                    % Padding inside box
                                    plot(m(j), result_g_3(i,j, k, l,f,g, e, o, 2), ...
                                    'p', ...                % Star marker
                                    'MarkerSize', 30, ...   % Large size
                                    'MarkerEdgeColor', 'k', ...
                                    'MarkerFaceColor', 'none', ... % Hollow
                                    'LineWidth', 3);
                                end

                                if selector(gamma(i), m(j), phi_h(k), r(l), theta(f), alpha(g), epsilon(e), beta(o)) == 2
                                    h_txt = text(m(j), result_g_3(i,j, k, l,f,g, e, o, 2) + 0.05, ...
                                                                                    'Scenario 2', ...
                                                                                    'Color', 'k', ...
                                                                                    'FontSize', 20, ...
                                                                                    'Interpreter', 'latex', ...
                                                                                    'BackgroundColor', 'w', ...      % Boxed background
                                                                                    'EdgeColor', 'k', ...            % Box border
                                                                                    'Margin', 4);                    % Padding inside box
                                    plot(m(j), result_g_3(i,j, k, l,f,g, e, o, 2), ...
                                    'p', ...                % Star marker
                                    'MarkerSize', 30, ...   % Large size
                                    'MarkerEdgeColor', 'k', ...
                                    'MarkerFaceColor', 'none', ... % Hollow
                                    'LineWidth', 3);
                                end
                            end
                        end
                    end


                end
            end

            xlabel('$\bar{m}$', 'Interpreter', 'latex', 'FontSize', 20);
            ylabel('$g_h$', 'Interpreter', 'latex', 'FontSize', 20);
            % ylim([0, 1]);
            grid on;
            yline(0.0, 'k--', 'LineWidth', 1.5);
            legend(legend_handles, legend_labels, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best'); % Set legend font size
            title('Growth Rate For Type H', 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
            hold off;

            sgtitle(sprintf('$\\gamma = %.2f, \\epsilon=%.2f, r=%.2f$', gamma(i), epsilon(e), r(l)), 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
            set(gcf, 'Position', [100, 100, 2200, 1000]); % Set figure size
            saveas(gcf, sprintf ([figure_save_path, '/TypeGrowthRate_gamma=%.2f_epsilon=%.2f_r=%.2f.png'], gamma(i), epsilon(e), r(l)));
            
        end
    end
end












%%
color_set = {'r', 'g', 'b', 'm', 'c', 'k'};  % one is orange


for i = 1:length(gamma)
    for e = 1:length(epsilon)
        for l = 1:length(r)

            figure;

            subplot(1, 2, 1);
            hold on;
            legend_handles = []; % Initialize an array to store scatter handles for the legend
            legend_labels = {};  % Initialize legend labels

            legend_num =0;
            for g = 1:length(alpha)
                for o = 1:length(beta)
                    legend_num = legend_num + 1;
                    color = color_set{legend_num};
                    legend_labels{legend_num} = sprintf('$\\alpha = %.2f, \\beta=%.2f$', alpha(g), beta(o)); % Add legend label for each alpha

                    % Create a dummy scatter plot for the legend
                    h = scatter(NaN, NaN, pointsize, color, 'filled');
                    legend_handles = [legend_handles, h];

                    for f = 1:length(theta)
                        for k = 1:length(phi_h)
                            for j = 1:length(m)

                                if result_flag(i,j, k, l,f,g, e, o, 1) == 2
                                    scatter(m(j), result_mu_hat(i,j, k, l,f,g, e, o, 1,1) + result_mu_hat(i,j, k, l,f,g, e, o, 1,2), pointsize, color, 'filled');
                                elseif result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                    scatter(m(j), result_mu_hat(i,j, k, l,f,g, e, o, 1,1) + result_mu_hat(i,j, k, l,f,g, e, o, 1,2), pointsize, 'o', 'MarkerEdgeColor', color, 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
                                elseif result_flag(i,j, k, l,f,g, e, o, 1) == 0
                                    scatter(m(j), result_mu_hat(i,j, k, l,f,g, e, o, 1,1) + result_mu_hat(i,j, k, l,f,g, e, o, 1,2), pointsize, 'x', 'MarkerEdgeColor', color, 'LineWidth', 1.5);
                                end
                                % if result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 1), pointsize, color, 'filled');
                                % % elseif result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                % %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 1), pointsize, 'o', 'MarkerEdgeColor', color);
                                % elseif result_flag(i,j, k, l,f,g, e, o, 1) == 0
                                %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 1), pointsize, 'x', 'MarkerEdgeColor', color, 'LineWidth', 1.5);
                                % end

                                if selector(gamma(i), m(j), phi_h(k), r(l), theta(f), alpha(g), epsilon(e), beta(o)) == 1
                                    h_txt = text(m(j), result_mu_hat(i,j, k, l,f,g, e, o, 1,1) + result_mu_hat(i,j, k, l,f,g, e, o, 1,2) + 0.4, ...
                                                                                    'Scenario 1', ...
                                                                                    'Color', 'k', ...
                                                                                    'FontSize', 20, ...
                                                                                    'Interpreter', 'latex', ...
                                                                                    'BackgroundColor', 'w', ...      % Boxed background
                                                                                    'EdgeColor', 'k', ...            % Box border
                                                                                    'Margin', 4);                    % Padding inside box
                                    plot(m(j), result_mu_hat(i,j, k, l,f,g, e, o, 1,1) + result_mu_hat(i,j, k, l,f,g, e, o, 1,2), ...
                                    'p', ...                % Star marker
                                    'MarkerSize', 30, ...   % Large size
                                    'MarkerEdgeColor', 'k', ...
                                    'MarkerFaceColor', 'none', ... % Hollow
                                    'LineWidth', 3);
                                end

                                if selector(gamma(i), m(j), phi_h(k), r(l), theta(f), alpha(g), epsilon(e), beta(o)) == 2
                                    h_txt = text(m(j), result_mu_hat(i,j, k, l,f,g, e, o, 1,1) + result_mu_hat(i,j, k, l,f,g, e, o, 1,2) + 0.4, ...
                                                                                    'Scenario 2', ...
                                                                                    'Color', 'k', ...
                                                                                    'FontSize', 20, ...
                                                                                    'Interpreter', 'latex', ...
                                                                                    'BackgroundColor', 'w', ...      % Boxed background
                                                                                    'EdgeColor', 'k', ...            % Box border
                                                                                    'Margin', 4);                    % Padding inside box
                                    plot(m(j), result_mu_hat(i,j, k, l,f,g, e, o, 1,1) + result_mu_hat(i,j, k, l,f,g, e, o, 1,2), ...
                                    'p', ...                % Star marker
                                    'MarkerSize', 30, ...   % Large size
                                    'MarkerEdgeColor', 'k', ...
                                    'MarkerFaceColor', 'none', ... % Hollow
                                    'LineWidth', 3);
                                end
                            end
                        end
                    end

                end
            end

            xlabel('$\bar{m}$', 'Interpreter', 'latex', 'FontSize', 20);
            ylabel('$\hat{\mu}_{ll} + \hat{\mu}_{lh}$', 'Interpreter', 'latex', 'FontSize', 20);
            ylim([0, 1.4]);
            % yline(0.0, 'k--', 'LineWidth', 1.5);
            grid on;
            legend(legend_handles, legend_labels, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best'); % Set legend font size
            title('Type L Total Steal Prob $\hat{\mu}_{ll} + \hat{\mu}_{lh}$', 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
            hold off;


            subplot(1, 2, 2);
            hold on;
            legend_handles = []; % Initialize an array to store scatter handles for the legend
            legend_labels = {};  % Initialize legend labels
            legend_num =0;

            for g = 1:length(alpha)
                for o = 1:length(beta)
                    legend_num = legend_num + 1;

                    color = color_set{legend_num};
                    legend_labels{legend_num} = sprintf('$\\alpha = %.2f, \\beta=%.2f$', alpha(g), beta(o)); % Add legend label for each alpha

                    % Create a dummy scatter plot for the legend
                    h = scatter(NaN, NaN, pointsize, color, 'filled');
                    legend_handles = [legend_handles, h];

                    for f = 1:length(theta)
                        for k = 1:length(phi_h)
                            for j = 1:length(m)

                                
                                if result_flag(i,j, k, l,f,g, e, o, 1) == 2
                                    scatter(m(j), result_mu_hat(i,j, k, l,f,g, e, o, 2,1) + result_mu_hat(i,j, k, l,f,g, e, o, 2,2), pointsize, color, 'filled');
                                elseif result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                    scatter(m(j), result_mu_hat(i,j, k, l,f,g, e, o, 2,1) + result_mu_hat(i,j, k, l,f,g, e, o, 2,2), pointsize, 'o', 'MarkerEdgeColor', color, 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
                                elseif result_flag(i,j, k, l,f,g, e, o, 1) == 0
                                    scatter(m(j), result_mu_hat(i,j, k, l,f,g, e, o, 2,1) + result_mu_hat(i,j, k, l,f,g, e, o, 2,2), pointsize, 'x', 'MarkerEdgeColor', color, 'LineWidth', 1.5);
                                end

                                % if result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 2), pointsize, color, 'filled');
                                % % elseif result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                % %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 2), pointsize, 'o', 'MarkerEdgeColor', color, 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
                                % elseif result_flag(i,j, k, l,f,g, e, o, 1) == 0
                                %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 2), pointsize, 'x', 'MarkerEdgeColor', color, 'LineWidth', 1.5);
                                % end

                                if selector(gamma(i), m(j), phi_h(k), r(l), theta(f), alpha(g), epsilon(e), beta(o)) == 1
                                    h_txt = text(m(j), result_mu_hat(i,j, k, l,f,g, e, o, 2,1) + result_mu_hat(i,j, k, l,f,g, e, o, 2,2) + 0.15, ...
                                                                                    'Scenario 1', ...
                                                                                    'Color', 'k', ...
                                                                                    'FontSize', 20, ...
                                                                                    'Interpreter', 'latex', ...
                                                                                    'BackgroundColor', 'w', ...      % Boxed background
                                                                                    'EdgeColor', 'k', ...            % Box border
                                                                                    'Margin', 4);                    % Padding inside box
                                    plot(m(j), result_mu_hat(i,j, k, l,f,g, e, o, 2,1) + result_mu_hat(i,j, k, l,f,g, e, o, 2,2), ...
                                    'p', ...                % Star marker
                                    'MarkerSize', 30, ...   % Large size
                                    'MarkerEdgeColor', 'k', ...
                                    'MarkerFaceColor', 'none', ... % Hollow
                                    'LineWidth', 3);
                                end

                                if selector(gamma(i), m(j), phi_h(k), r(l), theta(f), alpha(g), epsilon(e), beta(o)) == 2
                                    h_txt = text(m(j), result_mu_hat(i,j, k, l,f,g, e, o, 2,1) + result_mu_hat(i,j, k, l,f,g, e, o, 2,2) + 0.15, ...
                                                                                    'Scenario 2', ...
                                                                                    'Color', 'k', ...
                                                                                    'FontSize', 20, ...
                                                                                    'Interpreter', 'latex', ...
                                                                                    'BackgroundColor', 'w', ...      % Boxed background
                                                                                    'EdgeColor', 'k', ...            % Box border
                                                                                    'Margin', 4);                    % Padding inside box
                                    plot(m(j), result_mu_hat(i,j, k, l,f,g, e, o, 2,1) + result_mu_hat(i,j, k, l,f,g, e, o, 2,2), ...
                                    'p', ...                % Star marker
                                    'MarkerSize', 30, ...   % Large size
                                    'MarkerEdgeColor', 'k', ...
                                    'MarkerFaceColor', 'none', ... % Hollow
                                    'LineWidth', 3);
                                end

                            end
                        end
                    end


                end
            end

            xlabel('$\bar{m}$', 'Interpreter', 'latex', 'FontSize', 20);
            ylabel('$\hat{\mu}_{hl} + \hat{\mu}_{hh}$', 'Interpreter', 'latex', 'FontSize', 20);
            % ylim([0, 0.7]);
            ylim([0, 1.4]);
            grid on;
            % yline(0.0, 'k--', 'LineWidth', 1.5);
            legend(legend_handles, legend_labels, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best'); % Set legend font size
            title('Type H Total Steal Prob $\hat{\mu}_{hl} + \hat{\mu}_{hh}$', 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
            hold off;

            sgtitle(sprintf('$\\gamma = %.2f, \\epsilon=%.2f, r=%.2f$', gamma(i), epsilon(e), r(l)), 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
            set(gcf, 'Position', [100, 100, 2200, 1000]); % Set figure size
            saveas(gcf, sprintf ([figure_save_path, '/StealProbability_gamma=%.2f_epsilon=%.2f_r=%.2f.png'], gamma(i), epsilon(e), r(l)));
            
        end
    end
end









%%
color_set = {'r', 'g', 'b', 'm', 'c', 'k'};  % one is orange


for i = 1:length(gamma)
    for e = 1:length(epsilon)
        for l = 1:length(r)

            figure;

            subplot(1, 2, 1);
            hold on;
            legend_handles = []; % Initialize an array to store scatter handles for the legend
            legend_labels = {};  % Initialize legend labels

            legend_num =0;
            for g = 1:length(alpha)
                for o = 1:length(beta)
                    legend_num = legend_num + 1;
                    color = color_set{legend_num};
                    legend_labels{legend_num} = sprintf('$\\alpha = %.2f, \\beta=%.2f$', alpha(g), beta(o)); % Add legend label for each alpha

                    % Create a dummy scatter plot for the legend
                    h = scatter(NaN, NaN, pointsize, color, 'filled');
                    legend_handles = [legend_handles, h];

                    for f = 1:length(theta)
                        for k = 1:length(phi_h)
                            for j = 1:length(m)

                                if result_flag(i,j, k, l,f,g, e, o, 1) == 2
                                    scatter(m(j), result_nu_hat(i,j, k, l,f,g, e, o, 1), pointsize, color, 'filled');
                                elseif result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                    scatter(m(j), result_nu_hat(i,j, k, l,f,g, e, o, 1), pointsize, 'o', 'MarkerEdgeColor', color, 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
                                elseif result_flag(i,j, k, l,f,g, e, o, 1) == 0
                                    scatter(m(j), result_nu_hat(i,j, k, l,f,g, e, o, 1), pointsize, 'x', 'MarkerEdgeColor', color, 'LineWidth', 1.5);
                                end
                                % if result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 1), pointsize, color, 'filled');
                                % % elseif result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                % %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 1), pointsize, 'o', 'MarkerEdgeColor', color);
                                % elseif result_flag(i,j, k, l,f,g, e, o, 1) == 0
                                %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 1), pointsize, 'x', 'MarkerEdgeColor', color, 'LineWidth', 1.5);
                                % end

                                if selector(gamma(i), m(j), phi_h(k), r(l), theta(f), alpha(g), epsilon(e), beta(o)) == 1
                                    h_txt = text(m(j), result_nu_hat(i,j, k, l,f,g, e, o, 1) + 0.2, ...
                                                                                    'Scenario 1', ...
                                                                                    'Color', 'k', ...
                                                                                    'FontSize', 20, ...
                                                                                    'Interpreter', 'latex', ...
                                                                                    'BackgroundColor', 'w', ...      % Boxed background
                                                                                    'EdgeColor', 'k', ...            % Box border
                                                                                    'Margin', 4);                    % Padding inside box
                                    plot(m(j), result_nu_hat(i,j, k, l,f,g, e, o, 1), ...
                                    'p', ...                % Star marker
                                    'MarkerSize', 30, ...   % Large size
                                    'MarkerEdgeColor', 'k', ...
                                    'MarkerFaceColor', 'none', ... % Hollow
                                    'LineWidth', 3);
                                end

                                if selector(gamma(i), m(j), phi_h(k), r(l), theta(f), alpha(g), epsilon(e), beta(o)) == 2
                                    h_txt = text(m(j), result_nu_hat(i,j, k, l,f,g, e, o, 1) + 0.2, ...
                                                                                    'Scenario 2', ...
                                                                                    'Color', 'k', ...
                                                                                    'FontSize', 20, ...
                                                                                    'Interpreter', 'latex', ...
                                                                                    'BackgroundColor', 'w', ...      % Boxed background
                                                                                    'EdgeColor', 'k', ...            % Box border
                                                                                    'Margin', 4);                    % Padding inside box
                                    plot(m(j), result_nu_hat(i,j, k, l,f,g, e, o, 1), ...
                                    'p', ...                % Star marker
                                    'MarkerSize', 30, ...   % Large size
                                    'MarkerEdgeColor', 'k', ...
                                    'MarkerFaceColor', 'none', ... % Hollow
                                    'LineWidth', 3);
                                end
                            end
                        end
                    end

                end
            end

            xlabel('$\bar{m}$', 'Interpreter', 'latex', 'FontSize', 20);
            ylabel('$\hat{\nu}_{l}$', 'Interpreter', 'latex', 'FontSize', 20);
            ylim([0, 1]);
            grid on;
            % yline(0.0, 'k--', 'LineWidth', 1.5);
            legend(legend_handles, legend_labels, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best'); % Set legend font size
            title('Type L Reserve Prob $\hat{\nu}_{l}=1-\hat{\mu}_{ll}-\hat{\mu}_{hl} $', 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
            hold off;


            subplot(1, 2, 2);
            hold on;
            legend_handles = []; % Initialize an array to store scatter handles for the legend
            legend_labels = {};  % Initialize legend labels
            legend_num =0;

            for g = 1:length(alpha)
                for o = 1:length(beta)
                    legend_num = legend_num + 1;

                    color = color_set{legend_num};
                    legend_labels{legend_num} = sprintf('$\\alpha = %.2f, \\beta=%.2f$', alpha(g), beta(o)); % Add legend label for each alpha

                    % Create a dummy scatter plot for the legend
                    h = scatter(NaN, NaN, pointsize, color, 'filled');
                    legend_handles = [legend_handles, h];

                    for f = 1:length(theta)
                        for k = 1:length(phi_h)
                            for j = 1:length(m)

                                
                                if result_flag(i,j, k, l,f,g, e, o, 1) == 2
                                    scatter(m(j), result_nu_hat(i,j, k, l,f,g, e, o, 2), pointsize, color, 'filled');
                                elseif result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                    scatter(m(j), result_nu_hat(i,j, k, l,f,g, e, o, 2), pointsize, 'o', 'MarkerEdgeColor', color, 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
                                elseif result_flag(i,j, k, l,f,g, e, o, 1) == 0
                                    scatter(m(j), result_nu_hat(i,j, k, l,f,g, e, o, 2), pointsize, 'x', 'MarkerEdgeColor', color, 'LineWidth', 1.5);
                                end

                                % if result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 2), pointsize, color, 'filled');
                                % % elseif result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                % %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 2), pointsize, 'o', 'MarkerEdgeColor', color, 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
                                % elseif result_flag(i,j, k, l,f,g, e, o, 1) == 0
                                %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 2), pointsize, 'x', 'MarkerEdgeColor', color, 'LineWidth', 1.5);
                                % end

                                if selector(gamma(i), m(j), phi_h(k), r(l), theta(f), alpha(g), epsilon(e), beta(o)) == 1
                                    h_txt = text(m(j), result_nu_hat(i,j, k, l,f,g, e, o, 2) + 0.2, ...
                                                                                    'Scenario 1', ...
                                                                                    'Color', 'k', ...
                                                                                    'FontSize', 20, ...
                                                                                    'Interpreter', 'latex', ...
                                                                                    'BackgroundColor', 'w', ...      % Boxed background
                                                                                    'EdgeColor', 'k', ...            % Box border
                                                                                    'Margin', 4);                    % Padding inside box
                                    plot(m(j), result_nu_hat(i,j, k, l,f,g, e, o, 2), ...
                                    'p', ...                % Star marker
                                    'MarkerSize', 30, ...   % Large size
                                    'MarkerEdgeColor', 'k', ...
                                    'MarkerFaceColor', 'none', ... % Hollow
                                    'LineWidth', 3);
                                end

                                if selector(gamma(i), m(j), phi_h(k), r(l), theta(f), alpha(g), epsilon(e), beta(o)) == 2
                                    h_txt = text(m(j), result_nu_hat(i,j, k, l,f,g, e, o, 2) + 0.2, ...
                                                                                    'Scenario 2', ...
                                                                                    'Color', 'k', ...
                                                                                    'FontSize', 20, ...
                                                                                    'Interpreter', 'latex', ...
                                                                                    'BackgroundColor', 'w', ...      % Boxed background
                                                                                    'EdgeColor', 'k', ...            % Box border
                                                                                    'Margin', 4);                    % Padding inside box
                                    plot(m(j), result_nu_hat(i,j, k, l,f,g, e, o, 2), ...
                                    'p', ...                % Star marker
                                    'MarkerSize', 30, ...   % Large size
                                    'MarkerEdgeColor', 'k', ...
                                    'MarkerFaceColor', 'none', ... % Hollow
                                    'LineWidth', 3);
                                end

                            end
                        end
                    end


                end
            end

            xlabel('$\bar{m}$', 'Interpreter', 'latex', 'FontSize', 20);
            ylabel('$\hat{\nu}_{h}$', 'Interpreter', 'latex', 'FontSize', 20);
            ylim([0, 1]);
            grid on;
            % yline(0.0, 'k--', 'LineWidth', 1.5);
            legend(legend_handles, legend_labels, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best'); % Set legend font size
            title('Type H Reserve Prob $\hat{\nu}_{h}=1-\hat{\mu}_{hh}-\hat{\mu}_{lh}  $', 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
            hold off;

            sgtitle(sprintf('$\\gamma = %.2f, \\epsilon=%.2f, r=%.2f$', gamma(i), epsilon(e), r(l)), 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
            set(gcf, 'Position', [100, 100, 2200, 1000]); % Set figure size
            saveas(gcf, sprintf ([figure_save_path, '/ReserveProbability_gamma=%.2f_epsilon=%.2f_r=%.2f.png'], gamma(i), epsilon(e), r(l)));
            
        end
    end
end





%%
color_set = {'r', 'g', 'b', 'm', 'c', 'k'};  % one is orange


for i = 1:length(gamma)
    for e = 1:length(epsilon)
        for l = 1:length(r)

            figure;

            subplot(1, 2, 1);
            hold on;
            legend_handles = []; % Initialize an array to store scatter handles for the legend
            legend_labels = {};  % Initialize legend labels

            legend_num =0;
            for g = 1:length(alpha)
                for o = 1:length(beta)
                    legend_num = legend_num + 1;
                    color = color_set{legend_num};
                    legend_labels{legend_num} = sprintf('$\\alpha = %.2f, \\beta=%.2f$', alpha(g), beta(o)); % Add legend label for each alpha

                    % Create a dummy scatter plot for the legend
                    h = scatter(NaN, NaN, pointsize, color, 'filled');
                    legend_handles = [legend_handles, h];

                    for f = 1:length(theta)
                        for k = 1:length(phi_h)
                            for j = 1:length(m)

                                if result_flag(i,j, k, l,f,g, e, o, 1) == 2
                                    scatter(m(j), result_delta_hat(i,j, k, l,f,g, e, o, 1), pointsize, color, 'filled');
                                elseif result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                    scatter(m(j), result_delta_hat(i,j, k, l,f,g, e, o, 1), pointsize, 'o', 'MarkerEdgeColor', color, 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
                                elseif result_flag(i,j, k, l,f,g, e, o, 1) == 0
                                    scatter(m(j), result_delta_hat(i,j, k, l,f,g, e, o, 1), pointsize, 'x', 'MarkerEdgeColor', color, 'LineWidth', 1.5);
                                end
                                % if result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 1), pointsize, color, 'filled');
                                % % elseif result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                % %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 1), pointsize, 'o', 'MarkerEdgeColor', color);
                                % elseif result_flag(i,j, k, l,f,g, e, o, 1) == 0
                                %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 1), pointsize, 'x', 'MarkerEdgeColor', color, 'LineWidth', 1.5);
                                % end

                                if selector(gamma(i), m(j), phi_h(k), r(l), theta(f), alpha(g), epsilon(e), beta(o)) == 1
                                    h_txt = text(m(j), result_delta_hat(i,j, k, l,f,g, e, o, 1) - 0.2, ...
                                                                                    'Scenario 1', ...
                                                                                    'Color', 'k', ...
                                                                                    'FontSize', 20, ...
                                                                                    'Interpreter', 'latex', ...
                                                                                    'BackgroundColor', 'w', ...      % Boxed background
                                                                                    'EdgeColor', 'k', ...            % Box border
                                                                                    'Margin', 4);                    % Padding inside box
                                    plot(m(j), result_delta_hat(i,j, k, l,f,g, e, o, 1), ...
                                    'p', ...                % Star marker
                                    'MarkerSize', 30, ...   % Large size
                                    'MarkerEdgeColor', 'k', ...
                                    'MarkerFaceColor', 'none', ... % Hollow
                                    'LineWidth', 3);
                                end

                                if selector(gamma(i), m(j), phi_h(k), r(l), theta(f), alpha(g), epsilon(e), beta(o)) == 2
                                    h_txt = text(m(j), result_delta_hat(i,j, k, l,f,g, e, o, 1) - 0.4, ...
                                                                                    'Scenario 2', ...
                                                                                    'Color', 'k', ...
                                                                                    'FontSize', 20, ...
                                                                                    'Interpreter', 'latex', ...
                                                                                    'BackgroundColor', 'w', ...      % Boxed background
                                                                                    'EdgeColor', 'k', ...            % Box border
                                                                                    'Margin', 4);                    % Padding inside box
                                    plot(m(j), result_delta_hat(i,j, k, l,f,g, e, o, 1), ...
                                    'p', ...                % Star marker
                                    'MarkerSize', 30, ...   % Large size
                                    'MarkerEdgeColor', 'k', ...
                                    'MarkerFaceColor', 'none', ... % Hollow
                                    'LineWidth', 3);
                                end
                            end
                        end
                    end

                end
            end

            xlabel('$\bar{m}$', 'Interpreter', 'latex', 'FontSize', 20);
            ylabel('$\hat{\delta}_{l}$', 'Interpreter', 'latex', 'FontSize', 20);
            ylim([0, 1]);
            grid on;
            % yline(0.0, 'k--', 'LineWidth', 1.5);
            legend(legend_handles, legend_labels, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best'); % Set legend font size
            title('Type L External Innovation $\hat{\delta}_{l} $', 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
            hold off;


            subplot(1, 2, 2);
            hold on;
            legend_handles = []; % Initialize an array to store scatter handles for the legend
            legend_labels = {};  % Initialize legend labels
            legend_num =0;

            for g = 1:length(alpha)
                for o = 1:length(beta)
                    legend_num = legend_num + 1;

                    color = color_set{legend_num};
                    legend_labels{legend_num} = sprintf('$\\alpha = %.2f, \\beta=%.2f$', alpha(g), beta(o)); % Add legend label for each alpha

                    % Create a dummy scatter plot for the legend
                    h = scatter(NaN, NaN, pointsize, color, 'filled');
                    legend_handles = [legend_handles, h];

                    for f = 1:length(theta)
                        for k = 1:length(phi_h)
                            for j = 1:length(m)

                                
                                if result_flag(i,j, k, l,f,g, e, o, 1) == 2
                                    scatter(m(j), result_delta_hat(i,j, k, l,f,g, e, o, 2), pointsize, color, 'filled');
                                elseif result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                    scatter(m(j), result_delta_hat(i,j, k, l,f,g, e, o, 2), pointsize, 'o', 'MarkerEdgeColor', color, 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
                                elseif result_flag(i,j, k, l,f,g, e, o, 1) == 0
                                    scatter(m(j), result_delta_hat(i,j, k, l,f,g, e, o, 2), pointsize, 'x', 'MarkerEdgeColor', color, 'LineWidth', 1.5);
                                end

                                % if result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 2), pointsize, color, 'filled');
                                % % elseif result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                % %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 2), pointsize, 'o', 'MarkerEdgeColor', color, 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
                                % elseif result_flag(i,j, k, l,f,g, e, o, 1) == 0
                                %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 2), pointsize, 'x', 'MarkerEdgeColor', color, 'LineWidth', 1.5);
                                % end

                                if selector(gamma(i), m(j), phi_h(k), r(l), theta(f), alpha(g), epsilon(e), beta(o)) == 1
                                    h_txt = text(m(j), result_delta_hat(i,j, k, l,f,g, e, o, 2) - 0.2, ...
                                                                                    'Scenario 1', ...
                                                                                    'Color', 'k', ...
                                                                                    'FontSize', 20, ...
                                                                                    'Interpreter', 'latex', ...
                                                                                    'BackgroundColor', 'w', ...      % Boxed background
                                                                                    'EdgeColor', 'k', ...            % Box border
                                                                                    'Margin', 4);                    % Padding inside box
                                    plot(m(j), result_delta_hat(i,j, k, l,f,g, e, o, 2), ...
                                    'p', ...                % Star marker
                                    'MarkerSize', 30, ...   % Large size
                                    'MarkerEdgeColor', 'k', ...
                                    'MarkerFaceColor', 'none', ... % Hollow
                                    'LineWidth', 3);
                                end

                                if selector(gamma(i), m(j), phi_h(k), r(l), theta(f), alpha(g), epsilon(e), beta(o)) == 2
                                    h_txt = text(m(j), result_delta_hat(i,j, k, l,f,g, e, o, 2) - 0.4, ...
                                                                                    'Scenario 2', ...
                                                                                    'Color', 'k', ...
                                                                                    'FontSize', 20, ...
                                                                                    'Interpreter', 'latex', ...
                                                                                    'BackgroundColor', 'w', ...      % Boxed background
                                                                                    'EdgeColor', 'k', ...            % Box border
                                                                                    'Margin', 4);                    % Padding inside box
                                    plot(m(j), result_delta_hat(i,j, k, l,f,g, e, o, 2), ...
                                    'p', ...                % Star marker
                                    'MarkerSize', 30, ...   % Large size
                                    'MarkerEdgeColor', 'k', ...
                                    'MarkerFaceColor', 'none', ... % Hollow
                                    'LineWidth', 3);
                                end

                            end
                        end
                    end


                end
            end

            xlabel('$\bar{m}$', 'Interpreter', 'latex', 'FontSize', 20);
            ylabel('$\hat{\delta}_{h}$', 'Interpreter', 'latex', 'FontSize', 20);
            ylim([0, 1]);
            grid on;
            % yline(0.0, 'k--', 'LineWidth', 1.5);
            legend(legend_handles, legend_labels, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best'); % Set legend font size
            title('Type H External Innovation $\hat{\delta}_{h} $', 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
            hold off;

            sgtitle(sprintf('$\\gamma = %.2f, \\epsilon=%.2f, r=%.2f$', gamma(i), epsilon(e), r(l)), 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
            set(gcf, 'Position', [100, 100, 2200, 1000]); % Set figure size
            saveas(gcf, sprintf ([figure_save_path, '/ExternalInno_gamma=%.2f_epsilon=%.2f_r=%.2f.png'], gamma(i), epsilon(e), r(l)));
            
        end
    end
end






%%
color_set = {'r', 'g', 'b', 'm', 'c', 'k'};  % one is orange


for i = 1:length(gamma)
    for e = 1:length(epsilon)
        for l = 1:length(r)

            figure;

            subplot(1, 2, 1);
            hold on;
            legend_handles = []; % Initialize an array to store scatter handles for the legend
            legend_labels = {};  % Initialize legend labels

            legend_num =0;
            for g = 1:length(alpha)
                for o = 1:length(beta)
                    legend_num = legend_num + 1;
                    color = color_set{legend_num};
                    legend_labels{legend_num} = sprintf('$\\alpha = %.2f, \\beta=%.2f$', alpha(g), beta(o)); % Add legend label for each alpha

                    % Create a dummy scatter plot for the legend
                    h = scatter(NaN, NaN, pointsize, color, 'filled');
                    legend_handles = [legend_handles, h];

                    for f = 1:length(theta)
                        for k = 1:length(phi_h)
                            for j = 1:length(m)

                                if result_flag(i,j, k, l,f,g, e, o, 1) == 2
                                    scatter(m(j), result_lambda_hat(i,j, k, l,f,g, e, o, 1), pointsize, color, 'filled');
                                elseif result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                    scatter(m(j), result_lambda_hat(i,j, k, l,f,g, e, o, 1), pointsize, 'o', 'MarkerEdgeColor', color, 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
                                elseif result_flag(i,j, k, l,f,g, e, o, 1) == 0
                                    scatter(m(j), result_lambda_hat(i,j, k, l,f,g, e, o, 1), pointsize, 'x', 'MarkerEdgeColor', color, 'LineWidth', 1.5);
                                end
                                % if result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 1), pointsize, color, 'filled');
                                % % elseif result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                % %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 1), pointsize, 'o', 'MarkerEdgeColor', color);
                                % elseif result_flag(i,j, k, l,f,g, e, o, 1) == 0
                                %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 1), pointsize, 'x', 'MarkerEdgeColor', color, 'LineWidth', 1.5);
                                % end

                                if selector(gamma(i), m(j), phi_h(k), r(l), theta(f), alpha(g), epsilon(e), beta(o)) == 1
                                    h_txt = text(m(j), result_lambda_hat(i,j, k, l,f,g, e, o, 1) - 0.2, ...
                                                                                    'Scenario 1', ...
                                                                                    'Color', 'k', ...
                                                                                    'FontSize', 20, ...
                                                                                    'Interpreter', 'latex', ...
                                                                                    'BackgroundColor', 'w', ...      % Boxed background
                                                                                    'EdgeColor', 'k', ...            % Box border
                                                                                    'Margin', 4);                    % Padding inside box
                                    plot(m(j), result_lambda_hat(i,j, k, l,f,g, e, o, 1), ...
                                    'p', ...                % Star marker
                                    'MarkerSize', 30, ...   % Large size
                                    'MarkerEdgeColor', 'k', ...
                                    'MarkerFaceColor', 'none', ... % Hollow
                                    'LineWidth', 3);
                                end

                                if selector(gamma(i), m(j), phi_h(k), r(l), theta(f), alpha(g), epsilon(e), beta(o)) == 2
                                    h_txt = text(m(j), result_lambda_hat(i,j, k, l,f,g, e, o, 1) - 0.2, ...
                                                                                    'Scenario 2', ...
                                                                                    'Color', 'k', ...
                                                                                    'FontSize', 20, ...
                                                                                    'Interpreter', 'latex', ...
                                                                                    'BackgroundColor', 'w', ...      % Boxed background
                                                                                    'EdgeColor', 'k', ...            % Box border
                                                                                    'Margin', 4);                    % Padding inside box
                                    plot(m(j), result_lambda_hat(i,j, k, l,f,g, e, o, 1), ...
                                    'p', ...                % Star marker
                                    'MarkerSize', 30, ...   % Large size
                                    'MarkerEdgeColor', 'k', ...
                                    'MarkerFaceColor', 'none', ... % Hollow
                                    'LineWidth', 3);
                                end
                            end
                        end
                    end

                end
            end

            xlabel('$\bar{m}$', 'Interpreter', 'latex', 'FontSize', 20);
            ylabel('$\hat{\lambda}_{l}$', 'Interpreter', 'latex', 'FontSize', 20);
            ylim([0, 1]);
            grid on;
            % yline(0.0, 'k--', 'LineWidth', 1.5);
            legend(legend_handles, legend_labels, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best'); % Set legend font size
            title('Type L Internal Innovation $\hat{\lambda}_{l} $', 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
            hold off;


            subplot(1, 2, 2);
            hold on;
            legend_handles = []; % Initialize an array to store scatter handles for the legend
            legend_labels = {};  % Initialize legend labels
            legend_num =0;

            for g = 1:length(alpha)
                for o = 1:length(beta)
                    legend_num = legend_num + 1;

                    color = color_set{legend_num};
                    legend_labels{legend_num} = sprintf('$\\alpha = %.2f, \\beta=%.2f$', alpha(g), beta(o)); % Add legend label for each alpha

                    % Create a dummy scatter plot for the legend
                    h = scatter(NaN, NaN, pointsize, color, 'filled');
                    legend_handles = [legend_handles, h];

                    for f = 1:length(theta)
                        for k = 1:length(phi_h)
                            for j = 1:length(m)

                                
                                if result_flag(i,j, k, l,f,g, e, o, 1) == 2
                                    scatter(m(j), result_lambda_hat(i,j, k, l,f,g, e, o, 2), pointsize, color, 'filled');
                                elseif result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                    scatter(m(j), result_lambda_hat(i,j, k, l,f,g, e, o, 2), pointsize, 'o', 'MarkerEdgeColor', color, 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
                                elseif result_flag(i,j, k, l,f,g, e, o, 1) == 0
                                    scatter(m(j), result_lambda_hat(i,j, k, l,f,g, e, o, 2), pointsize, 'x', 'MarkerEdgeColor', color, 'LineWidth', 1.5);
                                end

                                % if result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 2), pointsize, color, 'filled');
                                % % elseif result_flag(i,j, k, l,f,g, e, o, 1) == 1
                                % %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 2), pointsize, 'o', 'MarkerEdgeColor', color, 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
                                % elseif result_flag(i,j, k, l,f,g, e, o, 1) == 0
                                %     scatter(m(j), result_s(i,j, k, l,f,g, e, o, 2), pointsize, 'x', 'MarkerEdgeColor', color, 'LineWidth', 1.5);
                                % end

                                if selector(gamma(i), m(j), phi_h(k), r(l), theta(f), alpha(g), epsilon(e), beta(o)) == 1
                                    h_txt = text(m(j), result_lambda_hat(i,j, k, l,f,g, e, o, 2) - 0.2, ...
                                                                                    'Scenario 1', ...
                                                                                    'Color', 'k', ...
                                                                                    'FontSize', 20, ...
                                                                                    'Interpreter', 'latex', ...
                                                                                    'BackgroundColor', 'w', ...      % Boxed background
                                                                                    'EdgeColor', 'k', ...            % Box border
                                                                                    'Margin', 4);                    % Padding inside box
                                    plot(m(j), result_lambda_hat(i,j, k, l,f,g, e, o, 2), ...
                                    'p', ...                % Star marker
                                    'MarkerSize', 30, ...   % Large size
                                    'MarkerEdgeColor', 'k', ...
                                    'MarkerFaceColor', 'none', ... % Hollow
                                    'LineWidth', 3);
                                end

                                if selector(gamma(i), m(j), phi_h(k), r(l), theta(f), alpha(g), epsilon(e), beta(o)) == 2
                                    h_txt = text(m(j), result_lambda_hat(i,j, k, l,f,g, e, o, 2) - 0.2, ...
                                                                                    'Scenario 2', ...
                                                                                    'Color', 'k', ...
                                                                                    'FontSize', 20, ...
                                                                                    'Interpreter', 'latex', ...
                                                                                    'BackgroundColor', 'w', ...      % Boxed background
                                                                                    'EdgeColor', 'k', ...            % Box border
                                                                                    'Margin', 4);                    % Padding inside box
                                    plot(m(j), result_lambda_hat(i,j, k, l,f,g, e, o, 2), ...
                                    'p', ...                % Star marker
                                    'MarkerSize', 30, ...   % Large size
                                    'MarkerEdgeColor', 'k', ...
                                    'MarkerFaceColor', 'none', ... % Hollow
                                    'LineWidth', 3);
                                end

                            end
                        end
                    end


                end
            end

            xlabel('$\bar{m}$', 'Interpreter', 'latex', 'FontSize', 20);
            ylabel('$\hat{\lambda}_{h}$', 'Interpreter', 'latex', 'FontSize', 20);
            ylim([0, 1]);
            grid on;
            % yline(0.0, 'k--', 'LineWidth', 1.5);
            legend(legend_handles, legend_labels, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best'); % Set legend font size
            title('Type H Internal Innovation $\hat{\lambda}_{h} $', 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
            hold off;

            sgtitle(sprintf('$\\gamma = %.2f, \\epsilon=%.2f, r=%.2f$', gamma(i), epsilon(e), r(l)), 'Interpreter', 'latex', 'FontSize', 20); % Set title font size
            set(gcf, 'Position', [100, 100, 2200, 1000]); % Set figure size
            saveas(gcf, sprintf ([figure_save_path, '/InternalInno_gamma=%.2f_epsilon=%.2f_r=%.2f.png'], gamma(i), epsilon(e), r(l)));
            
        end
    end
end

%%


color_set = {'r', 'g', 'b', 'm', 'c', 'k'};

for i = 1:length(gamma)
    for e = 1:length(epsilon)
        for l = 1:length(r)

            figure;

            legend_handles = [];
            legend_labels = {};

            for g = 1:length(alpha)
                subplot(1, 2, g);
                hold on;
                legend_handles = [];
                legend_labels = {};

                for o = 1:length(beta)
                    % Assign color based on g and o
                    if g == 1
                        color = color_set{o};      % o=1: 'r', o=2: 'g'
                    elseif g == 2
                        color = color_set{o+2};    % o=1: 'b', o=2: 'm'
                    end

                    legend_labels{o} = sprintf('$\\alpha = %.2f, \\beta=%.2f$', alpha(g), beta(o));
                    h = scatter(NaN, NaN, pointsize, color, 'filled');
                    legend_handles = [legend_handles, h];

                    for f = 1:length(theta)
                        for k = 1:length(phi_h)
                            plot(m, result_Expected_q(i, :, k, l, f, g, e, o, 1), 'Color', color, 'LineWidth', 2);
                            plot(m, result_Expected_q(i, :, k, l, f, g, e, o, 2), 'Color', color, 'LineWidth', 2, 'LineStyle', '--');
                        end
                    end
                end

                xlabel('$\bar{m}$', 'Interpreter', 'latex', 'FontSize', 20);
                ylabel('$E(q)$', 'Interpreter', 'latex', 'FontSize', 20);
                if gamma(i) == 0.1
                    ylim([0, 7]);
                elseif gamma(i) == 0.2
                    ylim([0, 5]);
                end
                grid on;
                legend(legend_handles, legend_labels, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best');
                title('Expected Quality: Solid for Type L, Dashed for Type H', 'Interpreter', 'latex', 'FontSize', 20);
                hold off;
            end
            sgtitle(sprintf('$\\gamma = %.2f, \\epsilon=%.2f, r=%.2f$', gamma(i), epsilon(e), r(l)), 'Interpreter', 'latex', 'FontSize', 20); % Set title font size

            set(gcf, 'Position', [100, 100, 2200, 1000]);
            saveas(gcf, sprintf([figure_save_path, '/q_expected_gamma=%.2f_epsilon=%.2f_r=%.2f.png'], gamma(i), epsilon(e), r(l)));
        end
    end
end


%%


color_set = {'r', 'g', 'b', 'm', 'c', 'k'};

for i = 1:length(gamma)
    for e = 1:length(epsilon)
        for l = 1:length(r)

            figure;

            legend_handles = [];
            legend_labels = {};

            for g = 1:length(alpha)
                subplot(1, 2, g);
                hold on;
                legend_handles = [];
                legend_labels = {};

                for o = 1:length(beta)
                    % Assign color based on g and o
                    if g == 1
                        color = color_set{o};      % o=1: 'r', o=2: 'g'
                    elseif g == 2
                        color = color_set{o+2};    % o=1: 'b', o=2: 'm'
                    end

                    legend_labels{o} = sprintf('$\\alpha = %.2f, \\beta=%.2f$', alpha(g), beta(o));
                    h = scatter(NaN, NaN, pointsize, color, 'filled');
                    legend_handles = [legend_handles, h];

                    for f = 1:length(theta)
                        for k = 1:length(phi_h)
                            plot(m, result_10quantile_q(i, :, k, l, f, g, e, o, 1), 'Color', color, 'LineWidth', 2);
                            plot(m, result_10quantile_q(i, :, k, l, f, g, e, o, 2), 'Color', color, 'LineWidth', 2, 'LineStyle', '--');
                        end
                    end
                end

                xlabel('$\bar{m}$', 'Interpreter', 'latex', 'FontSize', 20);
                ylabel('$E(q)$', 'Interpreter', 'latex', 'FontSize', 20);
                if gamma(i) == 0.1
                    ylim([0, 1.5]);
                elseif gamma(i) == 0.2
                    ylim([0, 1.5]);
                end
                grid on;
                legend(legend_handles, legend_labels, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best');
                title('10$\%$ Quantile: Solid for Type L, Dashed for Type H', 'Interpreter', 'latex', 'FontSize', 20);
                hold off;
            end
            sgtitle(sprintf('$\\gamma = %.2f, \\epsilon=%.2f, r=%.2f$', gamma(i), epsilon(e), r(l)), 'Interpreter', 'latex', 'FontSize', 20); % Set title font size

            set(gcf, 'Position', [100, 100, 2200, 1000]);
            saveas(gcf, sprintf([figure_save_path, '/q_10percent_gamma=%.2f_epsilon=%.2f_r=%.2f.png'], gamma(i), epsilon(e), r(l)));
        end
    end
end




%%


color_set = {'r', 'g', 'b', 'm', 'c', 'k'};

for i = 1:length(gamma)
    for e = 1:length(epsilon)
        for l = 1:length(r)

            figure;

            legend_handles = [];
            legend_labels = {};

            for g = 1:length(alpha)
                subplot(1, 2, g);
                hold on;
                legend_handles = [];
                legend_labels = {};

                for o = 1:length(beta)
                    % Assign color based on g and o
                    if g == 1
                        color = color_set{o};      % o=1: 'r', o=2: 'g'
                    elseif g == 2
                        color = color_set{o+2};    % o=1: 'b', o=2: 'm'
                    end

                    legend_labels{o} = sprintf('$\\alpha = %.2f, \\beta=%.2f$', alpha(g), beta(o));
                    h = scatter(NaN, NaN, pointsize, color, 'filled');
                    legend_handles = [legend_handles, h];

                    for f = 1:length(theta)
                        for k = 1:length(phi_h)
                            plot(m, result_50quantile_q(i, :, k, l, f, g, e, o, 1), 'Color', color, 'LineWidth', 2);
                            plot(m, result_50quantile_q(i, :, k, l, f, g, e, o, 2), 'Color', color, 'LineWidth', 2, 'LineStyle', '--');
                        end
                    end
                end

                xlabel('$\bar{m}$', 'Interpreter', 'latex', 'FontSize', 20);
                ylabel('$q$', 'Interpreter', 'latex', 'FontSize', 20);
                if gamma(i) == 0.1
                    ylim([0, 5.5]);
                elseif gamma(i) == 0.2
                    ylim([0, 5.5]);
                end
                grid on;
                legend(legend_handles, legend_labels, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best');
                title('50$\%$ Quantile: Solid for Type L, Dashed for Type H', 'Interpreter', 'latex', 'FontSize', 20);
                hold off;
            end
            sgtitle(sprintf('$\\gamma = %.2f, \\epsilon=%.2f, r=%.2f$', gamma(i), epsilon(e), r(l)), 'Interpreter', 'latex', 'FontSize', 20); % Set title font size

            set(gcf, 'Position', [100, 100, 2200, 1000]);
            saveas(gcf, sprintf([figure_save_path, '/q_50percent_gamma=%.2f_epsilon=%.2f_r=%.2f.png'], gamma(i), epsilon(e), r(l)));
        end
    end
end




%%


color_set = {'r', 'g', 'b', 'm', 'c', 'k'};

for i = 1:length(gamma)
    for e = 1:length(epsilon)
        for l = 1:length(r)

            figure;

            legend_handles = [];
            legend_labels = {};

            for g = 1:length(alpha)
                subplot(1, 2, g);
                hold on;
                legend_handles = [];
                legend_labels = {};

                for o = 1:length(beta)
                    % Assign color based on g and o
                    if g == 1
                        color = color_set{o};      % o=1: 'r', o=2: 'g'
                    elseif g == 2
                        color = color_set{o+2};    % o=1: 'b', o=2: 'm'
                    end

                    legend_labels{o} = sprintf('$\\alpha = %.2f, \\beta=%.2f$', alpha(g), beta(o));
                    h = scatter(NaN, NaN, pointsize, color, 'filled');
                    legend_handles = [legend_handles, h];

                    for f = 1:length(theta)
                        for k = 1:length(phi_h)
                            plot(m, result_90quantile_q(i, :, k, l, f, g, e, o, 1), 'Color', color, 'LineWidth', 2);
                            plot(m, result_90quantile_q(i, :, k, l, f, g, e, o, 2), 'Color', color, 'LineWidth', 2, 'LineStyle', '--');
                        end
                    end
                end

                xlabel('$\bar{m}$', 'Interpreter', 'latex', 'FontSize', 20);
                ylabel('$q$', 'Interpreter', 'latex', 'FontSize', 20);
                if gamma(i) == 0.1
                    ylim([1, 12]);
                elseif gamma(i) == 0.2
                    ylim([1, 12]);
                end
                grid on;
                legend(legend_handles, legend_labels, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best');
                title('90$\%$ Quantile: Solid for Type L, Dashed for Type H', 'Interpreter', 'latex', 'FontSize', 20);
                hold off;
            end
            sgtitle(sprintf('$\\gamma = %.2f, \\epsilon=%.2f, r=%.2f$', gamma(i), epsilon(e), r(l)), 'Interpreter', 'latex', 'FontSize', 20); % Set title font size

            set(gcf, 'Position', [100, 100, 2200, 1000]);
            saveas(gcf, sprintf([figure_save_path, '/q_90percent_gamma=%.2f_epsilon=%.2f_r=%.2f.png'], gamma(i), epsilon(e), r(l)));
        end
    end
end








%%


p_scenario1.gamma = 0.2;
p_scenario1.m_bar = 0.35;
p_scenario1.phi_h = 1.05;
p_scenario1.r = 0.15;
p_scenario1.theta = 5;
p_scenario1.alpha = 0.9;
p_scenario1.epsilon = 0.0;
p_scenario1.beta = 8;
p_scenario1.sigma=3;


p_scenario2.gamma = 0.2;
p_scenario2.m_bar = 0.35;
p_scenario2.phi_h = 1.05;
p_scenario2.r = 0.15;
p_scenario2.theta = 5;
p_scenario2.alpha = 0.9;
p_scenario2.epsilon = 1.0;
p_scenario2.beta = 8;
p_scenario2.sigma=3;


selector_plot(p_scenario1 );


selector_plot(p_scenario2);










%%

file_name = 'epsilon=1.00, gamma=0.20, sigma=3, phi_h=1.05, m_bar=0.90, r=0.15, alpha=0.90, theta=5.00, beta=8.00';
% data = load(['./data_beta/', file_name, '.mat'] , 'p', 'np' , 'eqm_save', 'iter_history');

% try
data = load(['./data_beta/', file_name, '.mat'] , 'p', 'np' , 'eqm_save', 'iter_history');
% catch ME
%     fprintf('[LOAD ERROR] Failed to load: %s.mat\nError: %s\n\n', file_name, ME.message);
%     % bad_files{end+1} = file_name;
%     % continue;
% end

eqm = data.eqm_save;
np = data.np;
p = data.p;

% eqm.calE
% eqm.calL

coef_l = eqm.omega_tilde(1,end) / (np.q(end)^(-p.theta - 1));
coef_h = eqm.omega_tilde(2,end) / (np.q(end)^(-p.theta - 1));
eqm.omega_tilde_interp_l = @(q) (q > max(np.q)) .* coef_l .* q.^(-p.theta-1) + ...
    (q >= min(np.q) & q <= max(np.q)) .* max(0, interp1(np.q, eqm.omega_tilde(1,:), q, 'spline', 'extrap')); % outside is set  for convergence of calE
eqm.omega_tilde_interp_h = @(q) (q > max(np.q)) .* coef_h .* q.^(-p.theta-1) + ...
    (q >= min(np.q) & q <= max(np.q)) .* max(0, interp1(np.q, eqm.omega_tilde(2,:), q, 'spline', 'extrap'));
eqm.lambda_interp_l = @(q) max(0, min(1, interp1(np.q, eqm.lambda(1,:), q, 'spline', 'extrap')));
eqm.lambda_interp_h = @(q) max(0, min(1, interp1(np.q, eqm.lambda(2,:), q, 'spline', 'extrap')));

% eqm.calE(1)        = integral(  @(q)  q.^(p.sigma - 1) .* eqm.omega_tilde_interp_l(q) , np.q_min, Inf);  % need more than 100 like 200 poitns to achieve accuracy
% eqm.calE(2)        = integral(  @(q)  q.^(p.sigma - 1) .* eqm.omega_tilde_interp_h(q) , np.q_min, Inf);

% eqm.calL(1) = integral(  @(q)  q.^(p.sigma - 1) .* eqm.lambda_interp_l(q) .* eqm.omega_tilde_interp_l(q) , np.q_min, Inf);
% eqm.calL(2) = integral(  @(q)  q.^(p.sigma - 1) .* eqm.lambda_interp_h(q) .* eqm.omega_tilde_interp_h(q) , np.q_min, Inf);



eqm.calE(1)        = integral(  @(q)  q.^(p.sigma - 1) .* eqm.omega_tilde_interp_l(q) , np.q_min, 3*  np.q_max);  % need more than 100 like 200 poitns to achieve accuracy
eqm.calE(2)        = integral(  @(q)  q.^(p.sigma - 1) .* eqm.omega_tilde_interp_h(q) , np.q_min, 3*  np.q_max);  % need

eqm.calL(1) = integral(  @(q)  q.^(p.sigma - 1) .* eqm.lambda_interp_l(q) .* eqm.omega_tilde_interp_l(q) , np.q_min, np.q_max);  % need
eqm.calL(2) = integral(  @(q)  q.^(p.sigma - 1) .* eqm.lambda_interp_h(q) .* eqm.omega_tilde_interp_h(q) , np.q_min, np.q_max);  % need


% eqm.calE(1) = trapz(np.q, np.q.^(p.sigma - 1) .* eqm.omega(1,:));
% eqm.calE(2) = trapz(np.q, np.q.^(p.sigma - 1) .* eqm.omega(2,:));

% eqm.calL(1) = trapz(np.q, np.q.^(p.sigma - 1) .* eqm.lambda(1,:) .*  eqm.omega_tilde(1,:));
% eqm.calL(2) = trapz(np.q, np.q.^(p.sigma - 1) .* eqm.lambda(2,:) .*  eqm.omega_tilde(2,:));

eqm.Q_type_sigma1(1) = eqm.calE(1) / ( eqm.s(1) * p.phi(1)^(p.sigma - 1)  );
eqm.Q_type_sigma1(2) = eqm.calE(2) / ( eqm.s(2) * p.phi(2)^(p.sigma - 1)  );
eqm.Q_sigma1 = eqm.Q_type_sigma1(1) + eqm.Q_type_sigma1(2);



eqm.Q_type_sigma1


constant = 5;

coef_l = eqm.omega_tilde(1,end) / (np.q(end)^(-p.theta - 1));
coef_h = eqm.omega_tilde(2,end) / (np.q(end)^(-p.theta - 1));
eqm.omega_tilde_interp_l = @(q) (q > max(np.q)) .* coef_l .* q.^(-p.theta-1) + ...
    (q >= min(np.q) & q <= max(np.q)) .* max(0, interp1(np.q, eqm.omega_tilde(1,:), q, 'spline', 'extrap')); % outside is set  for convergence of calE
eqm.omega_tilde_interp_h = @(q) (q > max(np.q)) .* coef_h .* q.^(-p.theta-1) + ...
    (q >= min(np.q) & q <= max(np.q)) .* max(0, interp1(np.q, eqm.omega_tilde(2,:), q, 'spline', 'extrap'));
eqm.lambda_interp_l = @(q) max(0, min(1, interp1(np.q, eqm.lambda(1,:), q, 'spline', 'extrap')));
eqm.lambda_interp_h = @(q) max(0, min(1, interp1(np.q, eqm.lambda(2,:), q, 'spline', 'extrap')));

% eqm.calE(1)        = integral(  @(q)  q.^(p.sigma - 1) .* eqm.omega_tilde_interp_l(q) , np.q_min, Inf);  % need more than 100 like 200 poitns to achieve accuracy
% eqm.calE(2)        = integral(  @(q)  q.^(p.sigma - 1) .* eqm.omega_tilde_interp_h(q) , np.q_min, Inf);

% eqm.calL(1) = integral(  @(q)  q.^(p.sigma - 1) .* eqm.lambda_interp_l(q) .* eqm.omega_tilde_interp_l(q) , np.q_min, Inf);
% eqm.calL(2) = integral(  @(q)  q.^(p.sigma - 1) .* eqm.lambda_interp_h(q) .* eqm.omega_tilde_interp_h(q) , np.q_min, Inf);



eqm.calE(1)        = integral(  @(q)  q.^(p.sigma - 1) .* eqm.omega_tilde_interp_l(q) , np.q_min, constant* np.q_max);  % need more than 100 like 200 poitns to achieve accuracy
eqm.calE(2)        = integral(  @(q)  q.^(p.sigma - 1) .* eqm.omega_tilde_interp_h(q) , np.q_min, constant* np.q_max);  % need

eqm.calL(1) = integral(  @(q)  q.^(p.sigma - 1) .* eqm.lambda_interp_l(q) .* eqm.omega_tilde_interp_l(q) , np.q_min, np.q_max);  % need
eqm.calL(2) = integral(  @(q)  q.^(p.sigma - 1) .* eqm.lambda_interp_h(q) .* eqm.omega_tilde_interp_h(q) , np.q_min, np.q_max);  % need


% eqm.calE(1) = trapz(np.q, np.q.^(p.sigma - 1) .* eqm.omega(1,:));
% eqm.calE(2) = trapz(np.q, np.q.^(p.sigma - 1) .* eqm.omega(2,:));

% eqm.calL(1) = trapz(np.q, np.q.^(p.sigma - 1) .* eqm.lambda(1,:) .*  eqm.omega_tilde(1,:));
% eqm.calL(2) = trapz(np.q, np.q.^(p.sigma - 1) .* eqm.lambda(2,:) .*  eqm.omega_tilde(2,:));

eqm.Q_type_sigma1(1) = eqm.calE(1) / ( eqm.s(1) * p.phi(1)^(p.sigma - 1)  );
eqm.Q_type_sigma1(2) = eqm.calE(2) / ( eqm.s(2) * p.phi(2)^(p.sigma - 1)  );
eqm.Q_sigma1 = eqm.Q_type_sigma1(1) + eqm.Q_type_sigma1(2);

eqm.Q_type_sigma1


