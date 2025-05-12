function np = NumericalParameters(func,p)

    np.q_max = 20;
    np.q_min = 0.01;
    np.q_num = 200;
    np.q_grid_order = 3;
    np.q = linspace(0,1, np.q_num);
    np.q = np.q.^np.q_grid_order; % q grid
    np.q = np.q_min + (np.q_max - np.q_min) * np.q; % q grid
    % np.dq = func.forward_diff(np.q); % q grid spacing
    % np.dq = np.q(2) - np.q(1); % q grid spacing

    % np.y_max = 10;
    % np.y_min = 1;
    % np.y_num = np.q_num;
    % np.y_grid_order = 1;
    % np.y = linspace(0,1,np.y_num);
    % np.y = np.y.^np.y_grid_order; % y grid
    % np.y = np.y_min + (np.y_max - np.y_min) * np.y; % y grid
    % % np.dy = func.forward_diff(np.y); % y grid spacing
    % np.dy = np.y(2) - np.y(1); % y grid spacing

    % np.x_max = 10;
    % np.x_min = p.alpha;
    % np.x_num = np.q_num;
    % np.x_grid_order = 1;
    % u = linspace(0,1,np.x_num); 
    % u = u.^np.x_grid_order; 
    % np.x = np.x_min + (np.x_max - np.x_min) * u;
    % % np.dx = func.forward_diff(np.x); % x grid spacing
    % np.dx = np.x(2) - np.x(1); % x grid spacing

    % np.iter_step_max = 1000; % max number of iterations
    % np.iter_diff_max = 1e-6; % max difference between iterations

    % np.f_X = p.theta * p.alpha^p.theta * np.x.^( -p.theta - 1); % f_X function)
    % np.f_Y = p.theta  * np.y.^( -p.theta - 1); % f_Y function

    % figure;
    % hold on;
    % plot(np.x, np.f_X, 'r', 'LineWidth', 2); % plot f_X function
    % plot(np.y, np.f_Y, 'b', 'LineWidth', 2); % plot f_Y function
    % legend('$f_X(x)$', '$f_Y(y)$', 'Interpreter', 'latex', 'Location', 'best');
    % title('Density Functions');
    % xlabel('x, y');
    % ylabel('Density');



end