clc;
clear;

% Parameters
num_points = 100;
CFL = 0.8; % Courant number
a = 1;
end_time = 1;

% Define the moving initial condition function
% Modified to ensure periodic wrapping
initial_condition = @(x, t) 0 * (x <= 0.3 +t)+ 1 * (x > 0.3+t & x < 0.7+t) + 0 * (x >= 0.7+t);
x_values = linspace(0, 2, num_points + 1);
h = 1 / num_points;
k = CFL * h / a;

% Initialize variables and plot initial conditions
time = 0;
U = zeros(3, num_points + 1);
U = [initial_condition(x_values, time); initial_condition(x_values, time); initial_condition(x_values, time)];
U_temp = U;


% Plotting section
figure;
%plot(x_values, initial_condition(x_values, end_time), 'k--'); % Initial condition plot
hold on
theoretical_plot = plot(x_values, initial_condition(x_values, end_time), 'k-', 'LineWidth', 2);  % Dynamic theoretical plot
LW_plot = plot(x_values + a * end_time, U(1, :), 'bo');  % Lax-Wendroff plot
LF_plot = plot(x_values + a * end_time, U(2, :), 'r.');  % Lax-Friedrichs plot
WB_plot = plot(x_values + a * end_time, U(3, :), 'g^');  % Warming-Beam plot
hold off
xlim([0.5, 2.5]); % Adjust the x-axis limits accordingly
ylim([-0.1, 1.1]);
legend('Theoretical', 'Lax-Wendroff', 'Lax-Friedrichs', 'Warming-Beam');
title_text = title(sprintf('t = %0.3f', end_time));
xlabel('x + at');
ylabel('u');


% Coefficients for Warming-Beam scheme
lambda = CFL;
lambda1 = lambda - 1;
lambda2 = lambda - 2;

% Time evolution loop
while (time + k) < end_time
    for j = 3:num_points % Avoid indexing errors with Warming-Beam
        % Lax-Wendroff scheme
        U_temp(1, j) = U(1, j) - 0.5 * CFL * (U(1, j + 1) - U(1, j - 1)) + 0.5 * (CFL)^2 * (U(1, j + 1) - 2 * U(1, j) + U(1, j - 1));
        % Lax-Friedrichs scheme
        U_temp(2, j) = 0.5 * (U(2, j - 1) + U(2, j + 1)) - 0.5 * CFL * (U(2, j + 1) - U(2, j - 1));
        % Warming-Beam update
        U_temp(3, j) = (0.5 * lambda * lambda1 * U(3, j - 2)) - (lambda * lambda2 * U(3, j - 1)) + (0.5 * lambda2 * lambda1 * U(3, j));
        %U_temp(3, j) = U(3, j) - 0.5 * CFL * (U(3, j + 1) - U(3, j - 1)) + 0.5 * CFL^2 * (U(3, j - 1) - 2 * U(3, j) + U(3, j + 1));

    end
    
    % Apply periodic boundary conditions
    U_temp(:, 1) = U_temp(:, num_points);
    U_temp(:, 2) = U_temp(:, num_points - 1); % Handling for j=3 case in Warming-Beam
    U_temp(:, num_points + 1) = U_temp(:, 2);

    % Update variables for the next iteration
    U = U_temp;
    time = time + k;
    set(theoretical_plot, 'YData', initial_condition(x_values, time));  
    set(LW_plot, 'YData', U(1, :));
    set(LF_plot, 'YData', U(2, :));
    set(WB_plot, 'YData', U(3, :));
    set(title_text, 'String', sprintf('t = %0.3f', time));

    %initial_condition = @(x, t) 0 * (x <= 0.3 + t*0.5) + 1 * (x > 0.3 + t*0.5 & x < 0.7 + t*0.5) + 0 * (x >= 0.7 + t*0.5);
    x_values = linspace(0, 2, num_points + 1);


end