clc;
clear;

% Parameters
L = 1.5;            % Domain length
M = 75;             % Number of cells
dx = L / M;         % Grid spacing
CFL = 0.8;          % CFL coefficient
dt = CFL * dx;      % Time step size
T = 0.5;              % Total simulation time

% Define the initial condition
initial_condition = @(x) (-0.5 * (x <= 0.5)) + (1 * (x > 0.5 & x <= 1));

% Initialize variables
x = linspace(0, L, M);  % Spatial grid
u = initial_condition(x); % Initial solution
exact = @(x) (-0.5 * (x <= 0.25)) + (((2*x)-1) .* (x >= 0.25 & x <= 1)) + (1 * (x >= 1 & x <= 1.25)) + (0 * (x >=1.25));
uexact = exact(x);

% Plot initial condition
figure;
plot(x, uexact, 'b-', 'LineWidth', 1.5);
hold on;
xlim([0, L]);
ylim([-1, 2]);
xlabel('x');
ylabel('u');
title('Numerical and Exact Solutions');
grid on;

% Main loop for time evolution using Lax-Friedrichs scheme
t = 0;
while t < T
    % Compute fluxes at cell interfaces using Lax-Friedrichs scheme
    fluxes = zeros(1, M-1);
    for i = 1:M-1
        % Calculate left and right states at cell interfaces
        u_L = u(i);
        u_R = u(i+1);
        
        % Calculate flux at the interface using Lax-Friedrichs method
        f_L = 0.5 * (u_L^2);
        f_R = 0.5 * (u_R^2);
        
        % Compute numerical flux using Lax-Friedrichs method
        fluxes(i) = 0.5 * (f_L + f_R) - 0.5 * CFL * (u_R - u_L);
    end
    
    % Apply boundary conditions
    u_new = zeros(size(u));
    u_new(1) = u(1); % Left boundary condition
    u_new(end) = 0;  % Right boundary condition (y = 0 for x > 1.5)
    
    % Update solution using Lax-Friedrichs scheme
    for i = 2:M-1
        u_new(i) = u(i) - (dt / dx) * (fluxes(i) - fluxes(i-1));
    end
    
    % Update time
    t = t + dt;
    
    % Plot updated solution
    
    title(sprintf('t = %.3f', t));
    
    
    % Update solution for next iteration
    u = u_new;
    
    % Plot the exact solution
    exact_solution = @(x, t) (-0.5 * (x <= 0.5 + t)) + (1 * (x > 0.5 + t & x <= 1 + t));
  
end
plot(x, u_new, 'ro', 'LineWidth', 1.5);
legend('Initial condition', 'Numerical Solution');