%-----THIS IS THE PROGRAM TO FIND THE EXACT SOLUTION FOR TEST-1------------
%------- OF A RIEMANN PROBLEM AND HELPS TO FIND THE DENSITY,---------------
%----------PRESSURE, VELOCITY AND INTERNAL ENERGY USING -------------------
%--------------------NEWTON RAPHSON METHOD---------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%---Test 4 is the right half of the Woodward and Colella problem; its
%---solution contains a left shock, a contact discontinuity and a
%---right rarefaction, graph appearing on running the code shows solution profiles
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function TEST4()

    % Define and globalizing initial conditions
    global rho_l P_l u_l rho_r P_r u_r gamma mu a_l a_r
    
    % Initial conditions
    rho_l = 1;
    P_l = 0.01;
    u_l = 0;

    rho_r = 1;
    P_r = 100;
    u_r = 0;

    gamma = 1.4;
    mu = sqrt((gamma - 1) / (gamma + 1));

    % Speed of sound
    a_l = sqrt(gamma * P_l / rho_l);
    a_r = sqrt(gamma * P_r / rho_r);

    % Call Test_4 function to obtain data
    time = 0.035;
    data = Test_4(time);

    % Plot the data
    plotData(data);

end

function plotData(data)
    % Plotting
    figure,
    subplot(2, 2, 1),
    plot(data.x, data.rho, '--r', 'LineWidth', 2);
    xlabel('x (m)');
    ylabel('Density (kg/m^3)');
    title('Density vs Position');
    grid on;

    subplot(2, 2, 2),
    plot(data.x, data.P, '-.g', 'LineWidth', 2);
    xlabel('x (m)');
    ylabel('Pressure (Pa)');
    title('Pressure vs Position');
    grid on;

    subplot(2, 2, 3),
    plot(data.x, data.u, ':b', 'LineWidth', 2);
    xlabel('x (m)');
    ylabel('Velocity (m/s)');
    title('Velocity vs Position');
    grid on;

    subplot(2, 2, 4),
    plot(data.x, data.e, '-m', 'LineWidth', 2);
    xlabel('x (m)');
    ylabel('Specific Internal Energy (J/kg)');
    title('Internal Energy vs Position');
    grid on;
end

function [data] = Test_4(t)
    global rho_l P_l u_l rho_r P_r u_r gamma mu a_l a_r
    
    % Initial conditions
    x0 = 0;

    % Test 4
    P_star = fzero(@FunctionTest4, 46.4162);
    v_star = u_r - 2 * (a_r / (gamma - 1)) * (1 - power((P_star / P_r), ((gamma - 1) / (2 * gamma))));
    rho_star_R = rho_r * power((P_star / P_r), 1 / gamma);
    rho_star_L = rho_l * ((P_star / P_l) + mu^2) / (1 + mu^2 * (P_star / P_l));
    v_shock = v_star * ((rho_star_L / rho_l) / ((rho_star_L / rho_l) - 1));

    % Key Positions
    x1 = x0 + v_shock * t;
    x2 = x0 + v_star * t;
    x4 = x0 + a_r * t;
    % Determining x3
    c_3 = a_r + ((gamma - 1) / 2) * v_star;
    x3 = x0 + (v_star + c_3) * t;

    disp("p*")
    disp(P_star);
    disp("u*")
    disp(v_star);
    disp("rho*L")
    disp(rho_star_L);
    disp("rho*R")
    disp(rho_star_R);

    % Start setting values
    n_points = 1000000;    % Set by user
    % Boundaries (can be set)
    x_min = -0.5;
    x_max = 0.5;

    x = linspace(x_min, x_max, n_points);
    data.x = x';
    data.rho = zeros(n_points, 1);   % Density
    data.P = zeros(n_points, 1);     % Pressure
    data.u = zeros(n_points, 1);     % Velocity
    data.e = zeros(n_points, 1);     % Internal energy

    for index = 1:n_points
        if data.x(index) < x1
            % Solution before x1 (before shock)
            data.rho(index) = rho_l;
            data.P(index) = P_l;
            data.u(index) = u_l;
        elseif (x1 <= data.x(index) && data.x(index) <= x2)
            % Solution between x1 and x2 (shock to contact)
            data.rho(index) = rho_star_L;
            data.P(index) = P_star;
            data.u(index) = v_star;
        elseif (x2 <= data.x(index) && data.x(index) <= x3)
            % Solution between x2 and x3 (contact to fan)
            data.rho(index) = rho_star_R;
            data.P(index) = P_star;
            data.u(index) = v_star;
        elseif (x3 <= data.x(index) && data.x(index) <= x4)
            % Solution between x3 and x4 (fan)
            c = -mu^2 * (u_r - (data.x(index) / t)) + (1 - mu^2) * a_r; 
            data.rho(index) = rho_r * power((c / a_r), 2 / (gamma - 1));
            data.P(index) = P_r * power((data.rho(index) / rho_r), gamma);
            data.u(index) = (1 - mu^2) * ((data.x(index) / t) - a_r + ((gamma - 1) / 2) * u_r);
        elseif x4 < data.x(index)
            % Solution after x4 (after fan)
            data.rho(index) = rho_r;
            data.P(index) = P_r;
            data.u(index) = u_r;
        end
        data.e(index) = data.P(index) / ((gamma - 1) * data.rho(index));
    end
end

function y = FunctionTest4(P)
    global rho_l P_l u_l rho_r P_r u_r gamma mu a_l a_r
    
    % Test 4 Sod Problem
    A_R = 2 / ((gamma + 1) * rho_r);
    B_R = (mu^2 * P_r);
    A_L = 2 / ((gamma + 1) * rho_l);
    B_L = (mu^2 * P_l);
      
    y = (P-P_l)*power((A_L/(P+B_L)),0.5) + ((2*a_r)/(gamma-1))*(-1+power((P/P_r),((gamma-1)/(2*gamma)))) +u_r-u_l;
end