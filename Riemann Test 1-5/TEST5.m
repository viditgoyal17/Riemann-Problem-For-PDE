    %-----THIS IS THE PROGRAM TO FIND THE EXACT SOLUTION FOR TEST-1------------
%------- OF A RIEMANN PROBLEM AND HELPS TO FIND THE DENSITY,---------------
%----------PRESSURE, VELOCITY AND INTERNAL ENERGY USING -------------------
%--------------------NEWTON RAPHSON METHOD---------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--- Test 5 is made up of the right and left shocks emerging
%---from the solution to tests 3 and 4 respectively; its solution represents the
%---colision of these two strong shocks and consists of a left facing shock (travelling
%---very slowly to the right), a right travelling contact discontinuity and a right
%---travelling shock wave. Graph appearing on running the code shows
%---solution profiles.
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%

function TEST5()

    % Define and globalizing initial conditions
    global rho_l P_l u_l rho_r P_r u_r gamma mu a_l a_r
    
    % Initial conditions
    rho_l = 5.99924;
    P_l = 460.894;
    u_l = 19.5975;

    rho_r = 5.99242;
    P_r = 46.0950;
    u_r = -6.19633;

    gamma = 1.4;
    mu = sqrt((gamma - 1) / (gamma + 1));

    % Speed of sound
    a_l = sqrt(gamma * P_l / rho_l);
    a_r = sqrt(gamma * P_r / rho_r);

    % Call Test_5 function to obtain data
    time = 0.035;
    data = Test_5(time);

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

function [data] = Test_5(t)
    
    global rho_l P_l u_l rho_r P_r u_r gamma mu a_l a_r
    
    % Initial conditions
    x0 = 0;
    A_R = 2/((gamma+1)*rho_r);
B_R = mu*mu*P_r;
A_L = 2/((gamma+1)*rho_l);
B_L = mu*mu*P_l;

    % Test 5
    P_star = fzero(@FunctionTest5, 1241.21);
  v_star = u_l - (P_star-P_l)*power((A_L/(P_star+B_L)),0.5);
v_shock_r = u_r + a_r*power((((gamma+1)/(2*gamma))*(P_star/P_r)) + ((gamma-1)/(2*gamma)) ,0.5);
v_shock_l = u_l - a_l*power((((gamma+1)/(2*gamma))*(P_star/P_l)) + ((gamma-1)/(2*gamma)) ,0.5);
rho_star_R = rho_r*(( (P_star/P_r) + mu^2 )/(1 + mu*mu*(P_star/P_r)));
rho_star_L = rho_l*(( (P_star/P_l) + mu^2 )/(1 + mu*mu*(P_star/P_l)));

    % Key Positions
    x1 = x0 - v_shock_l * t;
    x2 = x0 + v_star * t;
    x3 = x0 + v_shock_r * t;

    disp("p*")
    disp(P_star);
    disp("u*")
    disp(v_star);
    disp("rho*L")
    disp(rho_star_L);
    disp("rho*R")
    disp(rho_star_R);

    % Start setting values
    n_points = 1000;    % Set by user
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
            % Solution before x1 (before left shock)
            data.rho(index) = rho_l;
            data.P(index) = P_l;
            data.u(index) = u_l;
        elseif (x1 <= data.x(index) && data.x(index) <= x2)
            % Solution between x1 and x2 (between left shock and contact)
            data.rho(index) = rho_star_L;
            data.P(index) = P_star;
            data.u(index) = v_star;
        elseif (x2 <= data.x(index) && data.x(index) <= x3)
            % Solution between x2 and x3 (between contact and right shock)
            data.rho(index) = rho_star_R;
            data.P(index) = P_star;
            data.u(index) = v_star;
        elseif x3 < data.x(index)
            % Solution after x3 (after right shock)
            data.rho(index) = rho_r;
            data.P(index) = P_r;
            data.u(index) = u_r;
        end
        data.e(index) = data.P(index) / ((gamma - 1) * data.rho(index));
    end
end

function y = FunctionTest5(P)
    global rho_l P_l u_l rho_r P_r u_r gamma mu a_l a_r
    
    % Test 5 Sod Problem
    A_R = 2 / ((gamma + 1) * rho_r);
    B_R = mu^2 * P_r;
    A_L = 2 / ((gamma + 1) * rho_l);
    B_L = mu^2 * P_l;
      
   y = (P - P_r)*power(((A_R/(B_R+P))),0.5) + (P - P_l)*power(((A_L/(B_L+P))),0.5) + u_r-u_l;
end