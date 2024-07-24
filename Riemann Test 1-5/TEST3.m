%-----THIS IS THE PROGRAM TO FIND THE EXACT SOLUTION FOR TEST-1------------
%------- OF A RIEMANN PROBLEM AND HELPS TO FIND THE DENSITY,---------------
%----------PRESSURE, VELOCITY AND INTERNAL ENERGY USING -------------------
%--------------------NEWTON RAPHSON METHOD---------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%---Test 3 is a very severe test problem, the solution of which contains a
%---left rarefaction, a contact and a right shock; this test is actually the left half
%---of the blast wave problem of Woodward and Colella, graph appearing on running
%---the code shows solution profiles
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function TEST3()

    global rho_l P_l u_l rho_r P_r u_r gamma mu a_l a_r
    

    rho_l = 1;
    P_l = 1000;
    u_l = 0;

    rho_r = 1;
    P_r = 0.01;
    u_r = 0;

    gamma = 1.4;
    mu = sqrt((gamma - 1) / (gamma + 1));


    a_l = sqrt(gamma * P_l / rho_l);
    a_r = sqrt(gamma * P_r / rho_r);
    time = 0.012;
    data = Test_3(time);
    plotData(data);

end

function plotData(data)

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

function [data] = Test_3(t)

    global rho_l P_l u_l rho_r P_r u_r gamma mu a_l a_r
    

    x0 = 0;

    
    P_star = fzero(@FunctionTest3, 500.005);
    v_star = u_l + 2 * (a_l / (gamma - 1)) * (1 - power((P_star / P_l), ((gamma - 1) / (2 * gamma))));
    rho_star_R = rho_r * ((P_star / P_r) + mu^2) / (1 + mu^2 * (P_star / P_r));
    v_shock = v_star * ((rho_star_R / rho_r) / ((rho_star_R / rho_r) - 1));
    rho_star_L = rho_l * power((P_star / P_l), 1 / gamma);

    x1 = x0 - a_l * t; %Here u_l is zero so it was removed from the program to improve efficiency
    x3 = x0 + v_star * t;
    x4 = x0 + v_shock * t;
    c_2 = a_l - ((gamma - 1) / 2) * v_star;
    x2 = x0 + (v_star - c_2) * t;

    disp("p*")
    disp(P_star);
    disp("u*")
    disp(v_star);
    disp("rho*L")
    disp(rho_star_L);
    disp("rho*R")
    disp(rho_star_R);


    n_points = 1000000;    
    x_min = -0.5;
    x_max = 0.5;

    x = linspace(x_min, x_max, n_points);
    data.x = x';
    data.rho = zeros(n_points, 1);   
    data.P = zeros(n_points, 1);     
    data.u = zeros(n_points, 1);     
    data.e = zeros(n_points, 1);     

    for index = 1:n_points
        if data.x(index) < x1
            % Solution before x1 (before fan
            data.rho(index) = rho_l;
            data.P(index) = P_l;
            data.u(index) = u_l;
        elseif (x1 <= data.x(index) && data.x(index) <= x2)
            % Solution between x1 and x2 (fan)
            c = mu^2 * ((x0 - data.x(index)) / t) + (1 - mu^2) * a_l; 
            data.rho(index) = rho_l * power((c / a_l), 2 / (gamma - 1));
            data.P(index) = P_l * power((data.rho(index) / rho_l), gamma);
            data.u(index) = (1 - mu^2) * ((-(x0 - data.x(index)) / t) + a_l);
        elseif (x2 <= data.x(index) && data.x(index) <= x3)
            % Solution between x2 and x3 (between fan and contact)
            data.rho(index) = rho_star_L;
            data.P(index) = P_star;
            data.u(index) = v_star;
        elseif (x3 <= data.x(index) && data.x(index) <= x4)
            % Solution between x3 and x4 (between contact and shock)
            data.rho(index) = rho_star_R;
            data.P(index) = P_star;
            data.u(index) = v_star;
        elseif x4 < data.x(index)
            % Solution after x4 (after shock)
            data.rho(index) = rho_r;
            data.P(index) = P_r;
            data.u(index) = u_r;
        end
        data.e(index) = data.P(index) / ((gamma - 1) * data.rho(index));
    end
end

function y = FunctionTest3(P)
    global rho_l P_l u_l rho_r P_r u_r gamma mu a_l a_r
    
    % Test 3 Sod Problem
    A_R = 2 / ((gamma + 1) * rho_r);
    B_R = (mu^2 * P_r);
    A_L = 2 / ((gamma + 1) * rho_l);
    B_L = (mu^2 * P_l);
    
    y = (P - P_r)*((1-mu*mu)*(rho_r*(P + mu*mu*P_r))^-1)^(0.5)-(power(P_l , (gamma-1)/(2*gamma))-power(P , (gamma-1)/(2*gamma)))*(((1-mu*mu*mu*mu)*P_l^(1/gamma))*(mu*mu*mu*mu*rho_l)^-1)^(0.5);
end