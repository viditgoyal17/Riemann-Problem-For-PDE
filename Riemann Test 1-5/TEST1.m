%-----THIS IS THE PROGRAM TO FIND THE EXACT SOLUTION FOR TEST-1------------
%------- OF A RIEMANN PROBLEM AND HELPS TO FIND THE DENSITY,---------------
%----------PRESSURE, VELOCITY AND INTERNAL ENERGY USING -------------------
%--------------------NEWTON RAPHSON METHOD---------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%-----Test 1 is the so called Sod test problem. This is a very mild test 
%-----and its solution consists of a left rarefaction, a contact and a
%-----right shock. On running the program the graph which appears shows
%-----solution profiles for density, velocity, pressure and specific 
%-----internal energy across the complete wave structure, at time 0.25 units
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


function TEST1()

    % Define and globalizing initial conditions
    global rho_l P_l u_l rho_r P_r u_r gamma mu a_l a_r

    rho_l = 1;
    P_l = 1;
    u_l = 0;

    rho_r = 0.125;
    P_r = 0.1;
    u_r = 0;

    gamma = 1.4;
    mu = sqrt((gamma-1)/(gamma+1));
    a_l = sqrt(gamma*P_l/rho_l);
    a_r = sqrt(gamma*P_r/rho_r);

  
    time = 0.25;
    data = Test_1(time);
    plotData(data);

end

function plotData(data)

    figure,
    subplot(2,2,1),
    plot(data.x, data.rho, '--m', 'LineWidth', 2);
    xlabel('x (m)');
    ylabel('Density (kg/m^3)');
    title('Plot of Density vs Position');
    grid on;

    subplot(2,2,2),
    plot(data.x, data.P, '-.c', 'LineWidth', 2);
    xlabel('x (m)');
    ylabel('Pressure (Pa)');
    title('Plot of Pressure vs Position');
    grid on;

    subplot(2,2,3),
    plot(data.x, data.u, ':b', 'LineWidth', 2);
    xlabel('x (m)');
    ylabel('Velocity (m/s)');
    title('Plot of Velocity vs Position');
    grid on;

    subplot(2,2,4),
    plot(data.x, data.e, '-g', 'LineWidth', 2);
    xlabel('x (m)');
    ylabel('Specific Internal Energy (J/kg)');
    title('Plot of Internal Energy vs Position');
    grid on;
end

function [data] = Test_1(t)
   
    global rho_l P_l u_l rho_r P_r u_r gamma mu a_l a_r
    
    % Initial conditions
    x0 = 0;

    % Test 1
    P_star = fzero(@FunctionTest1, 0.31527);
    v_star = 2*(sqrt(gamma)/(gamma - 1))*(1 - power(P_star/P_l, (gamma - 1)/(2*gamma)));
    rho_star_R = rho_r*((P_star/P_r) + mu^2 )/(1 + mu*mu*(P_star/P_r));
    v_shock = v_star*((rho_star_R/rho_r)/( (rho_star_R/rho_r) - 1));
    rho_middle = rho_l*power((P_star/P_l), 1/gamma);


    x1 = x0 - (a_l)*t;  %Here u_l is zero so it was removed from the program to improve efficiency
    x3 = x0 + v_star*t;
    x4 = x0 + v_shock*t;


    c_2 = a_l - ((gamma - 1)/2)*v_star;
    x2 = x0 + (v_star - c_2)*t;

    disp("p*")
    disp(P_star);
    disp("u*")
    disp(v_star);
    disp("rho*L")
    disp(rho_middle);
    disp("rho*R")
    disp(rho_star_R);

    n_points = 100000;   
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
            data.rho(index) = rho_l;
            data.P(index) = P_l;
            data.u(index) = u_l;
        elseif (x1 <= data.x(index) && data.x(index) <= x2)
            
            c = mu*mu*((x0 - data.x(index))/t) + (1 - mu*mu)*a_l; 
            data.rho(index) = rho_l*power((c/a_l), 2/(gamma - 1));
            data.P(index) = P_l*power((data.rho(index)/rho_l), gamma);
            data.u(index) = (1 - mu*mu)*((-(x0-data.x(index))/t) + a_l);
        elseif (x2 <= data.x(index) && data.x(index) <= x3)
            
            data.rho(index) = rho_middle;
            data.P(index) = P_star;
            data.u(index) = v_star;
        elseif (x3 <= data.x(index) && data.x(index) <= x4)
            
            data.rho(index) = rho_star_R;
            data.P(index) = P_star;
            data.u(index) = v_star;
        elseif x4 < data.x(index)

            data.rho(index) = rho_r;
            data.P(index) = P_r;
            data.u(index) = u_r;
        end
        data.e(index) = data.P(index)/((gamma - 1)*data.rho(index));
    end
end

function y = FunctionTest1(P)
    global rho_l P_l u_l rho_r P_r u_r gamma mu a_l a_r
    
    
    A_R = 2/((gamma+1)*rho_r);
    B_R = (mu*mu*P_r);
    A_L = 2/((gamma+1)*rho_l);
    B_L = (mu*mu*P_l);
    
    y = (P - P_r)*sqrt((A_R/(B_R+P))) + ((2*a_l)/(gamma-1))*(power((P/P_l),(gamma-1)/(2*gamma))-1) + u_r - u_l;

end