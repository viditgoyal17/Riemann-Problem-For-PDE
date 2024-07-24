%-----THIS IS THE PROGRAM TO FIND THE EXACT SOLUTION FOR TEST-1------------
%------- OF A RIEMANN PROBLEM AND HELPS TO FIND THE DENSITY,---------------
%----------PRESSURE, VELOCITY AND INTERNAL ENERGY USING -------------------
%--------------------NEWTON RAPHSON METHOD---------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%---Test 2, called the 123 problem, has solution consisting of two strong rarefactions and
%---a trivial stationary contact discontinuity; the pressure p∗ is very small (close
%---to vacuum) and this can lead to difficulties in the iteration scheme to find p∗
%---numerically. On running the program the graph which appears shows solution profiles.
%---Test 2 is also useful in assessing the performance of numerical methods for low density flows function TEST2()
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


    global RHO_LEFT P_LEFT U_LEFT RHO_RIGHT P_RIGHT U_RIGHT GAMMA GAMMA_1 a_LEFT a_RIGHT
    
    RHO_LEFT = 1;
    P_LEFT = 0.4;
    U_LEFT = -2;

    RHO_RIGHT = 1;
    P_RIGHT = 0.4;
    U_RIGHT = 2;

    GAMMA = 1.4;
    GAMMA_1 = sqrt((GAMMA - 1) / (GAMMA + 1));
    a_LEFT = sqrt(GAMMA * P_LEFT / RHO_LEFT);
    a_RIGHT = sqrt(GAMMA * P_RIGHT / RHO_RIGHT);


    time = 0.15;
    data = Test_2(time);
    plotData(data);



function plotData(data)
    figure,
    subplot(2, 2, 1),
    plot(data.x, data.rho, '--m', 'LineWidth', 2);
    xlabel('x (m)');
    ylabel('Density (kg/m^3)');
    title('Density vs Position');
    grid on;

    subplot(2, 2, 2),
    plot(data.x, data.P, '-.c', 'LineWidth', 2);
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
    plot(data.x, data.e, '-g', 'LineWidth', 2);
    xlabel('x (m)');
    ylabel('Specific Internal Energy (J/kg)');
    title('Internal Energy vs Position');
    grid on;
end

function [data] = Test_2(t)

    global RHO_LEFT P_LEFT U_LEFT RHO_RIGHT P_RIGHT U_RIGHT GAMMA GAMMA_1 a_LEFT a_RIGHT
    
    x0 = 0;
    P_STAR = fzero(@FunctionTest2, 0.00189);
    V_STAR = U_RIGHT + 2 * (a_RIGHT / (GAMMA - 1)) * (-1 + power((P_STAR / P_RIGHT), ((GAMMA - 1) / (2 * GAMMA))));
    RHO_STAR_RIGHT = (RHO_RIGHT) * power((P_STAR / P_RIGHT), 1 / GAMMA);
    RHO_STAR_LEFT = (RHO_LEFT) * power((P_STAR / P_LEFT), 1 / GAMMA);

   
    x1 = x0 + (U_LEFT - a_LEFT) * t;
    x3 = x0 + V_STAR * t;
    x5 = x0 + (U_RIGHT + a_RIGHT) * t;
    
    a_star_L = power((GAMMA * P_STAR / RHO_STAR_LEFT), 0.5);
    a_star_R = power((GAMMA * P_STAR / RHO_STAR_RIGHT), 0.5);
    x2 = x0 + (V_STAR - a_star_L) * t;
    x4 = x0 + (V_STAR + a_star_R) * t;

    disp("p*")
    disp(P_STAR);
    disp("u*")
    disp(V_STAR);
    disp("rho*L")
    disp(RHO_STAR_LEFT);
    disp("rho*R")
    disp(RHO_STAR_RIGHT);


    N = 1000000;    
    x_min = -0.5;
    x_max = 0.5;

    x = linspace(x_min, x_max, N);
    data.x = x';
    data.rho = zeros(N, 1);   
    data.P = zeros(N, 1);     
    data.u = zeros(N, 1);   
    data.e = zeros(N, 1);    

    for index = 1:N
        if data.x(index) < x1
            data.rho(index) = RHO_LEFT;
            data.P(index) = P_LEFT;
            data.u(index) = U_LEFT;
        elseif (x1 <= data.x(index) && data.x(index) < x2)

            c = GAMMA_1 * GAMMA_1 * (U_LEFT - ((data.x(index)) / t)) + (1 - GAMMA_1 * GAMMA_1) * a_LEFT; 
            data.rho(index) = RHO_LEFT * power((c / a_LEFT), 2 / (GAMMA - 1));
            data.P(index) = P_LEFT * power((data.rho(index) / RHO_LEFT), GAMMA);
            data.u(index) = (1 - GAMMA_1 * GAMMA_1) * (((data.x(index)) / t) + a_LEFT + ((GAMMA - 1) / 2) * U_LEFT);
        elseif (x2 <= data.x(index) && data.x(index) < x3)

            data.rho(index) = RHO_STAR_LEFT;
            data.P(index) = P_STAR;
            data.u(index) = V_STAR;
        elseif (x3 <= data.x(index) && data.x(index) < x4)
            
            data.rho(index) = RHO_STAR_RIGHT;
            data.P(index) = P_STAR;
            data.u(index) = V_STAR;
        elseif (x4 <= data.x(index) && data.x(index) < x5)
            
            c = -GAMMA_1 * GAMMA_1 * (U_RIGHT - ((data.x(index)) / t)) + (1 - GAMMA_1 * GAMMA_1) * a_RIGHT; 
            data.rho(index) = RHO_RIGHT * power((c / a_RIGHT), 2 / (GAMMA - 1));
            data.P(index) = P_RIGHT * power((data.rho(index) / RHO_RIGHT), GAMMA);
            data.u(index) = (1 - GAMMA_1 * GAMMA_1) * (((data.x(index)) / t) - a_RIGHT + ((GAMMA - 1) / 2) * U_RIGHT);
        elseif x5 < data.x(index)

            data.rho(index) = RHO_RIGHT;
            data.P(index) = P_RIGHT;
            data.u(index) = U_RIGHT;
        end
        data.e(index) = data.P(index) / ((GAMMA - 1) * data.rho(index));
    end
end

function y = FunctionTest2(P)
    global RHO_LEFT P_LEFT U_LEFT RHO_RIGHT P_RIGHT U_RIGHT GAMMA GAMMA_1 a_LEFT a_RIGHT
    
    
    A_R = 2 / ((GAMMA + 1) * RHO_RIGHT);
    B_R = (GAMMA_1 * GAMMA_1 * P_RIGHT);
    A_L = 2 / ((GAMMA + 1) * RHO_LEFT);
    B_L = (GAMMA_1 * GAMMA_1 * P_LEFT);
    
    y = ((2*a_LEFT)/(GAMMA-1))*(-1+power((P/P_LEFT),((GAMMA-1)/(2*GAMMA)))) + ((2*a_RIGHT)/(GAMMA-1))*(-1+power((P/P_RIGHT),((GAMMA-1)/(2*GAMMA)))) +U_RIGHT-U_LEFT;
end