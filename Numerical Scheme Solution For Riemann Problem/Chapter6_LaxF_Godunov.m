% 1D Shock Tube Simulation using Lax-Friedrich and Adjustable Time Stepping
close all;
clear;
clc;

% Initialization of Parameters
N       = 100 ;             % Number of grid points
gamma   = 1.4 ;
endTime = 0.2 ;             % Physical End Time
CFL     = 0.5 ;

% Grid Initialization
x       = linspace(0,1,N+1) ;
xc      = 0.5 * ( x(1:N) + x(2:N+1) ) ;
xc(1)   = 0 ;               % Xmin
xc(N)   = 1 ;               % Xmax
time    = 0 ;

% Initial Conditions
densityRight    = 0.125 ;       densityLeft    = 1.00 ;
pressureRight   = 0.1  ;       pressureLeft   = 1 ;

rho     = zeros(N,1) ;
p       = zeros(N,1) ;
u       = zeros(N,1) ;

for i =  1:N
    if i<=N/2
        rho(i)  = densityLeft  ;
        p(i)    = pressureLeft ;        
    else
        rho(i)  = densityRight  ;
        p(i)    = pressureRight ;
    end
end
e   = p/(gamma-1) + 0.5*rho.*u.*u ;


new_rho = rho ;
new_u   = u   ;
new_e   = e   ;

while time <= endTime
    
    for i=2:N-1
        p(i)    = (gamma-1)*(e(i) - 0.5*rho(i)*u(i)*u(i)) ;
    end
    a       = sqrt(gamma*p./rho) ;
    lambda   = max(a) ;
    max_vel = max(u) ;
    
    dt      = CFL/N/(max_vel+lambda) ;  % Adjustable Time Step
    time    = time + dt ;
    
    for i=2:N-1
        dx          = xc(i) - xc(i-1) ;
        
        momRight = rho(i+1)*u(i+1) ;     
        rhoRight = rho(i+1) ;      
        uRight = u(i+1) ;      
        pRight = p(i+1) ;
        momPoint = rho(i)*u(i); 	
        rhoPoint = rho(i)   ;      
        uPoint = u(i)   ;      
        pPoint = p(i)   ;
        momLeft       = rho(i-1)*u(i-1) ;    	
        rhoLeft = rho(i-1) ;      
        uLeft = u(i-1) ;      
        pLeft = p(i-1) ;
        
        velFluxRight  = rhoRight*uRight*uRight + pRight ;    
        eRight = e(i+1) ;
        velFluxPoint  = rhoPoint*uPoint*uPoint + pPoint ;    
        ePoint = e(i)   ;
        velFluxLeft  = rhoLeft*uLeft*uLeft + pLeft ;    
        eLeft = e(i-1) ;
        
        energyFluxRight   = uRight * ( eRight + pRight ) ;
        energyFluxPoint   = uPoint * ( ePoint + pPoint ) ;
        energyFluxLeft   = uLeft * ( eLeft + pLeft ) ;
        
        rhoFluxRight   = 0.5*( momPoint + momRight ) -0.5*lambda*( rhoRight - rhoPoint ) ;
        rhoFluxLeft   = 0.5*( momPoint + momLeft ) -0.5*lambda*( rhoPoint - rhoLeft ) ;
        
        velFluxRight   = 0.5*( velFluxPoint + velFluxRight ) -0.5*lambda*( momRight - momPoint ) ;
        velFluxLeft   = 0.5*( velFluxPoint + velFluxLeft ) -0.5*lambda*( momPoint - momLeft ) ;
        
        energyFluxRight    = 0.5*( energyFluxPoint + energyFluxRight ) -0.5*lambda*( eRight - ePoint );
        energyFluxLeft    = 0.5*( energyFluxPoint + energyFluxLeft ) -0.5*lambda*( ePoint - eLeft ) ;
        
        new_rho(i)  = rhoPoint - dt/dx * ( rhoFluxRight - rhoFluxLeft ) ;
        velFlux    = momPoint - dt/dx * ( velFluxRight - velFluxLeft ) ;
        
        new_u(i)    = velFlux/new_rho(i) ;
        new_e(i)    = ePoint - dt/dx * ( energyFluxRight - energyFluxLeft ) ;
        
    end
    
    rho     = new_rho ;
    u       = new_u ;
    e       = new_e ;
    
end



% Gudunov's Method
% Gudunov Scheme
L = 1;
dxGud = 0.001;
xGud = (0:dxGud:L);
NGud = (L/dxGud)+1;
TGud = 0.2;
dtGud = 0.0001;
nGud = (TGud/dtGud);
C = dtGud/dxGud;
gamma = 1.4;

% Initial conditions
for i = 1:NGud
    if xGud(i) <= 0.5
        rhoGud(i) = 1;
        uGud(i) = 0;
        pGud(i) = 2.5*(gamma - 1);
    else
        rhoGud(i) = 0.125;
        uGud(i) = 0;
        pGud(i) = 0.25*(gamma - 1);
    end
end
ET(:) = ((1/2).*rhoGud(:).*(uGud(:).^2)) + (pGud(:)/(gamma - 1));



for i = 1:nGud

    % 1) Construct Q and E vectors

    QGud(1,:) = rhoGud(:);
    QGud(2,:) = rhoGud(:).*uGud(:);
    QGud(3,:) = ET(:);

    EGud(1,:) = rhoGud(:).*uGud(:);
    EGud(2,:) = EGud(1,:).*uGud(1,:) + pGud(1,:);
    EGud(3,:) = ET(:).*uGud(:) + pGud(:).*uGud(:);

    % 2) Calculate Q(n+1) using Gudunov method
    Q1 = QGud;
    flux = zeros(3,NGud);
    for k = 1:NGud-1
        c = sqrt(gamma*pGud(k)/rhoGud(k));
        A = [abs(QGud(2,k)) abs(QGud(2,k)+c) abs(QGud(2,k)-c)];
        al = max(A);
        flux(:,k) = 0.5*(EGud(:,k) + EGud(:,k+1)) - 0.5*(al)*(Q1(:,k+1) - Q1(:,k));
    end
    for k = 2:NGud-1
        QGud(:,k) = Q1(:,k) - C*(flux(:,k) - flux(:,k-1));
    end
    
    % 3) Update solution (rho,u,p)
    rhoGud(:) = QGud(1,:);
    uGud(1,:) = QGud(2,:)./rhoGud(1,:);
    ET(:) = QGud(3,:);
    pGud(1,:) = (ET(1,:) - ((1/2).*rhoGud(1,:).*(uGud(1,:).^2)))*(gamma - 1);
    
end






time = 0.2;
data = Test_1(time);


figure(1)
hold on
%plot(density(:,1),density(:,2),'r--','MarkerSize',12,'MarkerIndices',1:10:length(velocity),'LineWidth',2);
plot(data.x+0.5,data.rho,'r-','LineWidth',2);
plot(xc,rho,'b-','MarkerSize',12,'MarkerIndices',1:10:141,'LineWidth',1.5);
plot(xGud, rhoGud, 'Color', [0.2, 0.8, 0.2], 'LineWidth', 1.5); % Green color
xlabel(' x ','FontSize',12,'FontWeight','bold');
ylabel(' Density ','FontSize',12,'FontWeight','bold');
legend('Exact','Lax Friedrich','Gudunov','Location','northeast','FontSize',11);
%print(gcf,'Density.jpg','-dpng','-r300');

figure(2)
hold on
%plot(pressure(:,1),pressure(:,2),'r--','MarkerSize',12,'MarkerIndices',1:10:length(velocity),'LineWidth',2);
plot(data.x+0.5,data.P,'r','LineWidth',2);
plot(xc,p,'b','MarkerSize',12,'MarkerIndices',1:10:141,'LineWidth',1.5);
plot(xGud, pGud, 'Color', [0.2, 0.8, 0.2], 'LineWidth', 1.5); % Green color
xlabel(' x ','FontSize',12,'FontWeight','bold');
ylabel(' Pressure ','FontSize',12,'FontWeight','bold');
legend('Exact','Lax Friedrich','Gudunov','Location','northeast','FontSize',11);
%print(gcf,'Pressure.jpg','-dpng','-r300');

figure(3)
hold on
%plot(velocity(:,1),velocity(:,2),'r--','MarkerSize',12,'MarkerIndices',1:10:length(velocity),'LineWidth',2);
plot(data.x+0.5,data.u,'r','LineWidth',2);
plot(xc,u,'b','MarkerSize',12,'MarkerIndices',1:10:141,'LineWidth',1.5);
plot(xGud, uGud, 'Color', [0.2, 0.8, 0.2], 'LineWidth', 1.5); % Green color
xlabel(' x ','FontSize',12,'FontWeight','bold');
ylabel(' Velocity ','FontSize',12,'FontWeight','bold');
legend('Exact','Lax Friedrich','Gudunov','Location','south','FontSize',11);
%print(gcf,'Velocity.jpg','-dpng','-r300');
