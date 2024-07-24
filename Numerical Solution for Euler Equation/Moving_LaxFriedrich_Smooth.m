clear
clc

%Define
xmin = -1;
xmax = 2;
N=100;
dt= 0.016;
t=0;
tmax=1;
v=1;

%Dicretise
dx=(xmax-xmin)/N;
x= xmin-dx : dx : xmax+dx;
cou=(v*dt)/dx;
l=cou;
l1=1-cou;
l2=cou+1;
l3=2-cou;

%Initial Conditions
u0= exp(-8*(x).^2);
u=u0;
unpl=u0;

%loop through time
nsteps = tmax/dt;

for n = 1:nsteps

    % Calculate boundary conditions
    u(1)=u(3);
    u(N+3)=u(N+1);

    % Calculate the FOU scheme
    for i=2 : N+2
        %Lax-Friedrich
        unpl(i) = (0.5*l2*u(i-1))+(0.5*l1*u(i+1));

        %Lax-Wendroff
        %unpl(i) = (0.5*l2*l*u(i-1))+(l2*l1*u(i))+(0.5*l1*l*u(i+1));

        %Warming-beam
        %unpl(i) = (-0.5*l1*l*u(i-2))+(l*l3*u(i-1))-(0.5*l1*l3*u(i));

    end

    % update t and u
    t = t+dt;
    u = unpl;

    % calculate exact solution 
    exact = exp(-8*(x-v*t).^2);

    % plot solution
    plot(x,exact,'r-');
    hold on
    plot(x,u,'bo','markerfacecolor','b');
    hold off
    axis([xmin xmax -0.5 1.5])
    xlabel('x','fontsize',16)
    ylabel('U(t,x)','FontSize',16);
    title(sprintf('time = %1.3f',t),'FontSize',16);
    shg
    pause(dt);
end




