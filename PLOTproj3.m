clc
clear
% DATA
xLeft=0; xRight=2;
data=1;
g=1;
CFL=0.5;
switch data
    case 1
        h0=@(x) 1+0.5*sin(pi*x);
        m0=@(x) 0.25*h0(x);
        S=@(x,t) [0.5*pi*(-0.75)*cos(pi*(x-t));0.5*pi*(-0.25+0.25^2+h0(x-t)).*cos(pi*(x-t))];
end

%plot the solution for two different mesh grid
dx=0.01;
x=(xLeft:dx:xRight)';
N=length(x);
U=[h0(x),m0(x)];
m=2;
T=2.0;
U=ShallowWaterWENO(U,dx,CFL,m,T);

%plot for two different mesh size
subplot(2,1,1)
plot(x,U(:,1),'r')
hold on
plot(x,h0(x-T),':b')
legend('h obtained with LF','true solution h');
title('solution with dx=0.01 and T=2')
subplot(2,1,2)
plot(x,U(:,2),'r')
hold on
plot(x,m0(x-T),':b')
legend('m obtained with LF','true solution m')
%print('33aLFm3','-dpdf')
