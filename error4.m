clc
clear
% DATA
xLeft=0; xRight=2;
data=2;
g=1;
CFL=0.5;
switch data
    case 1
        h0=@(x) 1+0.5*sin(pi*x);
        m0=@(x) 0.25*h0(x);
        S=@(x,t) [0.5*pi*(-0.75)*cos(pi*(x-t));0.5*pi*(-0.25+0.25^2+h0(x-t)).*cos(pi*(x-t))];
        %S=@(x,t) [(pi/2)*(u-1)*cos(pi*(x-t));(pi/2)*cos(pi*(x-t)).*(-u+u^2+g*hIC(x-t))];
    case 2
        h0=@(x) 1-0.1*sin(pi*x);
        m0=@(x) 0+0.*x;
        S=@(x,t) [0;0];
    case 3
        h0=@(x) 1-0.2*sin(2*pi*x);
        m0=@(x) 0.5+0.*x;
        S=@(x,t) [0.*x;0.*x];
    case 4
        h0=@(x) 1+0.*x;
        m0=@(x) -1.5.*(x<=1)+0.*(x>1);
        S=@(x,t) [0.*x;0.*x];
end


% compute reference solution 
dx=0.005;
x=(xLeft:dx:xRight)';
N=length(x);
U=[h0(x),m0(x)];
m=2;
T=2.0;
Uex=ShallowWaterWENO(U,dx,CFL,m,T);
xex=x;

%plot the solution for two different mesh grid
M=[100,200,300];
errorh=zeros(1,3);
errorm=zeros(1,3);
n=zeros(1,3);
for i=1:3
    dx=2/M(i);
    x=(xLeft:dx:xRight)';
    N=length(x);
    U=[h0(x),m0(x)];
    m=2;
    T=2.0;
    U=ShallowWaterWENO(U,dx,CFL,m,T);
    errorh(1,i)=max(abs(U(:,1)-interp1(xex,Uex(:,1),x,'linear')));
    errorm(1,i)=max(abs(U(:,2)-interp1(xex,Uex(:,2),x,'linear'))); 
    n(1,i)=N;
end

 %p1 = polyfit(log(n),log(errorh),1);
 %p2 = polyfit(log(n),log(errorm),1);
%plot for two different mesh size
subplot(1,2,1)
loglog(1./M,errorh,'r');
hold on
loglog(1./M,(1./M).^(3),'--k')
legend('error h','O(Ne-3)');
title('error of h at final time')
%title(['error h','. Slope = ',num2str(p1(1))],'FontSize',20);
subplot(1,2,2)
loglog(1./M,errorm,'r')
hold on
loglog(1./M,(1./M).^(3),'--k')
legend('error m','O(Ne-3)')
title('error m at final time')
%title(['error of m','. Slope = ',num2str(p2(1))],'FontSize',20);
%print('34a3errRoe','-dpdf')


