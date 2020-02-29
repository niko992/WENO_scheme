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

%error (compute with 3 different mesh size)
M=[100,200,300];
errorh=zeros(1,3);
errorm=zeros(1,3);

for i=1:3
    dx=2/M(i);
    x=(xLeft:dx:xRight)';
    N=length(x);
    U=[h0(x),m0(x)];
    m=2;
    T=2.0;
    U=ShallowWaterWENO(U,dx,CFL,m,T);

    errorh(1,i)=max(abs(U(:,1)-h0(x-T)));
    errorm(1,i)=max(abs(U(:,2)-m0(x-T))); 
end


%plot for two different mesh size
subplot(1,2,1)
loglog(1./M,errorh,'r');
hold on
loglog(1./M,(1./M).^(3),'--k')
legend('error h','O(Ne-3)');
title('error of h at final time')
subplot(1,2,2)
loglog(1./M,errorm,'r')
hold on
loglog(1./M,(1./M).^(3),'--k')
legend('error m','O(Ne-3)')
title('error of m at final time')
%print('33aerrLFm3','-dpdf')