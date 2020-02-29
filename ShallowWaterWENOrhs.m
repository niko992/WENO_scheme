function [dq] = ShallowWaterWENOrhs(x,u,h,k,m,Crec,dw,beta,maxvel,time)
N = length(x); dq = zeros(N,2);
ql = zeros(N,2); qr = zeros(N,2); qm = zeros(N,2); qp = zeros(N,2);

% Extend data and assign boundary conditions Fort periodic boundary
[xe,He] = extend(x,u(:,1),h,m,'P',0,'P',0);
[xe,Me] = extend(x,u(:,2),h,m,'P',0,'P',0);
% Extend data and assign boundary conditions Fort Open boundary
%[xe,He] = extend(x,u(:,1),h,m,'O',0,'O',0);
%[xe,Me] = extend(x,u(:,2),h,m,'O',0,'O',0);

% needed only for the first case---
g=1;
h0=@(x) 1+0.5*sin(pi*x);
m0=@(x) 0.25*h0(x);
S=@(x,t) [0.5*pi*(-0.75)*cos(pi*(x-t));0.5*pi*(-0.25+0.25^2+h0(x-t)).*cos(pi*(x-t))];

%define cell left and right interface values
Hl = zeros(N+2,1); Hr = zeros(N+2,1); Ml = zeros(N+2,1); Mr = zeros(N+2,1);

for i=1:N+2
 [Hl(i),Hr(i)] = WENO(xe(i:(i+2*(m-1))),He(i:(i+2*(m-1))),m,Crec,dw,beta);
 [Ml(i),Mr(i)] = WENO(xe(i:(i+2*(m-1))),Me(i:(i+2*(m-1))),m,Crec,dw,beta);
end

% Compute residual
qr = [Hr(2:N+1) Mr(2:N+1)]; ql = [Hl(2:N+1) Ml(2:N+1)];
qm = [Hr(1:N) Mr(1:N)]; qp = [Hl(3:N+2) Ml(3:N+2)];


% Lax Friedrich flux with source non zero
   dq = - (LaxFriedFLUX(qr,qp,k,g,h) - ...
               LaxFriedFLUX(qm,ql,k,g,h))/h+ S(x,time)';
% Lax Friedrich flux with source  zero
% dq = - (LaxFriedFLUX(qr,qp,k,g,h) - ...
%             LaxFriedFLUX(qm,ql,k,g,h))/h;

% Roe flux with source non zero
% dq = - (RoeFLUX(qr,qp,g) -...
%                RoeFLUX(qm,ql,g))/h + S(x,time)';

% Roe flux with source zero
 % dq = - (RoeFLUX(qr,qp,g) -...
 %               RoeFLUX(qm,ql,g))/h;
end