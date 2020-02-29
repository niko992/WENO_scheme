function U = ShallowWaterWENO(U,h,CFL,m,T)
FinalTime=T;
xLeft=0; xRight=2;
x=xLeft:h:xRight;
time=0; tstep=0;

% Initialize reconstruction weights
Crec = zeros(m+1,m);
for r=-1:m-1
    Crec(r+2,:) = ReconstructWeights(m,r);
end

% Initialize linear weights
dw = LinearWeights(m,0);

% Compute smoothness indicator matrices
beta = zeros(m,m,m);
for r=0:m-1
    xl = -1/2 + [-r:1:m-r];
    beta(:,:,r+1) = betarcalc(xl,m);
end

while (time<FinalTime)
    maxvel=max(U(:,2)./U(:,1)+sqrt(U(:,1)));
    k=CFL*h/maxvel;
    if (time+k>FinalTime) 
         k = FinalTime-time; 
     end
 

  rhsU  = ShallowWaterWENOrhs(x, U,h,k,m,Crec,dw,beta,maxvel,time); 
  U1 = U + k*rhsU;
  rhsU  = ShallowWaterWENOrhs(x,U1,h,k,m,Crec,dw,beta,maxvel,time+k); 
  U2 = (3*U + U1 + k*rhsU)/4;
  rhsU  = ShallowWaterWENOrhs(x,U2,h,k,m,Crec,dw,beta,maxvel,time+k/2); 
  U  = (U + 2*U2 + 2*k*rhsU)/3;
  time = time+k; tstep = tstep+1;
end
end