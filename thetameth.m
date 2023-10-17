function [x,u] = thetameth(I,n,u0,f,bc,nu,theta)
% THETAMETH Theta-metodo.
% [U,X]=THETAMETH(I,N,U0,F,BC,NU,THETA) risolve l'equazione del calore
% utilizzando il THETA-metodo in tempo ed elementi finiti lineari in spazio
nx=n(1); h=(I(2)-I(1))/nx; x=(I(1):h:I(2))';
bc=bc(:);
nt=n(2); Deltat=(I(4)-I(3))/nt;
e=ones(nx+1,1);
K=spdiags([(h/(6*Deltat)-nu*theta/h)*e, (2*h/(3*Deltat)+2*nu*theta/h)*e, ...
    (h/(6*Deltat)-nu*theta/h)*e],-1:1,nx+1,nx+1);
B=spdiags([(h/(6*Deltat)+nu*(1-theta)/h)*e, (2*h/(3*Deltat)-nu*2*(1-theta)/h)*e, ...
    (h/(6*Deltat)+nu*(1-theta)/h)*e],-1:1,nx+1,nx+1);
M=h*spdiags([e/6, e*2/3,e/6],-1:1,nx+1,nx+1);
K(1,1)=1;       K(1,2)=0;     B(1,1:2)= 0;
K(nx+1,nx+1)=1; K(nx+1,nx)=0; B(nx+1,nx:nx+1)=0;
[L,U]=lu(K);
t=I(3);
uold=u0(x); fold=M*f(x,t);
for time=I(3)+Deltat:Deltat:I(4)
    fnew=M*f(x,time);
    b=theta*fnew+(1-theta)*fold+B*uold;
    b(1)=bc(1); b(end)=bc(2);
    y=L \ b;
    u=U \ y;
    uold=u; fold=fnew;
end
return
