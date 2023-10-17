function [u_h,t_h] = Runge_Kutta_implicito(M,A,b,c,tf,y0,h,g)
%RUNGE-KUTTA IMPLICITO con A triangolare inferiore e 
% b e c vettori colonna dall'array di Butcher, 
% M e g t.c. y'(t) = My(t) + g(t) (M=matrice iniziale A mxm)

m=length(y0);
s=length(b);

if nargin==6
    g = @(t) zeros(m,1).*t;
end


u_h=y0;
t_h=0;

K=zeros(m,s);

for t=0:h:tf-h
    for j=1:s
        K(:,j)=(eye(m) - h*A(j,j)*M)\(h*M*K(:,[1:j-1 j+1:end])*(A(j,[1:j-1 j+1:end]))' + M*u_h(:,end) + g(t + c(j)*h));
    end
    u_h = [u_h, u_h(:,end) + h*sum(b.*K,2)];
end
    

end