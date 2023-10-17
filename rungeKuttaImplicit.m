function [tk, uk] = rungeKuttaImplicit(f,tMax,y0,h,c,b,A)

% SUPER SAYAN GOD RUNGE KUTTA BETTER THAN MATLAB I HATE MY LIFE
% Runge-Kutta method of custom s, just assign the c,b vectors (line or
% column doesn't matter) and the matrix A.



[righe, colonne] = size(A);

if righe ~= colonne
    error('Matrix A must be square!');
end 
if length(c) ~= length(b)
    error('Vectors b and c must be the same length!')
end
if length(c) ~= righe
    error('Vectors b and c must have the same length as the rows/colums of A!')
end
if sum((triu(A)-diag(diag(A))),'all') ~= 0
    error('The given matrix A is not compatible, only tril matrix (explicit Runge-Kutta) are accepted!');
end

[ri, co] = size(b);
if ri > co
    b = b';
end


t0 = 0;
tk = t0:h:tMax;

n = length(tk);
uk = zeros(size(y0,1), n);

uk(:,1) = y0;
k = zeros(size(y0,1),righe);

nMax = 100;
toll = 1e-5;


for i=2:n
    for j=1:righe
        guess = f(tk(i-1) + c(j)*h, uk(:,i-1) + h * sum(A(j,1:((j-1)*(j>1)+1*(j==1))).*k(:,1:((j-1)*(j>1)+1*(j==1)) ) ));
        
        phi = @(kInc) f(tk(i-1) + c(j)*h, uk(:,i-1) + h * (A(j,1:((j-1)*(j>1)+1*(j==1))).*(k(:,1:((j-1)*(j>1)+1*(j==1))))*(j~=1) + kInc*A(j,j)));

        uPF = ptofis_sys(guess, phi, nMax, toll);
        k(:,j) = uPF(:,end);
    end

    uk(:,i) = uk(:,i-1) + h * sum(b.*k,2);
end
uk=uk';
tk=tk';