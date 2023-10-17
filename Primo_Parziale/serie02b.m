clc
clear
n = 100;
R1 = 1;
R2 = 2;
A = diag([1 -R2*ones(1,n-1)]) + diag(R1*ones(1,n-1),-1);
A(1,2:end) = 1;
nnz(A);
As = sparse(A);
% [L_A,U_A,P,Q] = lu(As);
% spy(U_A)
Dinv = diag(1./diag(A));
B_j = eye(n)-Dinv*A;
T = tril(A);
B_gs = eye(n)-inv(T)*A;
rho_j = max(abs(eig(B_j)))
rho_gs = max(abs(eig(B_gs)));

% 5
b = [2 ones(1,n-1)]';
x0 = zeros(n,1);
toll = 1e-6;
nmax = 1000;
[x k] = jacobi(A,b,x0,toll,nmax);
%% 1.2
clear
clc
n = 7;
A = diag(9*ones(1,n)) + diag(-3*ones(1,n-1),1) + diag(-3*ones(1,n-1),-1) + diag(ones(1,n-2),-2) + diag(ones(1,n-2),2);
b = [7 4 5 5 5 4 7]';
for i = 1:n
    r = abs(A(i,:));
    if abs(A(i,i))<max(r)
        error('La matrice non è a dominanza stretta')
    end
end
fprintf('\nLa matrice è a dominanza stretta\n');
if A ~= A'
    error('La matrice non è simmetrica');
elseif eig(A)<0
    error('La matrice non è def positiva')
end
fprintf('\nLa matrice è simmetrica e definita positiva\n');
x0 = zeros(n,1);
toll = 1e-6;
nmax = 1000;
[x,k] = gauss_seidel(A,b,x0,toll,nmax);

[x,k] = jacobi(A,b,x0,toll,nmax);


%% 2.1
clear
clc
n = 50;
A = diag(4*ones(1,n)) + diag(-1*ones(1,n-1),1) + diag(-1*ones(1,n-1),-1) + diag(-1*ones(1,n-2),-2) + diag(-1*ones(1,n-2),2);
b = 0.2*ones(n,1);
x0 = zeros(n,1);
if A ~= A'
    error('La matrice non è simmetrica')
end
if eig(A)<0
    fprintf('\nLa matrice non è definita positiva\n');
end
condiz = max(abs(eig(A)))/min(abs(eig(A)))

P = eye(n);
toll = 1e-5;
nmax = 1e4;

alpha = 0.2;
[x, k1] = richardson(A,b,P,x0,toll,nmax,alpha);
B_it = eye(n)-alpha*A;
rho_1 = max(abs(eig(B_it)));

alpha = 0.3;
B_it = eye(n)-alpha*A;
rho_2 = max(abs(eig(B_it)));
[x, k2] = richardson(A,b,P,x0,toll,nmax,alpha);

alpha_opt = 2/(min(eig(inv(P)*A))+max(eig(inv(P)*A)));
B_it = eye(n)-alpha_opt*A;
rho_opt = max(abs(eig(B_it)));
[x, k_opt] = richardson(A,b,P,x0,toll,nmax,alpha_opt);
fprintf('k1 = %d  k2 = %d  kopt = %d',k1, k2, k_opt);
fprintf('\nrho_1 = %f  rho_2 = %f  rho_opt = %f', rho_1, rho_2, rho_opt);


P = tril(A);
alpha = 1;
B_alpha     = eye(n) - alpha * inv(P) * A;
fprintf('Raggio spettrale: %f\n', max(abs(eig(B_alpha))))
[x1, k1] = richardson(A,b,P,x0,toll,nmax,alpha);
k1
[x2,k2] = gauss_seidel(A,b,x0,toll,nmax);
k2

fprintf('\nscarto soluzioni: %e\n',max(abs(x1-x2)));

%% 2.6
clear
clc
n = 50;
A = diag(4*ones(1,n)) + diag(-1*ones(1,n-1),1) + diag(-1*ones(1,n-1),-1) + diag(-1*ones(1,n-2),-2) + diag(-1*ones(1,n-2),2);
b = 0.2*ones(n,1);
x0 = zeros(n,1);
toll = 1e-5;
nmax = 1e4;
P = diag(2*ones(1,n)) + diag(-1*ones(1,n-1),-1) + diag(-1*ones(1,n-1),1);
if P ~= P'
    error('Matrice P non simmetrica');
end
if eig(P)<0
    error('La matrice non è def. positiva');
end
v = eig(inv(P)*A);
alpha = 2/(max(v)+min(v));
B = eye(n)-alpha*inv(P)*A;
fprintf('\nRaggio spettrale della matrice di iterazione è: %f',max(abs(eig(B))))
fprintf('\nIl numero di condizionamento è: %f\n',cond(inv(P)*A))
[x k] = richardson(A,b,P,x0,toll,nmax,alpha);
k
[x k] = richardson(A,b,P,x0,toll,nmax);
k

%% GRADIENTE
clear
clc
A1 = [6.8 -2.4; -2.4 8.2];
b1 = [4 8]';
[V,D]=eig(A1);
phi = @(x,y,A,b) 0.5*A(1,1)*x.^2 + 0.5*(A(1,2) + A(2,1))*x.*y + 0.5*A(2,2)*y.^2 - b(1).*x-b(2).*y;
x = linspace(-10,10,100);
[X,Y]=meshgrid(x,x);

% surf(X,Y,phi(X,Y,A1,b1));
% contour(X,Y,phi(X,Y,A1,b1),20);
 A2 = V*diag([1 D(2,2)])*V';
% surf(X,Y,phi(X,Y,A2,b1));
% contour(X,Y,phi(X,Y,A2,b1));


x0 = [-9 -9]';
[ric_alpha1, x1, k1] = richardson_it(A2,b1,eye(2),x0,1e-7,1000,0.05);
[ric_alpha2, x2, k2] = richardson_it(A2,b1,eye(2),x0,1e-7,1000,0.24);
[grad, x3, k3] = richardson_it(A2,b1,eye(2),x0,1e-7,1000);

subplot(2,2,1)
contour(X,Y,phi(X,Y,A2,b1));
title(['Alpha=0.05, ' num2str(k1) ' iterazioni']);
axis('equal');
hold on;
plot(ric_alpha1(1,:),ric_alpha1(2,:),'-or','LineWidth',2);

subplot(2,2,2)
contour(X,Y,phi(X,Y,A2,b1));
title(['Alpha=0.2',  num2str(k2) ' iterazioni']);
axis('equal');
hold on;
plot(ric_alpha2(1,1:6),ric_alpha2(2,1:6),'-or','LineWidth',2);

subplot(2,2,3)
contour(X,Y,phi(X,Y,A2,b1));
title(['Gradiente' ,  num2str(k1) ' iterazioni']);
axis('equal');
hold on;
plot(grad(1,:),grad(2,:),'-or','LineWidth',2);

%% punto 3
A2 = V*diag([1 D(2,2)])*V';
b1 = [4 8]';
x0 = [-9 -9]';
phi = @(x,y,A,b) 0.5*A(1,1)*x.^2 + 0.5*(A(1,2) + A(2,1))*x.*y + 0.5*A(2,2)*y.^2 - b(1).*x-b(2).*y;
P = [1.0912 -0.8587; -0.8587 1.5921];
surf(X,Y,phi(X,Y,inv(P)*A2,inv(P)*b1));
contour(X,Y,phi(X,Y,inv(P)*A2,inv(P)*b1),20);

[x_rich3, x, k_rich3] = richardson_it(inv(P)*A2,inv(P)*b1,eye(2),x0,1e-7,1000);
[x_rich_prec, x, k_rich_prec] = richardson_it(A2,b1,P,x0,1e-7,1000);
figure
contour(X,Y,phi(X,Y,A2,b1));
title(['Gradiente']);
axis('equal');
hold on;
plot(x_rich3(1,:),x_rich3(2,:),'-or','LineWidth',2);
figure;
contour(X,Y,phi(X,Y,A2,b1));
title(['Gradiente']);
axis('equal');
hold on;
plot(x_rich_prec(1,:),x_rich_prec(2,:),'-or','LineWidth',2);

%% punto 4
A2 = V*diag([1 D(2,2)])*V';
b1 = [4 8]';
x0 = [-9 -9]';
phi = @(x,y,A,b) 0.5*A(1,1)*x.^2 + 0.5*(A(1,2) + A(2,1))*x.*y + 0.5*A(2,2)*y.^2 - b(1).*x-b(2).*y;
[x_grad_con,k] = conjgrad(A2,b1,x0,1e-7,1000);
x = linspace(-10,10,100);
[X,Y]=meshgrid(x,x);
surf(X,Y,phi(X,Y,A2,b1));
contour(X,Y,phi(X,Y,A2,b1),20);
hold on;
plot(x_grad_con(1,:),x_grad_con(2,:),'-or','LineWidth',2);

%% Es 3.2
toll = 1e-6;
nmax = 5000;
h = 1;
N = 5:20;
for n = 5:20
    A = diag(4*ones(1,n))+diag(ones(1,n-1),-1)+diag(2*ones(1,n-2),-2)+diag(ones(1,n-1),1)+diag(2*ones(1,n-2),2);
    P = tril(A);
    b = ones(n,1);
    x0 = zeros(n,1);
    
    [x, k] = richardson(A,b,P,x0,toll,nmax);
    it(h) = k;
    con(h)=cond(A);
    h = h+1;
    
end
subplot(2,1,1)
plot(N,it);
title('ITERAZIONI');
subplot(2,1,2)
plot(N,con);
title('CONDIZIONAMENTO');

% grad coniugato ci sono le radici in d
% grad precondizionato e non, non ci sono le radici
%%
A = [-5 0.01 1; 0 -1 -11; 0 0 1];
eig(A)
%%
n = 300;
A = diag(6*ones(n,1))+diag(-2*ones(n-1,1),-1)+diag(-2*ones(n-1,1),1)+diag(ones(n-2,1),-2)+diag(ones(n-2,1),2);
x = ones(n,1);
b = A * x;
Pj = diag(diag(A));
Pgs = tril(A);
Bj = eye(n) - Pj \ A;
Bgs = eye(n) - (Pgs) \ A;
rho_j = max(abs(eig(Bj)));
rho_gs = max(abs(eig(Bgs)));
%G_S converge più rapidamente
x0 = b;
tol = 1e-6;
nmax = 1000;
[sol k] = richardson(A,b,Pgs,x0,tol,nmax,1);
disp('#####Iterazioni:');
k
disp('#####Errore relativo');
err_rel = norm(x-sol)/norm(x)
disp('#####Residuo normalizzato');
res_norm = norm(b-A*sol)/norm(b)
disp('#####Errore stimato');
err_stima = cond(A)*res_norm

k = 20;
d = (cond(A)-1)/(cond(A)+1);
disp('Errore dopo 20 iterazioni');
err_k = (d^k)*norm(x-x0)

beta = [2 3 4 5];
for i = 1 : length(beta)
    P = diag(beta(i)*ones(n,1)) + diag(-ones(n-1,1),-1) + diag(-ones(n-1,1),1);
    d(i) = (cond(P\A)-1)/(cond(P\A)+1);
end
[d_migliore i] = min(d);
beta(i)
d_migliore
[sol, f, r, k] = pcg(A,b,1e-6,1000,[],[],b);
k
err_norm = sqrt((sol-x)'*A*(sol-x))

