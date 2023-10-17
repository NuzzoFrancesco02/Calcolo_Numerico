%% PRIMA ITINERE
%% 1
sol = x_float(2,5,[1 0 1 0 1],3,0)
%5.25
%%
x = 2;
cont = 0;
for i = 1 : 50
    
    x = (x+5/x)/2;
    cont = cont + 3;
end
cont
%% 2
n = 90;
3*(n-1) + 10*(2*(n-1) + 3*(n-1) + 1) %fatt e 10 sost in avanti e indietro
%% 3
syms g;
Q = [1/sqrt(2) 0; 0 1; 1/sqrt(2) 0];
R = [2 -g; 0 1];
b = [1 2 2]';
x_star = R\(Q'*b);
(simplify(x_star))
%% 5
A = diag([4 4 3 6 -1])+ diag([1 -1 1 9],-1);
A(1,2)=-1;
compass(eig(A))
shift_finder(A,3)
invpowershift(A,3.5,1e-8,1000,ones(5,1))
%%
ceil(log2((5)/(1e-3))-1)
%%
f = @(x) ((x.^4)./4)+((x.^2)./2)-3.*x+5;
df = Jac(f);
ddf = Jac(df);
sol = newton(0,5,1e-10,df,ddf)
sol(2)
sol(3)
sol(6)
int = [0:0.01:3];
plot(int,df(int))
grid on
%%
f = @(x) sin(pi.*x./3).*(x-3).^2;
df = Jac(f);
ddf = Jac(df);
dddf = Jac(ddf);
sol = newton(4,1,1e-10,f,df,3)
%%
phi = @(x) x -(1-exp(1-3.*x))./3;
dphi = Jac(phi);
ddphi = Jac(dphi);
%%
g = -1/6;
phi = @(x) g*(x^2 - 4*x -5) + 5;
ptofis(4.9,phi,100,1e-10)
%% ESERCIZIO
n = 225;
A = full(gallery('poisson', 15));
b = 2*ones(n,1);
% verifico che A'= A e che l'autovalore reale più piccolo sia positivo:
if isequal(A,A') && min(eig(A))>0
    disp('A è simmetrica e definita positiva');
else
    disp('A non è simmetrica o definita positiva');
end
%{ Il metodo computazionalmente più conveniente è la fattorizzazione di
%Cholesky, in quanto A si è verificato che è simmetrica e definita
%positiva e questo metodo diretto ha un costo computazionale di O((n^3)/3),
%la metà della fattorizzazione LU con un costo di O((2n^3)/3). Il costo
%dell'algoritmo è quindi di (n^3)/3 per la fattorizzazione, di n^2 per
%risolvere il sistema R'*y=b, e n^2 operazioni per risolvere il sistema
%triangolare R*x=y. Complessivamente si ha (n^3)/3 + 2*n^2, cioè O(n^3)/3)
%}

%%
R = chol(A); % trovo la matrice R tale che R'*R = A
y = fwsub(R',b); % risolvo il sist triangolare inferiore
sol = bksub(R,y);% risolvo il sist triangolare superiore
sol(10)
matl = A\b;
matl(10)
err_stim_chol = cond(A)*norm(b-A*sol)/norm(b) %9.0107e-13
%%
% calcolo il raggio spettrale della matrice di iterazione
P = diag(diag(A));
B = eye(n)-P\A;
rho_jac = max(abs(eig(B))) % 0.9808
% vedo il metodo di Jacobi come un caso generale di richardson stazionario
% con alpha ottimale = 1, infatti:
a_opt = 2/(max(eig(P\A))+min(eig(P\A))); 
K_jac = max(eig(P\A))/min(eig(P\A));
d_jac = (K_jac-1)/(K_jac+1);
k_jac = 100;
abb = d_jac^k_jac   % 0.1437
%%
phi = @(y) 0.5*y'*A*y - y'*b;
x0 = b;
x1 = richardson(A,b,eye(n),x0,1e-10,1);
phi(x0) % -780
phi(x1) % -1.6603e+03
%%
P1 = diagonals([-1 4 -1],n);
R2 = ichol(sparse(A));
P2 = R2'*R2;

K1 = max(eig(P1\A))/min(eig(P1\A)) % 58.2413
K2 = max(eig(P2\A))/min(eig(P2\A)) % 10.0966
% P2 ha un condizionamento spettrale più basso, 
% converge più velocemente
c = (sqrt(K2)-1)/(sqrt(2)+1);
k = 20;
% |xk-x|_A <= (2*c^k)/(1+c^(2*k))*|x0-x|_A
% abb = (2*c^k)/(1+c^(2*k))
abb = (2*c^k)/(1+c^(2*k)) % 0.2499

% DEFINIZIONE FUNZIONE
%{
function M = diagonals(v,n)
    u = length(v);
    M = zeros(n);
    if mod(u,2)~=0
        iter = ceil(u/2);
        ind1 = iter;
        ind2 = iter;
        for i = 0 : iter-1
            if i == 0
                M = M + diag(v(ind1)*ones(n-i,1),i);
            else
                M = M + diag(v(ind1)*ones(n-i,1),i);
                M = M + diag(v(ind2)*ones(n-i,1),-i);
            end
            ind1 = ind1 + 1;
            ind2 = ind2 -1;
        end
    else
        iter = ceil(u/2);
        v = [v(1:iter) 0 v(iter+1:end)];
        iter = iter + 1;
        ind1 = iter;
        ind2 = iter;
        for i = 0 : iter-1
            M = M + diag(v(ind1)*ones(n-i,1),i);
            M = M + diag(v(ind2)*ones(n-i,1),-i);
            ind1 = ind1 + 1;
            ind2 = ind2 -1;
        end
    end
end
%}

%%
x0 = ones(n,1);
l1 = max(eig(A));
[lambda,x,iter,lambda_vect]=eigpower(A,1e-10,1e5,x0);
err = abs(l1-lambda_vect);
p = stimap(err);
ord = ceil(p(end)) % 1
fatt = err(end)/(err(end-1)^ord) % 0.9260
%%
l_1 = qrbasic(A,1e-8,1);
l_1(2)
l_2 = qrbasic(A,1e-8,2);
l_2(2)
l_20 = qrbasic(A,1e-8,20);
l_20(2)

l = sort(eig(A),'descend');
l(2);

