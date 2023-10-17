%% TDE
%% 1
0.5*eps_m(2,7,'eps')
%% 2
clear
clc
A = 216;
x = A;
for i = 1 : 10
    x = (A/(3*x^2))+2*(x)/3;
end
x
%% 3
A = [ 4 -1 -1; -1 4 -2; -1 -2 6];
b = [1 1 1]';
P = diag([4 4 6]);
a_opt = 2/(max(eig(P\A))+min(eig(P\A)))
richardson(A,b,P,zeros(3,1),1e-10,5,a_opt)
%% 4
syms t;
A = [5 -1 0; -1 t 1; 0 0 3];
P = diag(diag(A));
B = eye(3)-P\A;
pretty(simplify(abs(eig(B))))
% |teta| > 1/5
%% 5
A = hilb(7);
s = 0.2;
x0 = ones(7,1);
for i = 0 : 2
    invpowershift(A,s,1e-10,i,x0)
end
%% 6
syms t
A = [3 t 1; -t 1 4; 0 0 9];
pretty(simplify(eig(A)))
% -1 < t < 1
%% 7
f = @(x) cos(pi.*x).^2;
df = Jac(f);
df(0.5)
% df = 0, p = 1
%% 8
f = @(x) (x-1).*log(x);
df = Jac(f);
ddf = Jac(df);
ddf(1); % ddf(1) ≠ 0 --> m = 2
x = newton(0.9,1e-10,1,f,df,2);
x(2)
%% 9 
% p = 2
mu = (1e-2)/((1e-1)^2);
e_k2 = mu*((1e-2)^2)
%% 10
% alpha1: mai
% alpha2: 1<teta<3
%% 1.1
%{ Il metodo di Thomas: per risolvere un sistema A*x = b in modo diretto con una
% matrice tridiagonale con il metodo di Thomas ha un costo di 8*n-7
% operazioni, molto più conveniente rispetto a LU (O((2*n^3)/3)), questo è
% dovuto alla presenza di molti elementi nulli. Thomas fattorizza la matrice
% A in due matrici bidiagonali (quella inferiore con elementi 1 sulla diagonale), 
% in questo modo le sostituzioni in avanti e all'indietro costano molto meno rispett0 alle
% matrici triangolari piene.
%}
%% 1.2
n = 1000;
A = diagonals([-1 2 -1], n);
x = ones(n,1);
b = A*x;
c = rand(size(b));
c = c./norm(c);
db = c.*1e-6;

% per il calcolo dell'errore stimato, dalla teoria utilizzo il
% condizionamento in norma due di A e il residuo relativo normalizzato in norma 2:
% r_norm = norm(b-A*(x+dx))/norm(b) --> norm(-db)/norm(b) -->
% norm(db)/norm(b)
K2 = cond(A); % poiché il condizionamento è elevato la stima dell'errore sulla
              % base del res normalizzato non è soddisfacente: stimo un
              % errore molto più grande anche se ho un res
              % normalizzato molto piccolo
err_stim = K2*norm(db)/norm(b)

dx = (A\(b+db))-x;
err_vero = norm(dx)/norm(x) % infatti l'errore vero è tre ordini di grandezza
                            % inferiore a quello stimato
%% 1.3
n = 1000;
A = diagonals([-1 2 -1], n);
x = ones(n,1);
b = A*x;
[sol_g,N_g] = richardson(A,b,eye(n),b,1e-2,1e3);
N_g
x1_g = sol_g(1)
res_norm_g = norm(b-A*sol_g)/norm(b)
err_rel_g = norm(x-sol_g)/norm(x)
%% 1.4
[sol_gj,~,res_norm_gj,N_gj] = pcg(A,b,1e-2,1e3,[],[],b);
N_gj
x1_gj = sol_gj(1)
res_norm_gj 
err_rel_gj = norm(x-sol_gj)/norm(x)
% il metodo del pcg converge molto più velocemente, come ci si aspetta
% dalla teoria, in quanto le direzioni di discesa ad ogni iterazioni sono
% ottimali a tutte le precedenti, non solo a quella subito prima. Per
% questo converge più velocemente, con meno iterazioni (dalla teoria
% sappiamo anche che converge al massimo con n iterazioni), e residuo e
% errore relativo sono più piccoli
%% 1.5
n = 1000;
A = diagonals([-1 2 -1], n);
x = ones(n,1);
b = A*x;
l_max0 = eigpower(A,1e-6,0,ones(n,1))
l_max1 = eigpower(A,1e-6,1,ones(n,1))
[l_maxN] = eigpower(A,1e-6,1e3,ones(n,1))

%% 1.6
[l_maxN,~,~,lambda_vect] = eigpower(A,1e-6,1e3,ones(n,1));
l_max = 2*(1+cos(pi/1001));
err = abs(lambda_vect-l_max);
[p,c] = stimap(err);
ord = ceil(p(end))
err(end)/(err(end-1)^ord)
%% 1.7
n = 1000;
A = diagonals([-1 2 -1], n);
F = @(x) exp(-x./10) + A*x - ones(n,1);
J = @(x) (-diag(exp(-x./10))./10) + A;
x0 = ones(n,1);
for i = 1 : 3
    x = gaussnewton2(F,J,x0,1e-10,i);
    x(1)
end
%% 
syms g;
A = [8 -5 g;-5 12 3;g 3 8];
% (simplify(det(A)))
% f = @(g) - 12*g^2 - 30*g + 496;
% fsolve(f,-8)
% fsolve(f,5)
pretty(simplify(eig(A)))
%%
xmin= 2 ^ (-3-1);
epsm = 2 ^ (1-4);
xmin*epsm 