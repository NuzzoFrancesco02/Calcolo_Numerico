%% INTERPOLAZIONE TRIGONOMETRICA
% Ritorna l'equazione in SIMBOLICO dell'interpolazione trigonometrica
function I_f = poly_trigo(f,n)
syms x;
h = 2*pi/(n+1);
if mod(n,2)==0
    M = n/2;
    mu = 0;
else
    M = (n-1)/2;
    mu = 1;
end
S = 0;
ck_tilde = [];
j = (0:n);
xj = j.*h;
for k = -(M+mu):(M+mu)
    ck = (sum(f(xj).*exp(-1i.*k.*j.*h)))/(n+1);
    if k == -(M+1) || k == (M+1)
        ck_tilde = ck/2;
    else
        ck_tilde = ck;
    end
    S =  S + simplify(ck_tilde.*(cos(k*x)+i*sin(k*x)));
end

digits(5)
Svpa = vpa(S,4);
I_f = Svpa;
