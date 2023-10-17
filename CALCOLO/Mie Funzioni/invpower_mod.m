function [lambda,x,min_eigenvec_comp] = invpower_mod(A,x0,nmax)
%INVPOWER    Approssima l'autovalore di modulo minimo
%            di una matrice.
%   [LAMBDA,X,MIN_EIGVEC_COMP] = INVPOWER_MOD(A, X0, NMAX) esegue NMAX iterazioni 
%   dell'algoritmo del metodo delle potenze inverse, a partire dal dato 
%   iniziale X0. Gli output sono;
%   LAMBDA: minimo autovalore calcolato dopo NMAX iterazioni,
%   X: autovettore associato a MU ottenuto dopo NMAX iterazioni,
%   MIN_EIGVEC_COMP: vettore di lunghezza NMAX+1. Detta y_i
%   l'approssimazione di X ottenuta alla i-esima iterazione,
%   MIN_EIGVEC_COMP(i) contiene la componente di y_i nella direzione
%   definita dall'autovettore associato all'autovalore minimo (in modulo)
%   della matrice A.

[V,D] = eig(A);
D = diag(D);
[~,ind] = min(abs(D));
eigvec = V(:,ind);
norm_eigvec2 = norm(eigvec)^2;

[n,m] = size(A);
if n ~= m, error('Solo per matrici quadrate'); end

flag = '';
if nargin < 3
	nmax = 100;
end
if nargin < 2
	x0 = ones(n,1);
	flag = ', default vector';
end

fprintf("\nMetodo delle potenze inverse modificato (nnmax = %d%s)\n", nmax, flag);

% calcolo la fattorizzazione LU una volta per tutte
[L,U,P]=lu(A);

% iterazione zero fuori dal ciclo while
iter = 0;
y = x0/norm(x0); % y0
min_eigenvec_comp = abs(y'*eigvec/norm_eigvec2);
lambda = y'*A*y; % mu0
x = x0;
while abs(lambda) ~= 0 && iter<nmax
   iter = iter + 1;
   % risolvo Ax^{(k)}=y^{(k-1)}
   z=fwsub(L,P*y);
   x=bksub(U,z);
   y= x/norm(x);
   min_eigenvec_comp = [min_eigenvec_comp, abs(y'*eigvec/norm_eigvec2)];
   lambdanew = y'*A*y;
   err = abs(lambdanew - lambda);
   lambda = lambdanew; 
end

fprintf('Il metodo delle potenze inverse modificato dopo %d iterazioni tende all''autovalore\n', iter);
lambda

return
end