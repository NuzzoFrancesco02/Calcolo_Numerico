function [mu,x,min_eigenvec_comp]=invpower_mod(A,nmax,x0)
%INVPOWER    Approssima l'autovalore di modulo minimo
%            di una matrice.
%   [MU,X,MIN_EIGVEC_COMP] = INVPOWER_MOD(A, NMAX, X0) esegue NMAX iterazioni 
%   dell'algoritmo del metodo delle potenze inverse, a partire dal dato 
%   iniziale X0. Gli output sono;
%   MU: minimo autovalore calcolato dopo NMAX iterazioni,
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

% calcolo la fattorizzazione LU una volta per tutte
[L,U,P]=lu(A);

% iterazione zero fuori dal ciclo while
iter = 0;
y = x0/norm(x0); % y0
min_eigenvec_comp = abs(y'*eigvec/norm_eigvec2);
mu = y'*A*y; % mu0

while abs(mu) ~= 0 && iter<nmax
   iter = iter + 1;
   % risolvo Ax^{(k)}=y^{(k-1)}
   z=fwsub(L,P*y);
   x=bksub(U,z);
   y= x/norm(x);
   min_eigenvec_comp = [min_eigenvec_comp, abs(y'*eigvec/norm_eigvec2)];
   munew = y'*A*y;
   err = abs(munew - mu);
   mu = munew; 
end


return
