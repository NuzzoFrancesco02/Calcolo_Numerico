%% Autovalore approssimato
%   l = eig_approx(A,x,k)
%   ATTENZIONE!!!: la prima approssimazione si ha con 0: l(0) = prima
%   approssimazione
%
%   # OUTPUT
%   l : autovalore relativo
%
%   # INPUT
%   l = eig_approx(A,x,k)
%   A : matrice
%   x : autovettore approssimato
%   k : ultima iterata 

function l = eig_approx(A,x,k)
if nargin == 2
    k = 0;
end
it = -1;
while it < k
    y = x/norm(x);
    l = y' * A * y;
    x = A * y;
    it = it + 1;
end