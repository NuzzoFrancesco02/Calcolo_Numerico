function [lambda]=rayleigh(A,x)
%funzione per l'approssimazione dell'autovalore
%utilizzando il quoziente di rayleigh avendo l'approssimazione
%dell'autovettore x della matrice A
lambda=(x'*A*x)/(x'*x);
end