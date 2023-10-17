%% lab 15, esercizio 2
clc
clear 
close all

%% punto 2
A = toeplitz (1:4);
invpower (A, 1e-6, 100, [1 2 3 4]')

invpower (A, 1e-6, 100, ones (4, 1))


% Utilizziamo la funzione invpower_mod, che non richiede com input la 
% tolleranza, bensi' esegue un numero fissato di iterazioni 
% (in questo caso 100).

[lambda, x, min_eigenvec_comp] = invpower_mod (A, 100, ones (4, 1));

lambda


% Nel seguente grafico si riporta l'andamento della componente di y_i 
% (ossia le approssimazioni dell'autovettore calcolate durante le 
% iterazioni di invpower) lungo la direzione definita dall'autovettore 
% associato al minimo autovalore di A.
semilogy(min_eigenvec_comp)