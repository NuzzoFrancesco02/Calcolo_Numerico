function[B]=matrice_iterazione(A,P)
    B=eye(size(A,1))-inv(P)*A;
end