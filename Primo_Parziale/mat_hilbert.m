function M = mat_hilbert(n)
%Crea una matrice hilbertiana di ordine n
M = zeros(n);
if(n > 0)
    for i = 1 : n
        for j = 1 : n
            M(i,j) = 1/(i+j-1);
        end
    end
else
    fprintf('\nDimensione matrice sbagliata!\n');
end
return