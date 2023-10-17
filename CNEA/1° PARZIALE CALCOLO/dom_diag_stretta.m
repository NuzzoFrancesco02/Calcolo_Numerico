% Controllo se la matrice è a dominanza diagonale stretta per righe o per

function [boolean_righe, boolean_colonne] = dom_diag_stretta(M)
n = length(diag(M));
diag_M = diag(M);
% Dominanza stretta per colonne
boolean_colonne = false;
    for i = 1 : n
        r(i) = abs(diag_M(i)) - sum(abs( M( [1:i-1,i+1:n], i )));
    end
    if ( r > 0)
        boolean_colonne = true;
    end
% Dominanza stretta per righe
boolean_righe = false;
    for i = 1 : n
        c(i) = abs(diag_M(i)) - sum(abs (M (i,[1:i-1,i+1:n])));
    end
    if (c > 0)
        boolean_righe = true;
    end


    if (boolean_righe == false && boolean_colonne == false)
        error ('La matrice non è a dominanza stretta né per righe né per colonne');
    end
end