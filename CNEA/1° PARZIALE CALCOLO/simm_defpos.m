%Controllo se la matrice è simmetrica e definita positiva
function boolean = simm_defpos(M)
    boolean = true;
    if (isequal(M,M'))
        % è simmetrica
        for i = 1 : length(diag(M))
            if(det(M((1 : i), (1 : i))) < 0)
                boolean = false;
                % trovo una sottomatrice principale con det negativo
            end
        end
    else 
        error('La matrice non è simmetrica!');
    end
    if (boolean == false)
        error('La matrice non è definita positiva!');
    end
end