function []=dom_stretta_colonne(A)
    for i=1:size(A,1)
        s=0;
        for j=1:size(A,1)
           s=s+abs(A(j,i));
        end
        s=s-abs(A(i,i));
        if(abs(A(i,i))<=s)
            error('la matrice non è a dominanza stretta per colonne\n')
        end
    end
    disp('la matrice  è a dominanza stretta per colonne')
end
