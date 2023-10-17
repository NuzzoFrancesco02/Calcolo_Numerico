function []=dom_stretta_righe(A)
    for i=1:size(A,1)
        s=0;
        for j=1:size(A,1)
           s=s+abs(A(i,j));
        end
        s=s-abs(A(i,i));
        if(abs(A(i,i))<=s)
            error('la matrice non è a dominanza stretta per righe\n')
        end
    end
    disp('la matrice è a dominanza stretta per righe')
end
