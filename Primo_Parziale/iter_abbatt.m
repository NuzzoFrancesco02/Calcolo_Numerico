%% Abbattimento metodi iterativi
% A : matrice
% abb : fattore di abbattimento
% str : 'rich_staz','rich_din' o 'grad','grad_conj'
% P : matrice di precondizionamento

function k = iter_abbatt(A,abb,str,P)
    if nargin == 3
        P = eye(size(A,1));
    end
    if strcmp(str,'rich_staz') 
        if ((A~=A') | (P ~=P') )
            error('Matrici non simmetriche');
        end
        if ((min(eig(A))<=0) || (min(eig(P))<=0))
            error('Matrici non definite positive')
        end
        K = max(eig(P\A))/min(eig(P\A));
        d = (K-1)/(K+1);
        k = log(1/abb)/log(d);
    end
    if (strcmp(str,'rich_din') || strcmp(str,'grad'))
        if ((A~=A') | (P ~=P' ))
            error('Matrici non simmetriche');
        end
        if (min(eig(A))<=0 || min(eig(P))<=0)
            error('Matrici non definite positive')
        end
        K = max(eig(P\A))/min(eig(P\A));
        d = (K-1)/(K+1);
        k = log(1/abb)/log(d);
    end
    if strcmp(str,'grad_conj')
        if (A~=A' | P ~=P' )
            error('Matrici non simmetriche');
        end
        if (min(eig(A))<=0 || min(eig(P))<=0)
            error('Matrici non definite positive')
        end
        K = max(abs(eig(P\A)))/min(abs(eig(P\A)))
        c = (sqrt(K)-1)/(sqrt(K)+1);
        k
    end
    k = ceil(k);
end
