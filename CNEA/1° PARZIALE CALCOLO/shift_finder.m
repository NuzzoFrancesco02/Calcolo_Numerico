%% s = shift_finder(A,l)
%
%   Per quale valore dello shift s posso approssimare l'autovalore l,
%     ATTENZIONE: 
%     • se s ∈ R, non ci devono essere autovalori simmetrici all'asse reale
%     • se s ∈ Im, non ci devono essere
%       autovalori simmetrici all'asse immaginario
%   # INPUT
%   A : matrice assegnata
%   l : autovalore da approssimare con shift inverso
%
% SE L'ESERCIZIO CHIEDE μ NON C'È BISOGNO DI METTERE i, μ È GIÀ PARTE
% IMMAGINARIA.

function  s = shift_finder(A,l)
    k = find(round(l,4)==round(eig(A),4),1);
    if isempty(k)
        error("L'autovalore selezionato non è autovalore della matrice!")
    end
    if imag(l) ~= 0 
        str = 'imag';
    else
        str = 'real';
    end
    [m,n] = size(A);
    if m ~= n 
        error('La matrice non è quadrata');
    end
    
        flag = eig(A);
        k = 0;
        i = 1;
        while  i <= m && k < m-1
            if flag(i) ~= l
               lambda(k+1) = flag(i);
               k = k + 1;
            end
            i = i + 1;
        end
        [~, i] = min(abs(lambda)); % |λ-s| è più vicino ad s se prendo λs più piccolo
        l_min = lambda(i);
        if strcmp(str,'real') & (real(l)==round(min(real(eig(A))),4) | real(l)==round(max(real(eig(A))),4))
        s = ( abs(l_min)^2 - abs(l)^2 )  / (2*(real(l_min)-real(l)));
            if real(l)<real(l_min)
                if l < s 
                    fprintf('\ns < %f  &  s ~= %d\n',s,l);
                else
                    fprintf('\ns < %f\n',s);
                end
            else
                if l > s
                    fprintf('\ns > %f  &  s ~= %d\n',s,l);
                else
                    fprintf('\ns > %f\n',s);
                end
            end
        x = sym(s)
        elseif strcmp(str,'imag') & (imag(l)==round(min(imag(eig(A)),4)) | imag(l)==round(max(imag(eig(A))),4))
            s = ( abs(l_min)^2 - abs(l)^2 )  / (2*(imag(l_min)-imag(l)));
            l_min = conj(l);
            if imag(l)<imag(l_min)
                    fprintf('\nμ*i< %f*i\n',s);
            else
                    fprintf('\nμ*i > %f*i\n',s);
            end
        x = sym(s)
        %elseif m ~= 3
            %error('Ti attacchi al cazzo!');
        elseif strcmp(str,'real')
            lambda = eig(A);
            [~,H] = sort(sign(real(lambda)).*abs(lambda));
            k = find(round(lambda,4)==round(l,4),1);
            j = find(H==k,1);
            l_min = lambda(H(j-1));
            s_min = ( abs(l_min)^2 - abs(l)^2 )  / (2*(real(l_min)-real(l)));
            l_max = lambda(H(j+1));
            s_max = ( abs(l_max)^2 - abs(l)^2 )  / (2*(real(l_max)-real(l)));
            fprintf('\n%f < s < %f  &  s ~= %d\n',s_min,s_max,l);
        elseif strcmp(str,'imag')
            lambda = eig(A);
            [~,H] = sort(sign(imag(lambda)).*abs(lambda));
            k = find(round(lambda,4)==round(l,4),1);
            j = find(H==k,1);
            l_min = lambda(H(j-1));
            s_min = ( abs(l_min)^2 - abs(l)^2 )  / (2*(imag(l_min)-imag(l)));
            l_max = lambda(H(j+1));
            s_max = ( abs(l_max)^2 - abs(l)^2 )  / (2*(imag(l_max)-imag(l)));
            fprintf('\n%f < s < %f  &  s ~= %d\n',s_min,s_max,l);
        end
    