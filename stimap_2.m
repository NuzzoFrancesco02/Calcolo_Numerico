function [p,c] = stimap_2(E,H)
    if length(E)~=length(H)
        error('Dimensioni errate')
    end
    p = [];
    c = [];
    for i = 2:length(E)
        p = [p log(E(i)/E(i-1))/log(H(i)/H(i-1))];
        c = [c E(i)/(E(i-1)^p(end))];
    end
    fprintf('\n\t###########');
    fprintf('\nOrdine di convergenza: %.4f',p(end));
    fprintf('\nFattore di convergenza: %e',c(end));
    fprintf('\n\t###########\n');