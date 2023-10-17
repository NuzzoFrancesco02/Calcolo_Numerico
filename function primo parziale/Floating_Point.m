function [xmat,C,em,xmin,xmax,m]=Floating_Point(B,t,L,U)
% function [xmat,C,em,xmin,xmax,m]=Floating_Point(B,t,L,U)
%
% Questa funzione permette di calcolare l'insieme dei numeri floating point 
% positivi associato ai 4 numeri fondamentali che lo descrivono F=(β,t,L,U)
% Per ragioni di costo si è deciso di far valere il codice fino ad un
% valore di t pari a 10
%
% --------------------------------- INPUT ---------------------------------
%* B = base (numero intero che determina il sistema numerico)
%* t = numero di cifre  della mantissa
%* L = esponente minimo
%* U = esponente massimo
%
% -------------------------------- OUTPUT ---------------------------------
%* x = vettore contenente i numeri dell'insieme F (ordinati in senso crescente)
%* xmat = matrice contenente i numeri dell'insieme F [ xmat(m,e) ]
%* C = cardinalità dei soli numeri positivi
%* em = epsilon macchina
%* xmin = valore minimo
%* xmax = valore massimo
%* m = mantissa (in base 10)
%
% -------------------------------------------------------------------------


%% Controlli

if(t>10 || t<1)
    error("Valori di t non validi")
end


%% Valori Notevoli

em=B^(1-t);
xmin=B^(L-1);
xmax=B^U*(1-B^(-t));
e=L:U;


%% Determinazione della Mantissa

m=[];
for a1=1:B-1
    if (t<2)
        m=[m a1];
            
    else
        for a2=0:B-1
            if (t<3)
                m=[m a1*B+a2];

            else
                for a3=0:B-1
                    if (t<4)
                        m=[m a1*B^2+a2*B+a3];
                    
                    else
                        for a4=0:B-1
                            if (t<5)
                                m=[m a1*B^3+a2*B^2+a3*B+a4];
                                
                            else
                                for a5=0:B-1
                                    if (t<6)
                                        m=[m a1*B^4+a2*B^3+a3*B^2+a4*B+a5];
                            
                                    else
                                        for a6=0:B-1
                                            if (t<7)
                                                m=[m a1*B^5+a2*B^4+a3*B^3+a4*B^2+a5*B+a6];
                                                
                                            else
                                                for a7=0:B-1
                                                    if (t<8)
                                                        m=[m a1*B^6+a2*B^5+a3*B^4+a4*B^3+a5*B^2+a6*B+a7];
                                                        
                                                    else
                                                        for a8=0:B-1
                                                            if (t<9)
                                                                m=[m a1*B^7+a2*B^6+a3*B^5+a4*B^4+a5*B^3+a6*B^2+a7*B+a8];                                                             
                                    
                                                            else
                                                                for a9=0:B-1
                                                                    if (t<10)
                                                                        m=[m a1*B^8+a2*B^7+a3*B^6+a4*B^5+a5*B^4+a6*B^3+a7*B^2+a8*B+a9]; 
                                                        
                                                                    else
                                                                        for a10=0:B-1
                                                                            m=[m a1*B^9+a2*B^8+a3*B^7+a4*B^6+a5*B^5+a6*B^4+a7*B^3+a8*B^2+a9*B+a10];
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end


%% Calcolo dell'Insieme Numerico

for i=1:length(m)
    for j=1:length(e)
        
        xmat(i,j)=m(i)*B^(e(j)-t);
        
    end
end


% Cardinalità
C=(B-1)*B^(t-1)*(U-L+1);


end