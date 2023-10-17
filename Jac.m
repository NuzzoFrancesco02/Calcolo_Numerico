function [J] = Jac(F,var)
% Questa funzione, presi come input una funzione F (sia simbolica che anonymous) e il vettore "var" (di tipo char) delle 
% variabili da cui dipende F (che si può anche non inserire), calcola la Jacobiana di F.
%
%
% SINTASSI: [J] = Jac(F,'[var1,var2,var3,...]')    dove var1,var2,var3,... sono le variabili; per esempio var = '[x,y,z]' 
% 
%   
%                                             !ATTENZIONE! 
%- Qualora la funzione F di input fosse un vettore di funzioni anonymous, la Jacobiana restituita come output 
%  sarà anch'essa un vettore/una matrice di funzioni anonymous.
%- Quaolora la funzione F di input fosse un vettore di funzioni simboliche (sym), la Jacobiana restituita come output 
%  sarà anch'essa un vettore/una matrice di funzioni simboliche (sym).
%
%
%                                             !ATTENZIONE!
% La funzione F di input DEVE essere un vettore colonna e NON deve essere un vettore nullo.
% 
% 
%                                             !ATTENZIONE!
% Se non si specificasse il secondo parametro di input della funzione, ossia le variabili da cui dipende F, si potrebbero
% verificare due tipologie di problemi:
%
% 1) Qualora si inserisse come input una funzione dipendente da più di una variabile, la Jacobiana restituita come output 
%    dipenderà dalle medesime variabili, ma ordinate secondo un ordine alfabetico (al più può variare la posizione delle 
%    colonne della Jacobiana). 
%    Dunque, nell'ipotesi che la funzione F di input dipenda dalle variabili "c,a,b", accade che la chiamata: 
%    [J] = Jac(F(c,a,b)) avrà output J = J(a,b,c) e NON J = J(c,a,b).
%
% 2) Qualora la funzione F di input dipenda da n variabili, ma nella sua scrittura compaiano n-m variabili (con m<n), allora 
%    se non si specificasse il secondo parametro di input, la Jacobiana restituita come output dipenderà da n-m variabili e 
%    non dalle n variabili dichiarate inizialmente.
%    Dunque, nell'ipotesi che la funzione F di input sia dichiarata nel seguente modo: F = @(x,y,z) [x+3*y; x-y], 
%    (con n=3 e m=1 per ricondursi al caso generale sopra citato) allora accade che la chiamata:
%    [J] = Jac(F) avrà output J = @(x,y)[1,3;1,-1] e NON J = @(x,y,z)[1,3,0;1,-1,0]
%   
%
%                                             !ATTENZIONE!
% Qualora si inserisca il secondo parametro di input, ossia le variabili da cui dipende F, è importante non commettere errori.
% Infatti se le variabili vengono invertite, oppure se non vengono inserite tutte, o se ne vengono inserite più di quelle 
% effettive l'output restituito potrebbe dare risultati inattesi.
%
%
%
%                                 STRUTTURA FUNZIONE "F" E JACOBIANA "J" 
%
%
%                  /               \                          /                           \
%                  | f1(x1,...,xn) |                          | df1/dx1 - - - - - df1/dxn |
%                  |       -       |                          |    -    -           -     |
%              F = |       -       |    -->    J(x1,...,xn) = |    -       -        -     |
%                  |       -       |                          |    -           -    -     |
%                  |       -       |                          |    -                -     |
%                  | fm(x1,...,xn) |                          | dfm/dx1 - - - - - dfm/dxn |
%                  \               /                          \                           /
%
%
%
% Esempi:
%
%   Esempio 1: Calcolo della Jacobiana di una F dipendente da 2 variabili
%       F = @(x,y) [exp(x)+cos(y);sin(x)*sin(y)];
%       J = Jac(F);           % ritorna il function_handle  @(x,y)[exp(x),-sin(y);cos(x).*sin(y),cos(y).*sin(x)]
%       
%       % In modo del tutto equivalente si può chiamare la funzione nel seguente modo:
%       J = Jac(F,'[x,y]');   % ritorna il function_handle  @(x,y)[exp(x),-sin(y);cos(x).*sin(y),cos(y).*sin(x)]
%
%
%   Esempio 2: Calcolo della Jacobiana di una F dipendente da 1 variabile (di fatto si calcola la derivata)
%       F = @(x) exp(x).*tan(x);
%       J = Jac(F);           % ritorna il function_handle  @(x)exp(x).*(tan(x)+tan(x).^2+1)
%
%
%   Esempio 3: Data una F dipendente da 2 parametri non scritti in ordine alfabetico, l'output di Jac cambierà a 
%       seconda di come si chiama la funzione (inserendo o meno il secondo parametro di input)
%       F = @(y,x) [x+2*y; 3*x+4*y];      
%       J = Jac(F);           % ritorna il function_handle   @(x,y)[1,2;3,4]
%                             % di fatto questa chiamata a Jac è equivalente a J = Jac(F,'[x,y]')
%       
%       % Invece, scrivendo le variabili da cui dipende F in modo ordinato:
%       J = Jac(F,'[y,x]');   % ritorna il function_handle  @(y,x)[2,1;4,3]
%
%
%   Esempio 4: Data una F dipendente da 3 parametri, ma scritta utilizzandone 2, l'output di Jac cambierà a seconda 
%       di come si chiama la funzione (inserendo o meno il secondo parametro di input)
%       F = @(x,y,z) [2*x+y; x-y];
%       J = Jac(F);           % ritorna il function_handle  @(x,y)[2,1;1,-1]
%                             % di fatto questa chiamata a Jac è equivalente a J = Jac(F,'[x,y]')
%
%       % Invece, scrivendo tutte le variabili dichiarate nella creazione di F:
%       J = Jac(F,'[x,y,z]'); % ritorna il function_handle  @(x,y,z)[2,1,0;1,-1,0]
%
% Autore: Emanuele Bianco - 17 maggio 2022
% emanuele2.bianco@mail.polimi.it 


    isAnonymous = 0;
    
    % Se l'input è una funzione, la converto in sym
    if (isa(F,"function_handle")==1)
        % Controllo se la funzione in input presenta errori di scrittura
        try
            cellVar = num2cell(ones(1,nargin(F)));          % Vettore cella di "1" di dimensione 1xn, con n = numero di variabili dichiarate, ossia quelle in @(....)
            F(cellVar{:});                                  % Valuto la funzione nel vettore cellVar per verificare che non compaiano errori
        catch
            error("Scrivere correttamente la funzione!");   % Qualora la valutazione della funzione nel vettore cellVar abbia generato un errore
        end
        % Creazione matrice simbolica
        F = str2sym(char(F));
        isAnonymous = 1;
    end


    % CONTROLLO CONSISTENZA F 
    if size(F,2)~=1
        error("La funzione F deve essere un vettore colonna!");
    end
    
    % ERRORE SE F É IL VETTORE NULLO
    if F==zeros(size(F,1),1)
        error("La funzione F non deve essere un vettore nullo!");
    end

    
    if (nargin == 1)
        % Calcolo Jacobiana --> variabile sym
        J = jacobian(F,symvar(F));  %symvar(F) è il vettore sym contenente le variabili da cui dipende F: [x, y, z, ...]
        J = simplify(J);
                                    
    
        % Se l'input è una funzione, converto la Jacobiana simbolica in una Jacobiana anonymous
        if isAnonymous == 1

            % Creo la stringa contenente le variabili da cui dipende F
            var = symvar(char(F));    % symvar(char(F)) è il vettore cella della forma: [{'x'}; {'y'}; {'z'}; ...]
    
            % Concateno la stringa delle variabili alla stringa della matrice Jacobiana convertita da sym in char
            strJacobian = ['@(', sprintf('%s,',var{1:end-1}), var{end}, ') ', char(J)];    % Assume la forma: ['@(x,y,...) riga1; riga2; ...']
                                                                                  
            % Converto il vettore char in una funzione anonymous (inserendo i "." prima di "/", "*", "^")
            J = str2func(vectorize(strJacobian));
        end


    else
        % Calcolo Jacobiana --> variabile sym
        J = jacobian(F,str2sym(var));  %str2sym(var) serve a convertire la stringa delle variabili di F in un vettore sym della forma: [x, y, z, ...]
        J = simplify(J);
                                      
    
        % Se l'input è una funzione, converto la Jacobiana simbolica in una Jacobiana anonymous
        if isAnonymous == 1

            % Creo la stringa contenente le variabili da cui dipende F (ordinate come l'utente ha specificato)                                                
            varAnonymous = ['@(' var(2:end-1) ') '];      % var ha la forma '[x,y,z,...]' ==> "salvo" la parte centrale escludendo le parentesi quadre

            % Concateno la stringa delle variabili alla stringa della matrice Jacobiana convertita da sym in char 
            strJacobian = [varAnonymous, char(J)];  % Assume la forma: ['@(x,y,...) riga1; riga2; ...']
                                                    
            % Converto il vettore char in una funzione anonymous (inserendo i "." prima di "/", "*", "^")
            J = str2func(vectorize(strJacobian));
        end
    end
end