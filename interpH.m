%% Calcola numero intervalli necessari per avere una certa tolleranza
% f: funzione
% tol : tolleranza che voglio raggiungere
% a : primo elemento dell'intervallo
% b : ultimo elemento dell'intervallo

function num_int = interpH(f,tol,a,b)

  I = [a:0.000001:b];
  f_1 = Jac(f,'[x]');
  f_2 = Jac(f_1,'[x]');
  max_f2 = max(abs(f_2(I)));
  H = sqrt((8*tol)./(abs(max_f2)));
  num_int = ceil((b-a)/(H));


end