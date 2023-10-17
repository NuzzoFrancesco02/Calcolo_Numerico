function alpha = vectorangle(v,w)

% Angolo tra vettori v e w in gradi

alpha = acos((v'*w)/(norm(v)*norm(w)));
alpha = rad2deg(alpha);
end