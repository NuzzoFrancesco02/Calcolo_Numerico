%% Lab 17 Esercizio 2
clc
clear
close all

%% punto 2
A1 = [30 2 3 13; 5 11 10 8; 9 7 6 12 ; 4 14 15 1];
D1 = qrbasic(A1, 1e-10, 1000)  % 41 iterazioni
A2 = [-30 2 3 13; 5 11 10 8; 9 7 6 12 ; 4 14 15 1];
D2 = qrbasic(A2, 1e-10, 1000)  % 843 iterazioni

for i=1:length(D1)-1
    v1(i) = abs(D1(i+1)/D1(i));
end
V1 = max(v1)

for i=1:length(D2)-1
    v2(i) = abs(D2(i+1)/D2(i));
end
V2 = max(v2)