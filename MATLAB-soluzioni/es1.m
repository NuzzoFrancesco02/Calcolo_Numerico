%% Laboratorio 15, Esercizio 1
clc
clear
close all

%% punto 1 
A=randi([0,1],100,100);
s=sum(A);
for j=1:size(A,1)
    A(j,:)=A(j,:)./s;
end


%% punto 2
B=[
   0   0   0   1   0
   1   0   0   0   0
   0   1   0   0   0
   0   1   0   0   1
   1   1   1   1   0];
s=sum(B);
for j=1:size(B,1)
    B(j,:)=B(j,:)./s;
end

%% punto 4
[lA,xA,iterA]=eigpower(A,1.e-06,100,1/100*ones(100,1));
lA

[lB,xB,iterB]=eigpower(B,1.e-06,100,1/5*ones(5,1));
lB
xB