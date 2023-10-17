%% lab 17 es 1
clc
clear
close all

%% punto 1
B = [10 -1 1 0; 1 1 -1 3; 2 0 2 -1; 3 0 1 5]
gershcircles (B)


%% punto 2
eig(B)

%% punto 3
[lambda,x,iter]=invpower(B,1e-6,1000,ones(4,1));

%% punto 4
% shift +i , -i
[lambda,x,iter]=invpowershift(B,i,1e-6,1000,ones(4,1));
[lambda,x,iter]=invpowershift(B,-i,1e-6,1000,ones(4,1));

% lambda max, parto dal bordo dei cerchi di gersh
[lambda,x,iter]=invpowershift(B,12,1e-6,1000,ones(4,1));
% confronto con potenze
[lambda,x,iter]=eigpower(B,1e-6,1000,ones(4,1));