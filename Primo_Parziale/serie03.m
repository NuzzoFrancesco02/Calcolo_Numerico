% serie03
%% 1
% 1
n = 100;
A = randi([0 1],n);
s = sum(A);
for i = 1 : n
    A(:,i) = A(:,i)/s(i);
end
% 2
B = [0 0 0 1 0; 1 0 0 0 0; 0 1 0 0 0; 0 1 0 0 1; 1 1 1 1 0];
s = sum(B);
for i = 1 : 5
    B(:,i) = B(:,i)/s(i);
end
[l1, x1, k1] = eigpower(A,1e-6,100,ones(n,1)./n);
[l2, x2, k2] = eigpower(B,1e-6,100,ones(5,1)./5);

%% 2
A = toeplitz(1:4);
[l1, x1, iter] = invpower(A, 1e-6, 100, [1 2 3 4]' );
[l2, x2, iter] = invpower(A, 1e-6, 100, [1 1 1 1]' ); % si hanno problemi quando l'inizializzazione 
% è perpendicolare al vettore che stiamo cercando e restituisce il secondo
% minore in modulo

%% 3
A = [10 -1 1 0; 1 1 -1 3;2 0 2 -1; 3 0 1 5];

gershcircles(A)
figure(4)
compass(eig(A))
l = invpower(A) % Non converge, due autovalori con modulo uguale
l1 = invpowershift(A,i) % Converge all'autovalore con P.I > 0
l2 = invpowershift(A,-i) % Converge all'autovalore con P.I < 0

%% 4
a = -30;
A = [a 2 3 13; 5 11 10 8; 9 7 6 12; 4 14 15 1];
L = eig(A)
% qr si può utilizzare solo se hanno moduli distinti, se i moduli di due
% autovalori adiacenti sono molto vicini fra loro il metod converge molto
% più lentamente
D = qrbasic(A,1e-10,1000)

%% 5
IM = imread("ellie.png");
IM = rgb2gray(IM);
imshow(IM);
A = double(IM);
k = 2;
[U,S,V] = svd(A);


for i = 1 : 1
    k = 500*i;
    Ak = U(:,1:k)*S(1:k,1:k)*(V(:,1:k)');
    %subplot(2,2,i);
    IM_k = uint8(Ak);
    imshow(IM_k);
end
m = size(A,1);
n = size(A,2);
size1= m*n;
size2 = k + k*m + k*n;
sigma = diag(S);
sigma_out = diag(S(k+1:end,k+1:end));
err = sqrt(sum(sigma_out.^2));
err_rel = err / sqrt(sum(sigma.^2))


