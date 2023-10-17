clear
clc

IM = imread("street1.jpg");	%carico l'immagine
IM = rgb2gray(IM);	%converto l'immagine in scala di grigi
%imshow(IM)		%visualizzo l'immagine
A = double(IM);	%converto in formato a virgola mobile

[n,m] = size(A);
p = min(n,m)
[U,S,V] = svd(A);
k = min(10,p);
Ak = U(:,1:k) * S(1:k,1:k) * (V(:,1:k))';
Ak = uint8(Ak);
%imshow(Ak)
dia = diag(S);
loglog(1:p, diag(S), '--b', [k k], [S(p,p) S(1,1)], '--r', 'LineWidth', 2);
grid on;
xlabel('i');
ylabel('\sigma_i');
title('\sigma_i vs i');

for i = 1:p-2
	uncompressed_size = m*n;
	compressed_size = m*i + i + i*n;
	compression_ratio(i) = uncompressed_size / compressed_size;
	relative_error(i) = sqrt(sum(diag(S(i+1:end, i+1:end)).^2) / sum(diag(S).^2));
end

figure(1);
clf;
semilogy(compression_ratio,'b');
hold on;
semilogy(relative_error,'r');