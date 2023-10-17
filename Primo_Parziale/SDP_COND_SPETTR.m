%% SDP COND SPETTR
it = 0;
toll = 1e-8;
err = toll + 1; % 1
nmax = 1000;

y_max = x0/norm(x0); % 3n
y_min = y_max;
K = 0;

R = chol(A); %% (n^3)/3;
while it < nmax && err > toll
    x_max = A*y_max;                %   2n^2-n
    y_max = x_max/norm(x_max);      %   3n
    % R'*R *x_min = y_min
    f = fwsub(R',y_min);            %   n^2
    x_min = bksub(R,f);             %   n^2
    y_min = x_min/norm(x_min);      %   3n
    Kv = K;
    K = (y_max'*A*y_max)/(y_min'*A*y_min); %    (2n-1) + (2n-1) + (2n^2-n)
    it = it + 1;                    %   1
    err = abs(K-Kv);                %   1
end
