%% x = eps_m(B,a,b,str)
% str = 'eps' ---> b = eps macchina
% str = 't' ---> b = t, cifre mantissa 
% str = 'xmin'
% str = 'xmax_U'
% str = 'xmax_t'
function x = eps_m(B,a,str)
if nargin == 3
    b = 0;
end
if strcmp(str,'t')
    x = ceil(1-log(a)/log(B));
elseif strcmp(str,'eps')
    x = B^(1-a);
elseif strcmp(str,'xmin')
    x = 1 + log(a)/log(B);
    if x > 0
        x = ceil(x);
    else
        x = floor(x);
    end
elseif strcmp(srt,'xmax')

end