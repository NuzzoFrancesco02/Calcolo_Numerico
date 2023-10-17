function a_bar = armonica(a)
  n = numel(a);
  a_bar = (1/n * sum(a.^(-1)))^(-1);
end     
