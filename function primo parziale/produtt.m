function out = produtt(v)

out = v(1);
for n = 2:length(v)
out = v(n)*out;
end