function [l,j] = ljindinv(a)
% for indexing wigner functions where the m index is ignored
% returns the l,j indices corresponding to the effective single index a

if (a<=0) 
    l=NaN;j=NaN;
else
    l = floor(sqrt(a-1));
    j = a - l^2-l-1;
end
end