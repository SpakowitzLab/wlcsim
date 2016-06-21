function M = shearWLCpropagator(k,L,lp,g,epari,ephi,eta,LMAX)
% get structure factor at a particular k value
% for a shearable WLC 
% = InverseLaplace(G(k,p))
% get the propagator matrix for a shearable WLC of length L

% get the matrix
A = shearWLCgetAmat(k,lp,g,epari,ephi,eta,LMAX);
% get the poles
evals = eig(-A);
[tmp1,tmp2] = sort(real(evals),'descend');
evals = evals(tmp2);

% get the exponential components
M = zeros(LMAX+1,LMAX+1);
for ec = 1:LMAX+1
    cA = cofactor(evals(ec)*eye(LMAX+1)+A);
    p = evals(ec);
    expterm = cA/(prod(p-evals(1:ec-1))*prod(p-evals(ec+1:end)));
    M = M + expterm*exp(evals(ec)*L);
end

end