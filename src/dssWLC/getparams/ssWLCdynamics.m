function [evals,evecs,pareval,Amat] = ssWLCdynamics(eb,g,epar,eperph,eta,L,xir,xiu,pmax)
% get the eigenvalues and eigenvectors for continuous shearable WLC
% dynamics, using mode decomposition of the coupled shear and bend
% pareval is the eigenvalue for stretch

if (nargin<9)
    pmax = 100;
end

Amat = zeros(2*pmax);
% for each p, lists X1,O2
for p = 1:pmax
    x1i = 2*(p-1)+1;    
    o2i = 2*(p-1)+2;
    
    coeff = pi^2*p^2/L^2;
    
    Amat(x1i,x1i) = -(eperph/xir*coeff + g^2*eperph/xiu);
    Amat(x1i,o2i) = eta*eb/xir*coeff+eta*eb*g^2/xiu;      
    Amat(o2i,o2i) = (-eb/xiu*coeff);
    Amat(o2i,x1i) = eta*eb/xiu*coeff;
    
    for pp = 1:pmax
        if (mod(pp+p,2)==1)
            x1ip = 2*(pp-1)+1;           
            o2ip = 2*(pp-1)+2;
            
            tmp = 4*g/L/xiu*(p*pp)/(p^2-pp^2);
            
            Amat(x1i,x1ip) = Amat(x1i,x1ip) +tmp*eta*eb;
            Amat(x1i,o2ip) = Amat(x1i,o2ip) -tmp*eb;                        
            Amat(o2i,o2ip) = Amat(o2i,o2ip) -tmp*eta*eb;
            Amat(o2i,x1ip) = Amat(o2i,x1ip) +tmp*eperph;
           % [p,pp,x1i,o2ip,o2i,x1ip, Amat(x1i,o2ip),Amat(o2i,x1ip)]
        end
    end
end
[V,E] = eig(Amat); E = diag(E);
[a,b] = sort(real(E));
evals = E(b);
evecs = V(:,b);

pareval = -epar/xir*pi^2/L^2;

end