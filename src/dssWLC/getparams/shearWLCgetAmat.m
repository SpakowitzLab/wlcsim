function A = shearWLCgetAmat(k,eb,g,epari,ephi,eta,LMAX)
% get the A matrix for a continuous shearable WLC at a particular k value
% uses parameters as defined in the manuscript

lp = eb/(1+eta^2*eb*ephi);

alj = @(l,j) sqrt((l-j)*(l+j)/((2*l-1)*(2*l+1)));
blj = @(l,j) alj(l,j)*alj(l-1,j);
clj = @(l,j) alj(l,j)^2 + alj(l+1,j)^2;
A = zeros(LMAX+1,LMAX+1);
for l = 0:LMAX
    lind=l+1;
    A(lind,lind) = l*(l+1)/2/lp +k^2/2*ephi*(1-clj(l,0)) + k^2/2*epari*clj(l,0);
    if (l < LMAX)
        A(lind,lind+1) = -i*g*k*alj(l+1,0) - i*k*eta*ephi*l*alj(l+1,0);
    end
    if (l>0)
        A(lind,lind-1) = -i*g*k*alj(l,0) + i*k*eta*ephi*(l+1)*alj(l,0);
    end
    if (l < LMAX-1)
        A(lind,lind+2) = -k^2/2*ephi*blj(l+2,0) + k^2/2*epari*blj(l+2,0);
    end
    if (l>1)
        A(lind,lind-2) = -k^2/3*ephi*blj(l,0) + k^2/2*epari*blj(l,0);    
    end
end

end