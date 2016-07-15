function [xivals,tmp] = expandFsph(del,eb,alpha,LMAX)
% for the dssWLC model structure factor and moment calculations,
% expand the function F(uhat) in terms of spherical harmonics of uhat
% del is the discretization length
% eb is the bending modulus
% eperpi is 1/shear modulus
% eta is the bend-shear coupling
% WARNING: this becomes unstable at higher indices for higher eperpih

a = eb/del;
b = alpha*eb/(1+alpha)/2/del;
tmp = b/a^2;
Ivals = zeros(LMAX+1,1);

%for n = 0:LMAX
%    Ivals(n+1) = quad(@(x) x.^n.*exp(a*x-b*x.^2),-1,1,1e-12);
%end
% 
% %%

cutoff = 0.05e-2;
%[b/a^2, cutoff]
% below this cutoff, use an expansion of the quadratic part of the integral
if (b/a^2 < cutoff)
    %%
    % get the J integrals (incomplete gamma functions)
    % rescale all J by multiplying by exp(-a)/a^m
    pmax = ceil(log(eps)/log(b/a^2)+5); % maximum power to go to
    indmax = 2*pmax+LMAX;
    Jvals = zeros(2*pmax+LMAX+1,1);
    Jvals(1) = (1 - exp(-2*a));
    for m = 1:indmax
        Jvals(m+1) = -m/a*Jvals(m) + (1-(-1)^m*exp(-2*a));
    end
    
    % get the I_n integrals, multiplied by exp(-a)
    ks = 0:pmax;
    for n = 0:LMAX
        Ivals(n+1) = sum(1/a*(-1).^(ks).*b.^ks.*Jvals(n+2*ks+1)'./factorial(ks));
    end
%     %%
else
     %%
%     for n = 0:LMAX
%         Ivals(n+1) = quad(@(x) x.^n.*exp(a*x-b*x.^2),-1,1,1e-16);
%     end
% end
    % more asymptotics to avoid bad erf values
    cutoff2 = 4.5;
    if ((a-2*b)/2/sqrt(b)<cutoff2)
        Ivals(1) = sqrt(pi/4/b)*exp(a^2/4/b+b-a)*(erf((a+2*b)/(2*sqrt(b))) - erf((a-2*b)/(2*sqrt(b))));
    else
        % use large a asymptotics
        erfas = @(x) -1/sqrt(pi)*(1/x -1/2/x^3+3/4/x^5-15/8/x^7+105/16/x^9-945/32/x^11 ...
            +10395/64/x^13-135135/128/x^15+2027025/256/x^17+34459425/512/x^19-654729075/1024/x^21 ...
            - 13749310575/2048/x^23);
        Ivals(1) = sqrt(pi/4/b)*(exp(-2*a)*erfas((a+2*b)/2/sqrt(b)) - erfas((a-2*b)/2/sqrt(b)));
    end
    %Ivals(1) = sqrt(pi/4/b)*exp(a^2/4/b+b-a)*(erf((a+2*b)/(2*sqrt(b))) - erf((a-2*b)/(2*sqrt(b))));
    Ivals(2) = -1/2/b *(1-exp(-2*a)-a*Ivals(1));
    %
    for n = 0:LMAX-2
        Ivals(n+3) = (1/(n+1)*(1-(-1)^(n+1)*exp(-2*a)) ...
            - a/(n+1)*Ivals(n+2) - Ivals(n+1))*(n+1)/(-2*b);
    end
    %%
end

%%
% rescale by 0 index value
Ivals = Ivals/Ivals(1);
%%
% polynomial forms for the legendre functions
% coefficients listed from x^0 to x^n powers
%Ppoly(i,j) has the coefficient of the x^j term in the i-th polynomial
Ppoly = zeros(LMAX+1,LMAX+1);
Ppoly(1,1) = 1;
Ppoly(2,2) = 1;
for n = 1:LMAX-1
    Ppoly(n+2,:) = ((2*n+1)*[0,Ppoly(n+1,1:end-1)]-n*Ppoly(n,:))/(n+1);
end

xivals = ones(LMAX,1);
% build up legendre projections from power function projections
for n = 1:LMAX
    xivals(n+1) = Ppoly(n+1,:)*Ivals;
end
