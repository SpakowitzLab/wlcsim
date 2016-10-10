function [ebeff,geff,epari,eperpi,etaeff,err,plen] = dssWLCminLpParams(del,alpha,Rcoeff)
% find the parameters for the effective shearable WLC
% for a given discretization
% alpha is fixed; push persistence length as close as possible to 1

if (nargin>2)
    R2l = Rcoeff(1); R2c = Rcoeff(2); R4l = Rcoeff(3); R4c = Rcoeff(4);
else
    % these are the moment components for a continuous plain worm-like chain
    R2l = 2/3;
    R2c = -2/3;
    R4l = -208/45;
    R4c = 856/135;
end

if (nargin<2)
    X =0.5;
end

% plen is the persistence length; we want to find the lowest persistence
% length for which a solution exists

tol = 1e-6; 
plenmin = 1;
maxtry = 100;
options = optimset('TolX',1e-8,'MaxFunEvals',1e3,'MaxIter',5e3,'TolFun',1e-8);
ebguess = 1.2;

% first go through linearly to find an upper bound 
% (lowest p where there is a solution
maxtry1 = 100;
stp = logspace(-6,log10(max(3,del)),maxtry1);
plen = plenmin;
M1func = @(eb) feval(@(x) x(2),expandFsph(del,eb,alpha,1));
for tryc = 1:maxtry1
    plen = plenmin + stp(tryc);
    
    M1 = exp(-del/(plen));     
    try
        [ebeff,fval,errval] = fzero(@(eb) M1func(eb) - M1, ebguess,options);
    catch err
        disp(sprintf('Failed to solve for eb, moving to next plen. %f', plen))
        continue
    end
    [geff,epari,eperpi,etaeff]= dssWLCgetParams(del,ebeff,alpha,R2l,R2c,R4l);
    %[plen, ebeff]
    %[geff epari eperpi etaeff]
    ind = find(epari>0 & eperpi>0 & etaeff > 0 & geff<1);
    
    if (length(ind)>0)
        plenmax = plen;
        if (tryc>1)
            plenmin = plenmin+stp(tryc-1);
        end
        break
    end
end

plen = (plenmin+plenmax)/2;
for tryc = 1:maxtry
    M1 = exp(-del/(plen));
    [ebeff,fval,errval] = fzero(@(eb) M1func(eb) - M1, ebguess,options);
    [geff,epari,eperpi,etaeff]= dssWLCgetParams(del,ebeff,alpha,R2l,R2c,R4l);
    ind = find(epari>0 & eperpi>0 & etaeff > 0);    
    
    if (length(ind) > 0)
        plenmax = plen;        
    else
        plenmin = plen;
    end
    plen = (plenmin+plenmax)/2;
    
    if (plenmax-plenmin) < tol
        break
    end
    
    %[tryc, plenmin, plenmax, plen, length(ind)]
end
if (tryc>=maxtry); disp('WARNING: failed to find minimal plen'); end

plen = plenmax;
M1 = exp(-del/(plen));

% now find the parameters given this persistence length
[ebeff,fval,errval] = fzero(@(eb) M1func(eb) - M1, ebguess,options);
[geff,epari,eperpi,etaeff]= dssWLCgetParams(del,ebeff,alpha,R2l,R2c,R4l);
ind = find(epari>0 & eperpi>0 & etaeff > 0);
geff = geff(ind); epari = epari(ind); eperpi = eperpi(ind); etaeff = etaeff(ind);
% calculate the error in the constant term of R^4, pick the solution with
% the smallest error
errcomp = zeros(size(ind));
for c3 = 1:length(ind)
    moms = dssWLCmoments(del,ebeff,geff(c3),epari(c3),eperpi(c3),etaeff(c3),1);
    errcomp(c3) = (moms(end)-R4c)/R4c;
end
[a,b] = min(abs(errcomp));

geff = geff(b); epari = epari(b); eperpi = eperpi(b); etaeff = etaeff(b);
err = errcomp(b);
end