function [lscale,params,Svals] = dssWLClengthScale(del,alpha,Ltot,klist,Splain,Ycouple,cutoff,LMAX,params)
% get the "length-scale of accuracy" for the discrete shearable WLC
% based on cutoff for structure factor error
% Splain is the structure factor used for comparison (must be same size as
% klist)
% alpha =eta^2*eb/eperp is the dimensionless parameter setting coupling
% strength
% if params is supplied use the given parameters
% (eb,gamma,epari,eperpi,eta)
% otherwise find the parameters to match the plain WLC for a given alpha

X = alpha/(1+alpha);

if (nargin<8)
    LMAX=10;
end
    
if (nargin<9)
    % get the WLC parameters
    [eb,g,epari,eperpi,eta,err,plen] = dssWLCminLpParams(del,alpha);
    params = [eb,g,epari,eperpi,eta,plen,err];
else
    eb = params(1); g = params(2); epari = params(3); eperpi = params(4); eta = params(5);
end

Svals = dssWLCstructFact(klist,del,eb,g,epari,eperpi,eta,Ycouple,Ltot,LMAX);

errvals = abs(Svals-Splain)./Splain;
indv = find(errvals>cutoff/100);
if (length(indv)<2)
    error(sprintf('Values out of range with alpha = %f', alpha))
end
if (errvals(indv(1))>cutoff)
    display(sprintf('Cannot get down to the cutoff with alpha=%f. Try extending to lower k', alpha))
    lscale=inf;
    params = []; Svals = [];
    return
end

errvals = real(errvals);
% find local maxima
[pks,locs] = findpeaks(errvals(indv));
% find the first max that's above cutoff
[a,b] = find(errvals(indv(locs))>cutoff);
if (length(b)<1)
    % no local max above cutoff; search from the end
    if (max(errvals)<cutoff) 
         error(sprintf('No values above cutoff, try increasing maximum k. alpha=%f',alpha))
    else
        kshear = fzero(@(x) interp1(klist(indv),real(errvals(indv)),x,'linear','extrap')-cutoff,klist(indv(end)));
    end
else
    % search below first local max above cutoff
    indv = indv(1:locs(b(1)));
     kshear = fzero(@(x) interp1(klist(indv),real(errvals(indv)),x,'linear','extrap')-cutoff,klist(indv(end)));
end
lscale = 1/kshear;
end