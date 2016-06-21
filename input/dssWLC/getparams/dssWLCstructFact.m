% stuff for testing only
% load('lengthscales.mat')
% dc = 15; del = dellist(dc);
% [lp,g,epari,eperpi,ec,err,plen] = discreteWLCshearParams(del,Xminvals(dc));
% eperpi = 1/(1/eperpi-ec^2/lp);
% eb = lp;
% eta = -ec/eb;
function [Svals,LMAX] = dssWLCstructFact(klist,del,eb,g,epari,eperpi,eta,Ycouple,Ltot,LMAXin)

%%
% ---------
% Find the structure factor for a dssWLC model
% ----------


pmax = 10; % maximum power to which exponential is expanded
nk = length(klist);

eperpih = 1/(1/eperpi + eta^2*eb);
alpha = eta^2*eb*eperpi;

%% 
% expand F in terms of spherical harmonics
if (eperpi==0)
    % non-shearable WLC, do direct projection using plane wave decomposition
    for c = 0:LMAXin
     xivals(c+1) = besseli(c+0.5,eb/del)/besseli(0.5,eb/del);
    end
else
    xivals = expandFsph(del,eb,alpha,LMAXin);
end

LMAX=LMAXin;
% cut off LMAX such that xi values are decreasing
% tmp = diff(abs(xivals)); 
% ind = find(tmp>0);
% if (length(ind)>0)
%     LMAX = ind(1)-1;
%     xivals = xivals(1:LMAX+1);
% else
%     LMAX=LMAXin;
% end
%%
Imax = ljind(LMAX,LMAX);
% ----
%% get the coefficients h_a,b^j for the spherical harmonic expansion of
% exp(H) just for 1 specific k value, in order to find the nonzero elements
% in the sparse matrices
%Hmat(a,b) has the coefficient for Y_a*(u0hat) Y_b(uhat)
% basically we're working in the linear space of pairs of Y functions

k = klist(end);

% get the coefficients for pairs of Y in the H matrix 
Hmat = sparse(zeros(Imax,Imax));
Hmat(ljind(1,0),1) = i*k*del*g*4*pi/sqrt(3);
Hmat(ljind(2,0),1) =-k^2*del/3*4*pi/sqrt(5)*(epari - eperpih);
Hmat(1,1) = -k^2*del/2*eperpih*4*pi + Hmat(ljind(2,0),1)*sqrt(5)/2;
Hmat(ljind(1,1),ljind(1,1)) = eta*eb*eperpih*i*k*4*pi/3;
Hmat(ljind(1,-1),ljind(1,-1)) = eta*eb*eperpih*i*k*4*pi/3;

% get the spherical harmonic expansion coefficients for exp(H) by expanding the exponential as a
% power series in the Ys (accurate for small k)

% Hpow contains the sph harmonic expansion of H to
% different powers
Hpow = {};
% 0th power
Hpow{1} = sparse(zeros(Imax,Imax)); 
Hpow{1}(1,1) = 4*pi;

% 1st power
Hpow{2} = Hmat;        
    
[a1list,b1list,vals1] = find(Hmat);
nv1= length(a1list);

% coefficients for the exp(H) matrix
expHmat = Hpow{1};
if (pmax>=2)
   expHmat = expHmat + Hpow{2};
end

for p = 2:pmax
    Hpow{p+1} = sparse(zeros(Imax,Imax));
    [a2list,b2list,vals2] = find(Hpow{p});
    nv2 = length(a2list);
    for c1 = 1:nv1
        aind1 = a1list(c1);
        bind1 = b1list(c1);
        for c2 = 1:nv2
            aind2 = a2list(c2);
            bind2 = b2list(c2);
            % product of 2 pairs of sph harmonics, expand into a asum over
            % single pairs
            Hpow{p+1} = Hpow{p+1} + vals1(c1)*vals2(c2)*transpose(Ycouple{aind1}(aind2,1:Imax))*Ycouple{bind1}(bind2,1:Imax);
        end
    end
    
    expHmat = expHmat + 1/factorial(p)*Hpow{p+1};
    
end


%% get M tensor elements, independent of k value (so long as nonzero 
% form of expHmat remains the same)
%Mmat{c}(l,l0) is the value for the c-th pair of a,j and b,j indices

[aind,bind,vals] = find(expHmat);

Mmat={};
for c = 1:length(aind)
    %[c, size(Mmat), length(aind)]
    a = aind(c); b = bind(c);
    [al,aj] =ljindinv(a); [bl,bj] = ljindinv(b);
    
    Mmat{c} = sparse(zeros(LMAX+1,LMAX+1));
    
    for l0 = 0:LMAX
        [lind,tmp,al0vals] = find(Ycouple{ljind(l0,0)}(1:Imax,a));
        
        for lc = 1:length(lind)
            [ll,tmp] = ljindinv(lind(lc));
            
            for lam = max(0,abs(bl-ll)):min(LMAX,bl+ll)
                tmp = Ycouple{b}(ljind(lam,0),ljind(ll,bj));
                if (abs(tmp)>0)
                    Mmat{c}(ll+1,l0+1) = Mmat{c}(ll+1,l0+1) + tmp*al0vals(lc)...
                        *xivals(lam+1)*sqrt((2*lam+1)/(2*ll+1));
                end
            end
        end
    end
end

%% get overall structure factor, recalculating h_a,b^j for each specific k
% value
for kc = 1:nk
    
    k = klist(kc);
    
    % get the coefficients for pairs of Y in the H matrix
    Hmat = sparse(zeros(Imax,Imax));
    Hmat(ljind(1,0),1) = i*k*del*g*4*pi/sqrt(3);
    Hmat(ljind(2,0),1) =-k^2*del/3*4*pi/sqrt(5)*(epari - eperpih);
    Hmat(1,1) = -k^2*del/2*eperpih*4*pi + Hmat(ljind(2,0),1)*sqrt(5)/2;
    Hmat(ljind(1,1),ljind(1,1)) = eta*eb*eperpih*i*k*4*pi/3;
    Hmat(ljind(1,-1),ljind(1,-1)) = eta*eb*eperpih*i*k*4*pi/3;
    
    % get the spherical harmonic expansion coefficients for exp(H) by expanding the exponential as a
    % power series in the Ys (accurate for small k)
    
    % Hpow contains the sph harmonic expansion of H to
    % different powers
    Hpow = {};
    % 0th power
    Hpow{1} = sparse(zeros(Imax,Imax));
    Hpow{1}(1,1) = 4*pi;
    
    % 1st power
    Hpow{2} = Hmat;
    
    [a1list,b1list,vals1] = find(Hmat);
    nv1= length(a1list);
    
    % coefficients for the exp(H) matrix
    expHmat = Hpow{1};
    if (pmax>=2)
        expHmat = expHmat + Hpow{2};
    end
    
    for p = 2:pmax
        Hpow{p+1} = sparse(zeros(Imax,Imax));
        [a2list,b2list,vals2] = find(Hpow{p});
        nv2 = length(a2list);
        for c1 = 1:nv1
            aind1 = a1list(c1);
            bind1 = b1list(c1);
            for c2 = 1:nv2
                aind2 = a2list(c2);
                bind2 = b2list(c2);
                % product of 2 pairs of sph harmonics, expand into a asum over
                % single pairs
                Hpow{p+1} = Hpow{p+1} + vals1(c1)*vals2(c2)*transpose(Ycouple{aind1}(aind2,1:Imax))*Ycouple{bind1}(bind2,1:Imax);
                %Hpow{p+1} = Hpow{p+1} + vals1(c1)*vals2(c2)*transpose(Ycoupleflip{aind1}(1:Imax,aind2))*Ycoupleflip{bind1}(1:Imax,bind2);
            end
        end
        
        expHmat = expHmat + 1/factorial(p)*Hpow{p+1};        
    end
    
    % overall matrix of g coefficients
    [aind,bind,vals] = find(expHmat);
    gcoeff = zeros(LMAX+1,LMAX+1);
    
    for c = 1:length(aind)
        gcoeff = gcoeff + Mmat{c}*vals(c);
    end
    
    % number of segments    
    nseg = Ltot/del;
        
    % now find the structure factor
    I = eye(LMAX+1);
    gtot = ((nseg+1)*I+gcoeff^(nseg+2) - gcoeff*(nseg+2))*inv(gcoeff-I)^2;        
    Svals(kc) = 2*gtot(1,1)/(nseg+1)^2;   
    
end