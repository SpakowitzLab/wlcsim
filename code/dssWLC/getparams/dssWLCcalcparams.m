function [eb,gam,epar,eperp,eta,alpha,zetau,delt] = dssWLCcalcparams(del,Nseg,varargin)
% calculate the appropriate parameters for a dssWLC simulation 
% for discretizing a continuous WLC with persistence length lp=1
% usage:
% dssWLCcalcparams(del,Nseg,'interpfile',filename,'alpha',alpha,'LMAX',LMAX
% ,'YcoupleFile',filename)
%
% Required arguments:
% del = segment length
% Nseg = total number of segments in the chain
%
% optional arguments:
%
% interpfile: if supplied then interpolate the energetic parameters from the
% the given space-delimited file. Otherwise, recalculate from scratch
%
% alpha: if single value is supplied then calculate the energetic
% parameters for the given alpha value. If pair of values is supplied then
% this gives the range of possible alpha values, and the code optimizes
% over the possible alphas within this range to get the best length scale
% of accuracy. Default: [1e-5,1.5]
%
% LMAX: if optimizing over alpha values, the maximum l-index of the
% matrices used to calculate structure factors; default=10
%
% YcoupleFile: if optimizing over alpha values, the file containing saved
% coupling coefficients for spherical harmonics. Default: YcoupleSave.mat
% if this file does not exist it will be created as the coupling
% coefficients are calculated.

p = inputParser;
addRequired(p,'del',@isnumeric);
addRequired(p,'Nseg',@(x) isnumeric(x) && x>0);
addParamValue(p,'interpfile',NaN,@ischar);
addParamValue(p,'alpha',[1e-5,2],@isnumeric);
addParamValue(p,'LMAX',10,@(x) (isnumeric(x) && x <= 20 && x>=0));
addParamValue(p,'YcoupleFile','YcoupleSave.mat',@ischar);
parse(p,del,Nseg,varargin{:});

if any(ismember('interpfile',p.UsingDefaults))
    % calculate parameter files directly
    if (any(ismember('alpha',p.UsingDefaults)) || length(p.Results.alpha)>1)
        disp('Calculating dssWLC energetic parameters, optimizing alpha for lowest length scale of accuracy')
        
        % load in the spherical harmonic coupling coefficients
        if (exist(p.Results.YcoupleFile,'file'))
            disp(sprintf('Loading Y coupling coefficients from %s', p.Results.YcoupleFile));
            load(p.Results.YcoupleFile,'Ycouple');
        else
            % recalculate & save if file does not exist
            disp(sprintf('Recalculating spherical harmonic coupling coefficients. Will save in %s and avoid recalculating in future.',p.Results.YcoupleFile))
            LMAX = 20;
            Ycouple = getYcouple(LMAX,p.Results.YcoupleFile);
        end
        
        nk = 100; Ltot = 1000; LMAX = p.Results.LMAX;
        % cutoff in structure factor error for calculating length scale
        cutoff = 1e-4;
        
        klist = logspace(log10(1/del/10),log10(1/del*2),nk);
        
        % get plain wlc structure factor
        display('Calculating plain WLC structure factor')
        nseg = Ltot/del;
        I = eye(LMAX+1);
        for kc = 1:nk
            Mplain = shearWLCpropagator(klist(kc),del,1,1,0,0,0,LMAX);
            Mtot = ((nseg+1)*I+Mplain^(nseg+2) - Mplain*(nseg+2))*inv(Mplain-I)^2;
            Splain(kc) = 2*Mtot(1,1)/(nseg+1)^2;
        end
        
        options = optimset('Display','iter','TolX',1e-4);
        alpharange = p.Results.alpha;
        display('Optimizing over alpha')
        [alpha,lscale] = fminbnd(@(alpha) dssWLClengthScale(del,alpha,Ltot,klist,Splain,Ycouple,cutoff,LMAX),...
            alpharange(1),alpharange(2),options);
        [eb,gam,epari,eperpi,eta,err,plen] = dssWLCminLpParams(del,alpha);
        epar = 1/epari;
        eperp = 1/eperpi;
    else
        alpha = p.Results.alpha
        disp(sprintf('Calculating energetic parameters using alpha=%f',alpha))
        [eb,gam,epari,eperpi,eta,err,plen] = dssWLCminLpParams(del,alpha);
        epar = 1/epari;
        eperp = 1/eperpi;
    end
else    
    % interpolate parameter values from the data file
    disp(sprintf('Interpolating energetic parameters from tabulated file: %s',p.Results.interpfile))
    data = dlmread('dssWLCparams.txt');
    params = interp1(data(:,1),data(:,2:end),del);
    eb = params(1); gam = params(2); epar = params(3); eperp = params(4); eta = params(5);
    alpha = eta^2*eb/eperp;
end

disp(sprintf('Energetic parameters (eb,g,epar,eperp,eta): %f %f %f %f %f',eb,gam,epar,eperp,eta))
disp('Calculating dynamic parameters: zeta_u and del_t')

% find appropriate xiu value
eperph = eperp + eta^2*eb;
xir = 1;
xiulist = logspace(-7,1,50);
L=del;
pval = 1;
tfast = zeros(size(xiulist));

for uc = 1:length(xiulist)
    xiu = xiulist(uc);
    [evals,evecs,pareval] = ssWLCdynamics(eb,gam,epar,eperph,eta,L,xir,xiu,50);
    tfast(uc) = -1/evals(end);
end

dt = diff(log10(tfast));
lxiu = interp1(dt,log10(xiulist(1:end-1)),(dt(1)+dt(end))/2);
zetau = 10^lxiu*del/(1+1/Nseg);

% get the appropriate delt
delt = 0.5/(eperp*gam^2*del)*zetau;

