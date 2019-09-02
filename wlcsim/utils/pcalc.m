function [X,P]=pcalc(r,NP,L0)
%% calculate end-to-end distribution of polymers

%% initialize
NB=length(r);
N=NB/NP;
L=(N-1)*L0;
LTOT=L;

% setup histograms
XMAX=1.0;
NH=100;
DX=XMAX/NH;
X=transpose((DX/2):DX:(XMAX-DX/2));

IL=round(L/LTOT*N);
L=(IL-1)/(N-1)*LTOT;
I1=((1:1:NP)-1)*N+1;
I2=((1:1:NP)-1)*N+IL;

%%%% start calculations %%%%
P=zeros(NH,1);
REND=sqrt(sum(power(r(I2,1:3)-r(I1,1:3),2),2))/L;
P=P+transpose(hist(REND,X));

P=P/DX/NP;
end