function [g,epari,eperpi,eta]= dssWLCgetParams(del,eb,alpha,R2l,R2c,R4l)
% -------
% calculate g,epari,eperpi for a dssWLC model
% where del, alpha, eb are fixed and we want to fit the moments R2l, R2c, R4l
% alpha = eta^2*eb/eperp

% solve for eb to match lpeff (using xi_1)
% solve for eperpih = 1/(eperp-hat) as a function of g to fit R2c
% solve for epari as a function of g to fit R2l
% solve for g to fit R4l
% -------

if max(abs(imag([del,eb,alpha])))>0 || eb<0 || del <0 || alpha<0 
    moms = 0; M1=0; epari = 0; eperpi=0; R2ctmp = 0; g =0;
    return
end

x = alpha/(1+alpha)*eb;
M0=1;
%%
LMAX=4;
[xivals,tmp] = expandFsph(del,eb,alpha,LMAX);

M1 = xivals(2); M2 = xivals(3); M3 = xivals(4); M4 = xivals(5);

%[M1 M2 M3 M4]
%return
% --------------
%% coupling coefficients
% go up one with D10
alj = @(l,j) sqrt(3*(l-j)*(l+j)/((2*l+1)*(2*l-1)));
% go up one in both l and j with D11
ahlj = @(l,j) sqrt(3/2)*sqrt((l+j)*(l+j-1)/((2*l-1)*(2*l+1)));
% stay in place with D20
blj = @(l,j) sqrt(5)*(l^2+l-3*j^2)/((2*l-1)*(2*l+3));
% go up 2 with D20, for j=0
bhlj = @(l,j) 3*sqrt(5)/2/(2*l-1)*sqrt((j-l)*(j-l+1)*(j+l)*(j+l-1)/((2*l+1)*(2*l-3)));
% go up one with D30
clj = @(l,j) -sqrt(7)*3/2*(1+5*j^2-l^2)/((2*l-3)*(2*l+3))*alj(l,j)/sqrt(3);
% stay put with D40
dlj = @(l,j) 9/4*(35*j^4 + 3*(l-1)*l*(l+1)*(l+2) - 5*j^2*(6*l^2+6*l-5))...
    /((2*l-3)*(2*l-1)*(2*l+3)*(2*l+5));
dhlj = @(l,j) 15/2*(l^2-l-7*j^2-2)/((2*l-5)*(2*l-1)*(2*l+3))*sqrt((j-l)*(j-l+1)*(j+l)*(j+l-1)/((2*l+1)*(2*l-3)));

%% get ec/eperp as function of gamma
A = 2/(M1-1)^2*i*del/sqrt(3);
B = -i/3/sqrt(3)*(2 + 2*sqrt(5)*M2*blj(1,1));
C = i*del/sqrt(3)*M1;
Z1 = R2c/A/B; Z2 = -C/B;

%% get epari as a function of gamma

%s00 = M0*(  - 2/3*(del/2*epari + del*eperpi)) ...
%    +2*ec^2/9*eperpi^2*sqrt(5)*M2*bhlj(2,0)*ahlj(2,0)*ahlj(1,1);


A = (M1-1)*R2c/del;
B = -1/del*(-2/9*ahlj(1,1)^2+2*sqrt(5)/9*M2*bhlj(2,0)*ahlj(2,0)*ahlj(1,1)-2*del/3/x);
C = del*alj(1,0)/3;
D = 1/3;
Z3 = -B*Z1^2/D;
Z4 = (-B*Z2^2-C)/D;
Z5 = (R2l - A - 2*B*Z1*Z2)/D;

%E = Z1/g + Z2*g;
%eperpi = E^2/x;
%epari = Z3/g^2 + Z4*g^2 + Z5;

%% to fit R4l, go through each of the individual steps and calculate
% coefficients of: 1/g^4, 1/g^3,1/g^2, 1/g, 1, g, g^2, g^3, g^4
Eg =  [0,0,0,Z1,0,Z2,0,0,0]; %E
Eg2 = [0,0,Z1^2,0,2*Z1*Z2,0,Z2^2,0,0]; %E^2
Eg3 = [0,Z1^3,0,3*Z1^2*Z2,0,3*Z1*Z2^2,0,Z2^3,0];%E^3
Eg4 = [Z1^4,0,4*Z1^3*Z2,0,6*Z1^2*Z2^2,0,4*Z1*Z2^3,0,Z2^4];%E^4
Eparg = [0,0,Z3,0,Z5,0,Z4,0,0]; %epari
EpargG = [0,0,0,Z3,0,Z5,0,Z4,0]; %epari*g
EpargG2 = [0,0,0,0,Z3,0,Z5,0,Z4]; %epari*g^2
EpargEg = [0,Z1*Z3, 0,Z2*Z3+Z1*Z5,0,Z2*Z5+Z1*Z4,0,Z2*Z4,0]; %E*epari
EpargEg2 = [Z1^2*Z3,0,Z1^2*Z5+2*Z1*Z2*Z3,0,Z2^2*Z3+Z1^2*Z4+2*Z1*Z2*Z5,0,Z2^2*Z5+2*Z1*Z2*Z4,0,Z2^2*Z4];
Eparg2 = [Z3^2,0,2*Z3*Z5,0,Z5^2+2*Z3*Z4,0,2*Z4*Z5,0,Z4^2];%epari^2
G = [0,0,0,0,0,1,0,0,0]; %g
G2 = [0,0,0,0,0,0,1,0,0]; %g^2
G3 = [0,0,0,0,0,0,0,1,0]; %g^3
G4 = [0,0,0,0,0,0,0,0,1]; %g^4
EgG = [0,0,0,0,Z1,0,Z2,0,0]; %E*g
Eg2G = [0,0,0,Z1^2,0,2*Z1*Z2,0,Z2^2,0]; %E^2*g
EgG2 = [0,0,0,0,0,Z1,0,Z2,0]; %E*g^2
Eg2G2 = [0,0,0,0,Z1^2,0,2*Z1*Z2,0,Z2^2]; %E^2*g^2

% single derivatives

% get s01 (step up from 0 to 1, X10 in the notes)
s01= -i/sqrt(27)*Eg*(2*M0 + 2*sqrt(5)*M2*blj(1,1)) + i*del/sqrt(3)*M1*G;

% get s10 (step down from 1 to 0, X01 in the notes)
s10 = i*del*G/sqrt(3)*M0;

% s12 (step up from 1 to 2, X21 in the notes)
s12 =  -i*Eg/(3*sqrt(5))*(2*sqrt(3)*alj(2,1)^2*M1 + 2*sqrt(7)*clj(2,1)*alj(2,1)*M3) ...
    + i*del*G/sqrt(3)*alj(2,0)*bhlj(2,0)*M2;

%s21 (step down from 2 to 1, X12 in notes)
s21 = -2*i*Eg/(3*sqrt(3))*(blj(1,1)*M0 + blj(1,1)^2*sqrt(5)*M2)...
    + i*del*G/sqrt(3)*blj(1,0)*M1;

%-------------
% double derivatives
s00 = M0*(-del^2*G2/3*alj(1,0) - 2*Eg2/9*ahlj(1,1)^2 - 2/3*(del/2*Eparg + del*Eg2/x)) ...
    +2*Eg2/9*sqrt(5)*M2*bhlj(2,0)*ahlj(2,0)*ahlj(1,1);

s11 = 4*del*EgG/9*alj(2,1)^2*(M0 + blj(1,1)*sqrt(5)*M2) ...
    - 2*ahlj(2,0)^2*Eg2/9*alj(2,0)*(alj(2,0)*M1 + clj(2,0)*M3*sqrt(7/3)) ...
    + (-del^2*G2/3*alj(2,0) + Eg2/9*2*ahlj(1,1)*ahlj(2,0) - 2/3/sqrt(5)*del*(Eparg-Eg2/x))*alj(2,0)*M1 ...
    + 2*Eg2/9*ahlj(2,0)*ahlj(1,1)*(alj(2,0)*M1+clj(2,0)*M3*sqrt(7/3))...
    +(-del^2*G2/3*alj(1,0) - 2*ahlj(1,1)^2*Eg2/9 - 2*del/3*(Eparg/2+Eg2/x))*M1;


s02 = -2/9*Eg2*ahlj(2,2)^2 *(dlj(2,2)*M4*sqrt(9/5) + blj(2,2)*M2+M0/sqrt(5)) ...
    - 2/9*Eg2*ahlj(2,0)^2*(dlj(2,0)*M4*sqrt(9/5)+blj(2,0)*M2+M0/sqrt(5)) ...
    +4/3/sqrt(15)*del*EgG*alj(2,1)*(alj(2,1)*M1*sqrt(3) + clj(2,1)*M3*sqrt(7)) ...
    +(-del^2*G2*alj(2,0)/3 + 2*Eg2/9*ahlj(1,1)*ahlj(2,0) - 2*del/3/sqrt(5)*(Eparg-Eg2/x))*M2;

s20 = (-del^2*G2*alj(2,0)/3 + 2*Eg2/9*ahlj(1,1)*ahlj(2,0) - 2*del/3/sqrt(5)*(Eparg-Eg2/x))*bhlj(2,0)*M0 ...
    - 2*Eg2/9*ahlj(2,0)^2*bhlj(2,0)^2*M2*sqrt(5);

% -------------
% triple derivatives

v0010 = -i*del^3*G3/3/sqrt(3)*(alj(2,0)^2+alj(1,0)^2) - i*6*del*Eg2G/9/sqrt(3)*(-ahlj(2,0)*alj(2,0)*ahlj(1,1)+ahlj(1,1)^2*alj(1,0)) ...
    - 6*i*del^2/3/sqrt(3)*(EpargG/2+Eg2G/x) - 6*i*del^2/3/sqrt(15)*(EpargG-Eg2G/x) *alj(2,0);

v2010 = -i*6*del*Eg2G/9/sqrt(3)*(ahlj(2,0)^2*alj(2,0) - ahlj(1,1)*alj(1,0)*ahlj(2,0));

v1111 = 3*i*Eg3/27*(-ahlj(2,0)^2-ahlj(1,1))^2+3*i*del^2*EgG2/9*(alj(1,0)-alj(2,0)*ahlj(2,0)) ...
    +6*i*del/9*(EpargEg/2+Eg3/x)-6*i*del/9/sqrt(5)*(EpargEg-Eg3/x)*ahlj(2,0);
v3111 = 3*i*Eg3/27*(-ahlj(2,0)^2-ahlj(1,1))*ahlj(2,0)*ahlj(3,1);

%MRT3(ljind(3,-1),ljind(1,-1))/(8*pi^2)

t01 = 2/sqrt(3)*v1111*(M0+blj(1,1)*M2*sqrt(5)) + 2/sqrt(3)*v3111*(bhlj(3,1)*M2*sqrt(5)+dhlj(3,1)*M4*sqrt(9)) ...
    +v2010*(alj(2,0)*M1+bhlj(3,0)*M3*sqrt(7/3)) + v0010*M1;

t10 = M0*v0010 + bhlj(2,0)*M2*sqrt(5)*v2010;

% -----------------------
% quadruple derivatives
v0000 = 6*Eg4/81*(ahlj(2,2)^2*ahlj(1,1))^2 ...
    + del^4*G4/9*(alj(2,0)^2*alj(1,0)+alj(1,0)^3) ...
    + 12*del^2*Eg2G2/27*(-alj(1,0)*ahlj(1,1)^2 + alj(2,0)*ahlj(2,0)*ahlj(1,1))*-ahlj(1,1)...
    + 12*del^2/9*(Eparg2/4+Eg4/x^2+EpargEg2/x) + 12*del^2/45*(Eparg2+Eg4/x^2-2*EpargEg2/x)*bhlj(2,0) ...
    + 12*del/3*(EpargG2/2+Eg2G2/x)*del^2/3*alj(1,0)...
    +12*del/3*(EpargEg2/2+Eg4/x)*2/9*ahlj(1,1)^2 ...
    + 12*del/3/sqrt(5)*(EpargG2-Eg2G2/x)*del^2/3*alj(2,0)*bhlj(2,0) ...
    -12*del/3/sqrt(5)*(EpargEg2-Eg4/x)*2/9*ahlj(2,0)*ahlj(1,1)*bhlj(2,0);
%gvec = [1/g^4;1/g^3;1/g^2;1/g;1;g;g^2;g^3;g^4];


v2000 = 6*Eg4/81 *ahlj(2,2)^2*ahlj(1,1)*(-ahlj(2,2)^2*ahlj(2,0) - ahlj(2,2)*ahlj(3,-1)*ahlj(3,1)) ...
 + 12*del^2*Eg2G2/27*(-alj(1,0)*ahlj(1,1)^2 + alj(2,0)*ahlj(2,0)*ahlj(1,1))*ahlj(2,0) ...
 -24*del/27*(EpargEg2/2+Eg4/x)*ahlj(2,0)*ahlj(1,1) + 24*del/27/sqrt(5)*(EpargEg2-Eg4/x)*ahlj(2,0)^2*bhlj(2,0);

v4000 = 6*Eg4/81*ahlj(2,2)^2*ahlj(1,1)*ahlj(2,2)*ahlj(3,-1)*ahlj(4,0);

q00 = M0*v0000 + M2*sqrt(5)*v2000 + M4*sqrt(9)*v4000;

X1 = -1/2*(M1-5)/(M1-1)^3*conv(conv(s01,s01),conv(s10,s10)); X1 = X1(13:21);
X2 = 1/2*(M1-3)/(M1-1)^2*conv(conv(s00,s01),s10); X2 = X2(9:17);
X3 = -1/2*conv(s00,s00); X3 = X3(5:13);
X4 = -1/((M1-1)^2*(M2-1))*conv(conv(s01,s12),conv(s21,s10)); X4 = X4(13:21);
X5 = 1/(M1-1)^2*conv(conv(s01,s11),s10); X5 = X5(9:17);
X6 = 1/((M1-1)*(M2-1))*conv(conv(s02,s21),s10); X6 = X6(9:17);
X7 = 1/((M1-1)*(M2-1))*conv(conv(s01,s12),s20); X7 = X7(9:17);
X8 = -1/(M2-1)*conv(s02,s20); X8 = X8(5:13);
X9 = 1/(1-M1)*conv(t01,s10); X9 = X9(5:13);
X10 = 1/(1-M1)*conv(s01,t10); X10 = X10(5:13);
X11 = q00;

% overall polynomial for R4l
R4lpoly = 1/del*(24*X1+24*X2+6*X3+24*X4+12*X5+12*X6+12*X7+6*X8+4*X9+4*X10+X11);
tmppoly = R4lpoly;

% solve for gamma
R4lpoly(5) = R4lpoly(5)-R4l;

groots = roots(fliplr(R4lpoly));

% want only the positive real solutions
gvals = groots(find(groots>0 & imag(groots)==0));

if (length(gvals) == 0)
    %disp('Warning: No positive real g solution')
    g = NaN; eperpi = NaN; epari = NaN; ec = NaN;
else
    g = gvals;
    E = Z1./gvals + Z2*gvals;
    eperpi = E.^2/x;
    ec = E./eperpi;
    epari = Z3./gvals.^2 + Z4*gvals.^2 + Z5;
end

%gvec = [1/g^4;1/g^3;1/g^2;1/g;1;g;g^2;g^3;g^4];

%s11*gvec
%[tmppoly*gvec,R4l]
%X11*gvec
%t01=t01*gvec
%t10 = t10*gvec
%q00 = q00*gvec

eperpi = 1./(1./eperpi-ec.^2./eb);
eta = -ec/eb;

end