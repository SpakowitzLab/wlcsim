function moms = dssWLCmoments(del,eb,g,epari,eperpi,eta,getR4c)
% --------------------------
% get the loworder moment components of the dssWLC with the given
% parameters
% R2l, R2c are linear and constant components of <R^2>, repectively
% R4l is the linear component of <R^4>
% if getR4c is present and true, also get the constant components of R^4
% ---------------------------

if (nargin<7)
    getR4c = 0;
end

alpha = eta^2*eb*eperpi;
ec=-eta*eb;
% inverse epsilon hat
eperphi = 1/(1/eperpi + eta^2*eb);

LMAX=4;
[xivals,tmp] = expandFsph(del,eb,alpha,LMAX);
M0=1;
M1 = xivals(2); M2 = xivals(3); M3 = xivals(4); M4 = xivals(5);

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
% ----------
% get individual diagram steps

% single derivatives

% get s01 (step up from 0 to 1, X10 in the notes)
s01= -i/sqrt(27)*ec*eperphi*(2*M0 + 2*sqrt(5)*M2*blj(1,1)) + i*del*g/sqrt(3)*M1;

% get s10 (step down from 1 to 0, X01 in the notes)
s10 = i*del*g/sqrt(3)*M0;

% s12 (step up from 1 to 2, X21 in the notes)
s12 =  -i*ec/(3*sqrt(5))*eperphi*(2*sqrt(3)*alj(2,1)^2*M1 + 2*sqrt(7)*clj(2,1)*alj(2,1)*M3) ...
    + i*del*g/sqrt(3)*alj(2,0)*bhlj(2,0)*M2;

tmp1 = i*del*g/sqrt(3)*alj(2,0)*bhlj(2,0)*M2;
%s21 (step down from 2 to 1, X12 in notes)
s21 = -2*i*ec/(3*sqrt(3))*eperphi*(blj(1,1)*M0 + blj(1,1)^2*sqrt(5)*M2)...
    + i*del*g/sqrt(3)*blj(1,0)*M1;
tmp2= i*del*g/sqrt(3)*blj(1,0)*M1;
%-------------
% double derivatives
s00 = M0*(-del^2*g^2/3*alj(1,0) - 2*ec^2/9*eperphi^2*ahlj(1,1)^2 - 2/3*(del/2*epari + del*eperphi)) ...
    +2*ec^2/9*eperphi^2*sqrt(5)*M2*bhlj(2,0)*ahlj(2,0)*ahlj(1,1);

s11 = 4*del*g*ec/9*eperphi*alj(2,1)^2*(M0 + blj(1,1)*sqrt(5)*M2) ...
    - 2*ahlj(2,0)^2*ec^2*eperphi^2/9*alj(2,0)*(alj(2,0)*M1 + clj(2,0)*M3*sqrt(7/3)) ...
    + (-del^2*g^2/3*alj(2,0) + ec^2/9*eperphi^2*2*ahlj(1,1)*ahlj(2,0) - 2/3/sqrt(5)*del*(epari-eperphi))*alj(2,0)*M1 ...
    + 2*ec^2/9*eperphi^2*ahlj(2,0)*ahlj(1,1)*(alj(2,0)*M1+clj(2,0)*M3*sqrt(7/3))...
    +(-del^2*g^2/3*alj(1,0) - 2*ahlj(1,1)^2*ec^2/9*eperphi^2 - 2*del/3*(epari/2+eperphi))*M1;

s02 = -2/9*ec^2*eperphi^2*ahlj(2,2)^2 *(dlj(2,2)*M4*sqrt(9/5) + blj(2,2)*M2+M0/sqrt(5)) ...
    - 2/9*ec^2*eperphi^2*ahlj(2,0)^2*(dlj(2,0)*M4*sqrt(9/5)+blj(2,0)*M2+M0/sqrt(5)) ...
    +4/3/sqrt(15)*del*g*ec*eperphi*alj(2,1)*(alj(2,1)*M1*sqrt(3) + clj(2,1)*M3*sqrt(7)) ...
    +(-del^2*g^2*alj(2,0)/3 + 2*ec^2/9*eperphi^2*ahlj(1,1)*ahlj(2,0) - 2*del/3/sqrt(5)*(epari-eperphi))*M2;

s20 = (-del^2*g^2*alj(2,0)/3 + 2*ec^2/9*eperphi^2*ahlj(1,1)*ahlj(2,0) - 2*del/3/sqrt(5)*(epari-eperphi))*bhlj(2,0)*M0 ...
    - 2*ec^2/9*eperphi^2*ahlj(2,0)^2*bhlj(2,0)^2*M2*sqrt(5);

% -------------
% triple derivatives

v0010 = -i*del^3*g^3/3/sqrt(3)*(alj(2,0)^2+alj(1,0)^2) - i*6*del*g*ec^2/9/sqrt(3)*eperphi^2*(-ahlj(2,0)*alj(2,0)*ahlj(1,1)+ahlj(1,1)^2*alj(1,0)) ...
    - 6*i*del^2*g/3/sqrt(3)*(epari/2+eperphi) - 6*i*del^2*g/3/sqrt(15)*(epari-eperphi) *alj(2,0);
v2010 = -i*6*del*g*ec^2/9/sqrt(3)*eperphi^2*(ahlj(2,0)^2*alj(2,0) - ahlj(1,1)*alj(1,0)*ahlj(2,0));
v1111 = 3*i*ec^3/27*eperphi^3*(-ahlj(2,0)^2-ahlj(1,1))^2+3*i*del^2*g^2*ec/9*eperphi*(alj(1,0)-alj(2,0)*ahlj(2,0)) ...
    +6*i*ec*del/9*eperphi*(epari/2+eperphi)-6*i*ec*del/9/sqrt(5)*eperphi*(epari-eperphi)*ahlj(2,0);
v3111 = 3*i*ec^3/27*eperphi^3*(-ahlj(2,0)^2-ahlj(1,1))*ahlj(2,0)*ahlj(3,1);

%MRT3(ljind(3,-1),ljind(1,-1))/(8*pi^2)

t01 = 2/sqrt(3)*v1111*(M0+blj(1,1)*M2*sqrt(5)) + 2/sqrt(3)*v3111*(bhlj(3,1)*M2*sqrt(5)+dhlj(3,1)*M4*sqrt(9)) ...
    +v2010*(alj(2,0)*M1+bhlj(3,0)*M3*sqrt(7/3)) + v0010*M1;

t10 = M0*v0010 + bhlj(2,0)*M2*sqrt(5)*v2010;

% -----------------------
% quadruple derivatives
v0000 = 6*ec^4/81*eperphi^4*(ahlj(2,2)^2*ahlj(1,1))^2 ...
    + del^4*g^4/9*(alj(2,0)^2*alj(1,0)+alj(1,0)^3) ...
    + 12*del^2*g^2*ec^2/27*eperphi^2*(-alj(1,0)*ahlj(1,1)^2 + alj(2,0)*ahlj(2,0)*ahlj(1,1))*-ahlj(1,1)...
    + 12*del^2/9*(epari/2+eperphi)^2 + 12*del^2/45*(epari-eperphi)^2*bhlj(2,0) ...
    + 12*del/3*(epari/2+eperphi)*(del^2*g^2/3*alj(1,0)+2*ec^2/9*eperphi^2*ahlj(1,1)^2) ...
    + 12*del/3/sqrt(5)*(epari-eperphi)*(del^2*g^2/3*alj(2,0)-2*ec^2/9*eperphi^2*ahlj(2,0)*ahlj(1,1))*bhlj(2,0);

v2000 = 6*ec^4/81*eperphi^4 *ahlj(2,2)^2*ahlj(1,1)*(-ahlj(2,2)^2*ahlj(2,0) - ahlj(2,2)*ahlj(3,-1)*ahlj(3,1)) ...
 + 12*del^2*g^2*ec^2/27*eperphi^2*(-alj(1,0)*ahlj(1,1)^2 + alj(2,0)*ahlj(2,0)*ahlj(1,1))*ahlj(2,0) ...
 -24*del/27*(epari/2+eperphi)*ec^2*eperphi^2*ahlj(2,0)*ahlj(1,1) + 24*del*ec^2/27/sqrt(5)*eperphi^2*(epari-eperphi)*ahlj(2,0)^2*bhlj(2,0);

v4000 = 6*ec^4/81*eperphi^4*ahlj(2,2)^2*ahlj(1,1)*ahlj(2,2)*ahlj(3,-1)*ahlj(4,0);

%MRT4(ljind(4,0),1)/(8*pi^2)

q00 = M0*v0000 + M2*sqrt(5)*v2000 + M4*sqrt(9)*v4000;

% % -------------------------
% %first level diagrams

 R2l = 2*s01*s10/(M1-1)/del - s00/del;
 R2c = 2*s01*s10/(M1-1)^2; 
% 
% ------------------
% second level diagrams

X1 = -1/2*(M1-5)/(M1-1)^3*s01^2*s10^2;
X2 = 1/2*(M1-3)/(M1-1)^2*s00*s01*s10;
X3 = -1/2*s00^2;
X4 = -1/((M1-1)^2*(M2-1))*s01*s12*s21*s10;
X5 = 1/(M1-1)^2*s01*s11*s10;
X6 = 1/((M1-1)*(M2-1))*s02*s21*s10;
X7 = 1/((M1-1)*(M2-1))*s01*s12*s20;
X8 = -1/(M2-1)*s02*s20;
X9 = 1/(1-M1)*t01*s10;
X10 = 1/(1-M1)*s01*t10;
X11 = q00;

X1c = 3/(M1-1)^4*s01^2*s10^2;
X2c = - 1/(M1-1)^3*s00*s01*s10;
X4c = - (3+M2*(M1-2)-2*M1+M1*(-2+M2+M1))/((M1-1)^4*(M2-1)^2)*s01*s12*s21*s10;
X5c = 2/(M1-1)^3*s01*s11*s10;
X6c = (-2+M1+M2)/((M1-1)^2*(M2-1)^2)*s02*s21*s10;
X7c = (-2+M1+M2)/((M1-1)^2*(M2-1)^2)*s01*s12*s20;
X8c = -1/(M2-1)^2*s02*s20;
X9c = -1/(1-M1)^2*t01*s10;
X10c = -1/(1-M1)^2*s01*t10;

R4l = 1/del*(24*X1+24*X2+6*X3+24*X4+12*X5+12*X6+12*X7+6*X8+4*X9+4*X10+X11);

if (getR4c)
    R4c = (24*X1c+24*X2c+24*X4c+12*X5c+12*X6c+12*X7c+6*X8c+4*X9c+4*X10c);
    moms = [R2l,R2c,R4l,R4c];
else
    moms = [R2l,R2c,R4l];
end
end