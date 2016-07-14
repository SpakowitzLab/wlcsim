function Ycouple = getYcouple(LMAX,savefile)
% --------------------
%% tabulate coupling coefficients for the spherical harmonics
% LMAX is the maximal L index
% optional: if savefile supplied, will save results in a .mat file for
% future use

% specifically, Ycouple contains the coefficients for decomposing a 
% product of spherical harmonics in terms of single spherical harmonics
% Y_l1^m1 Y_l2^m2 = sum(Ycouple{ind1}(ind2,ind3) *Y_l3^(j1+j2))
% where ind1 = l1^2+j1+l1+1
% WARNING: this is not at all efficient, since it doesn't use recursion
% formulas for the wigner  3j coefficients
% WARNING: this makes use of the vanilla factorial function in matlab
% and can thus lead to numerical errors for indices > 21


tic
% maximal indices
Imax = ljind(LMAX,LMAX);
Imax3 = ljind(2*LMAX,2*LMAX);
% set up an array of sparse matrices
W3Jsparse = cell(Imax,1);
for ind1 = 1:Imax
    W3Jsparse{ind1} = sparse(zeros(Imax,Imax3));
    Ycouple{ind1} = sparse(zeros(Imax,Imax3));
end

%% pretabulate factorials
factsave = zeros(4*LMAX+1,1);
for c = 1:4*LMAX+1
    factsave(c) = factorial(c);
end
%%
% get the Wigner3j coefficients
for l1 = 0:LMAX
    l1
    for j1 = -l1:l1
        ind1 = ljind(l1,j1);
        [l1,j1]
        
        for l2 = 0:LMAX            
            for j2 = -l2:l2
                ind2 = ljind(l2,j2);
                
                for l3 = abs(l1-l2):l1+l2
                    j3 = j1+j2;
                    if abs(j3)>l3
                        continue
                    end
                    ind3 = ljind(l3,j3);
                    if (mod(l1+l2+l3,2)>0); continue; end
                    
                    W3Jsparse{ind1}(ind2,ind3) = Wigner3j(l1,l2,l3,j1,j2,-j3,factsave);                     
                end
            end
        end
    end
end
toc

% -----------------------
%% tabulate coupling coefficients
tic
Imax = ljind(LMAX,LMAX);
    
for l1 = 0:LMAX
    l1
    for j1 = -l1:l1
        ind1 = ljind(l1,j1);       
        for l2 = 0:LMAX
            for j2 = -l2:l2
                ind2 = ljind(l2,j2);
                
                for l3 = abs(l1-l2):min(LMAX,l1+l2)
                    j3 = j1+j2;
                    if abs(j3)>l3
                        continue
                    end
                    ind3 = ljind(l3,j3);
                    if (mod(l1+l2+l3,2)>0); continue; end
                    
                    C1 = W3Jsparse{ljind(l1,0)}(ljind(l2,0),ljind(l3,0));                    
                    C2 = W3Jsparse{ind1}(ind2,ind3);
                  
                    Ycouple{ind1}(ind2,ind3) = sqrt((2*l1+1)*(2*l2+1)*(2*l3+1)/(4*pi))...
                        *(-1)^(j1+j2)*C1*C2;                                                         
                end
            end
        end
    end
end
toc

if (nargin>1)
    save(savefile,'LMAX','Imax','Imax3','Ycouple','W3Jsparse')
end
end