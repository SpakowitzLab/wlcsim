PROGRAM test_interp
IMPLICIT NONE
DOUBLE PRECISION RBIN(3)
Double PRECISION LBox
INTEGER IX(2), IY(2), IZ(2)
DOUBLE PRECISION WX(2), WY(2), WZ(2)
INTEGER NBINX
Double precision DEL
INTEGER confineType

RBIN(1)=1.5000-0.00001
RBIN(2)=2.0000
RBIN(3)=2.0000

LBox=4.0
NBINX=4
confineType=3
DEL=1.00000000

print*, "hello world"
call MC_interp(IX,IY,IZ,WX,WY,WZ,RBIN,LBOX,confineType,DEL)

print*, "WX",WX(1),WX(2)
print*, "WY",WY(1),WY(2)
print*, "WZ",WZ(1),WZ(2)

print*, "IX",IX(1),IX(2)
print*, "IY",IY(1),IY(2)
print*, "IZ",IZ(1),IZ(2)
end program
