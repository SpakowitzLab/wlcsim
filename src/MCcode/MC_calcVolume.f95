!-----------------------------------------------------------!
!
!     Calculate volume of each Bin within the confinement
!
!            Started by Quinn 3/17/15
!
!
!  
! confineType  |  Discription
! _____________|_________________________________
!    0         |  No confinement
!    1         |  Betwene two plates in Z direction at 0 and LBox
!    2         |  Cube of size LBox**3,  range: 0-LBox
!    3         |  Circle of radius LBox/2 inside box of size LBox

SUBROUTINE MC_caclVolume(confineType,NBINX,DEL, LBox, &
                         Vol,rand_stat)


!use mt19937, only : grnd, init_genrand, rnorm, mt, mti
use mersenne_twister
use setPrecision
IMPLICIT NONE

INTEGER confineType  ! Specifier for type of confinement
DOUBLE PRECISION LBox ! Side length of box
INTEGER I      ! for loops
INTEGER ix,iy,iz      ! location of conder
DOUBLE PRECISION x,y,z  
INTEGER :: STATUS = 0
INTEGER NBINX  ! length of side of box as an integer number of bins
DOUBLE PRECISION Rsqrd 
DOUBLE PRECISION Vol(NBINX**3)  ! output: volume of bins
DOUBLE PRECISION V      
DOUBLE PRECISION corner(8,3)
INTEGER nc
INTEGER, PARAMETER:: npts = 10000 
DOUBLE PRECISION DEL           ! side length of bins
DOUBLE PRECISION rsq, minr
type(random_stat) rand_stat !for random numer generator
real urand(3)

if (abs(DEL*NBINX-LBOX).gt.0.000001_dp) then 
    print*, "DEL=", DEL
    print*, "NBINX=",NBINX
    print*, "LBOX=",LBOX
    print*, "Error in MC_calcvolume, make box integer lenth*DEL"
    STOP 1
endif

if (confineType.EQ.0) then
    print*, "Don't call MC_calcVolume with this type of boundary"
    STATUS=1
    STOP 1
elseif(confineType.EQ.1) then
    print*, "Don't call MC_calcVolume with this type of boundary"
    STATUS=1
    STOP 1
elseif(confineType.EQ.2) then
    print*, "Don't call MC_calcVolume with this type of boundary"
    STATUS=1
    STOP 1
elseif(confineType.EQ.3) then
    Rsqrd=(LBox/2.0_dp)**2
    Do ix=1,NBINX
        Do iy=1,NBINX
            do iz=1,NBINX  
                x=DEL*ix
                y=DEL*iy
                z=DEL*iz
                corner(1,1)=x;    corner(1,2)=y;    corner(1,3)=z
                corner(2,1)=x;    corner(3,2)=y-DEL;corner(3,3)=z-DEL
                corner(3,1)=x;    corner(2,2)=y;    corner(2,3)=z-DEL
                corner(4,1)=x-DEL;corner(4,2)=y-DEL;corner(4,3)=z-DEL
                corner(5,1)=x-DEL;corner(5,2)=y;    corner(5,3)=z-DEL
                corner(6,1)=x;    corner(6,2)=y-DEL;corner(6,3)=z
                corner(7,1)=x-DEL;corner(7,2)=y;    corner(7,3)=z
                corner(8,1)=x-DEL;corner(8,2)=y-DEL;corner(8,3)=z
                
                nc=0
                minr=Rsqrd+2*DEL ! just big
                do I=1,8
                    rsq=((corner(I,1)-LBox/2.0_dp)**2+ & 
                        (corner(I,2)-LBox/2.0_dp)**2+ &
                        (corner(I,3)-LBox/2.0_dp)**2)
                    minr=min(minr,rsq)
                    if (rsq.gt.Rsqrd) then
                        nc=nc+1
                    endif   
                enddo
                if (nc.eq.0) then
                    ! inside
                    V=DEL**3
                elseif ((nc.eq.8).and.(minr.gt.Rsqrd+DEL*0.5_dp)) then
                    V=0.0_dp
                else
                    V=0.0_dp
                    do I=1,npts
                        call random_number(urand,rand_stat)
                        rsq=( (corner(1,1)-urand(1)*DEL-LBox/2.0_dp)**2+ &
                              (corner(1,2)-urand(2)*DEL-LBox/2.0_dp)**2+ &
                              (corner(1,3)-urand(3)*DEL-LBox/2.0_dp)**2)
                        if (rsq.lt. Rsqrd)  then
                            V=V+1.0_dp
                        endif
                    enddo
                    V=(Del**3)*V/dble(npts)
                endif
                Vol(ix+(iy-1)*NBINX+(iz-1)*NBINX**2)=V
            enddo
        enddo              
    enddo 
else 
   print*, "Undefined confine Type"
   stop 1
endif




END
