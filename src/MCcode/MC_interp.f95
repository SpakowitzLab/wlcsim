!-----------------------------------------------------------------
!
!     This program linearly interpolates a bead at RBIN into
!     8 bins indexed by IX, IY, IZ with weights WX, WY, WZ
!
!        Addapted from MC_int.f95
!        By Quinn MacPherson in summer 2016
!
!---------------------------------------------------------------
subroutine interp(confineType,RBIN,LBOX,NBINX,DEL,IX,IY,IZ,WX,WY,WZ)
use setPrecision
IMPLICIT NONE
integer, intent (in) :: confineType
DOUBLE PRECISION, intent(inout) :: RBIN(3) ! position or posiion within bin 
DOUBLE PRECISION, intent(in) :: LBOX(3) ! Side length of box
INTEGER, intent(in) :: NBINX(3)      ! number of discritations in each direction
DOUBLE PRECISION, intent(in) :: DEL  ! size of discritation
INTEGER, intent(out) :: IX(2)  ! Output
INTEGER, intent(out) :: IY(2) ! Output
INTEGER, intent(out) :: IZ(2)  ! Output
DOUBLE PRECISION, intent(out) :: WX(2) ! Output
DOUBLE PRECISION, intent(out) :: WY(2) ! Output
DOUBLE PRECISION, intent(out) :: WZ(2) ! Output
SELECT CASE (confineType)
CASE (0) ! Box from 0-LBOX, Bins split by boundaries
    ! Periodic BC
    RBIN(1)=RBIN(1)-floor(RBIN(1)/LBOX(1))*LBOX(1)
    RBIN(2)=RBIN(2)-floor(RBIN(2)/LBOX(2))*LBOX(2)
    RBIN(3)=RBIN(3)-floor(RBIN(3)/LBOX(3))*LBOX(3)

    ! Binning  
    IX(1)=ceiling(RBIN(1)/DEL)
    IY(1)=ceiling(RBIN(2)/DEL)
    IZ(1)=ceiling(RBIN(3)/DEL)
    
    IX(2)=IX(1)-1
    IY(2)=IY(1)-1
    IZ(2)=IZ(1)-1
    
    ! Calculate the bin weighting
    WX(2)=(DEL*IX(1)-RBIN(1))/DEL   
    WX(1)=1.0_dp-WX(2)              
    WY(2)=(DEL*IY(1)-RBIN(2))/DEL   
    WY(1)=1.0_dp-WY(2)              
    WZ(2)=(DEL*IZ(1)-RBIN(3))/DEL   
    WZ(1)=1.0_dp-WZ(2)              

    ! Periodic BC on Bins:
    IX(1)=IX(1)-floor(REAL((IX(1)-1))/REAL(NBINX(1))) * NBINX(1)
    IX(2)=IX(2)-floor(REAL((IX(2)-1))/REAL(NBINX(1))) * NBINX(1)
    IY(1)=IY(1)-floor(REAL((IY(1)-1))/REAL(NBINX(2))) * NBINX(2)
    IY(2)=IY(2)-floor(REAL((IY(2)-1))/REAL(NBINX(2))) * NBINX(2)
    IZ(1)=IZ(1)-floor(REAL((IZ(1)-1))/REAL(NBINX(3))) * NBINX(3)
    IZ(2)=IZ(2)-floor(REAL((IZ(2)-1))/REAL(NBINX(3))) * NBINX(3)
CASE (1)
    ! Periodic BC
    RBIN(1)=RBIN(1)-floor(RBIN(1)/LBOX(2))*LBOX(1)
    RBIN(2)=RBIN(2)-floor(RBIN(2)/LBOX(1))*LBOX(2)

    ! Binning  
    IX(1)=ceiling(RBIN(1)/DEL)
    IY(1)=ceiling(RBIN(2)/DEL)
    IZ(1)=nint(RBIN(3)/DEL)+1 ! Note 1.0 so that box centers are on half intigers 
    
    IX(2)=IX(1)-1
    IY(2)=IY(1)-1
    IZ(2)=IZ(1)-1
    
    ! Calculate the bin weighting
    WX(2)=(DEL*IX(1)-RBIN(1))/DEL   ! WX(2)=(RBIN(1)-IX(1)*DEL)/(IX(2)*DEL-IX(1)*DEL)
    WX(1)=1.0_dp-WX(2)              ! WX(1)=(IX(2)*DEL-RBIN(1))/(IX(2)*DEL-IX(1)*DEL)
    WY(2)=(DEL*IY(1)-RBIN(2))/DEL   ! WY(2)=(RBIN(2)-IY(1)*DEL)/(IY(2)*DEL-IY(1)*DEL)
    WY(1)=1.0_dp-WY(2)              ! WY(1)=(IY(2)*DEL-RBIN(2))/(IY(2)*DEL-IY(1)*DEL)
    WZ(2)=(DEL*IZ(1)-0.5_dp*DEL-RBIN(3))/DEL   ! WZ(2)=(RBIN(3)-IZ(1)*DEL)/(IZ(2)*DEL-IZ(1)*DEL)
    WZ(1)=1.0_dp-WZ(2)                   ! WZ(1)=(IZ(2)*DEL-RBIN(3))/(IZ(2)*DEL-IZ(1)*DEL)

    if ((WZ(1).lt.0).OR.(WZ(2).lt.0)) then
        print*, "negitive W"
        stop 1
    endif

    ! Periodic BC on Bins:
    IX(1)=IX(1)-floor(REAL((IX(1)-1))/REAL(NBINX(1))) * NBINX(1)
    IX(2)=IX(2)-floor(REAL((IX(2)-1))/REAL(NBINX(1))) * NBINX(1)
    IY(1)=IY(1)-floor(REAL((IY(1)-1))/REAL(NBINX(2))) * NBINX(2)
    IY(2)=IY(2)-floor(REAL((IY(2)-1))/REAL(NBINX(2))) * NBINX(2)
CASE (2) ! Box confinement
    ! Binning  
    IX(1)=nint(RBIN(1)/DEL)+1 
    IY(1)=nint(RBIN(2)/DEL)+1
    IZ(1)=nint(RBIN(3)/DEL)+1 ! Note +1 because fortran starts a 1 
    
    IX(2)=IX(1)-1
    IY(2)=IY(1)-1
    IZ(2)=IZ(1)-1
     
    ! Calculate the bin weighting
    WX(2)=(DEL*IX(1)-0.5_dp*DEL-RBIN(1))/DEL
    WX(1)=1.0_dp-WX(2)             
    WY(2)=(DEL*IY(1)-0.5_dp*DEL-RBIN(2))/DEL 
    WY(1)=1.0_dp-WY(2)                           
    WZ(2)=(DEL*IZ(1)-0.5_dp*DEL-RBIN(3))/DEL    
    WZ(1)=1.0_dp-WZ(2)                  
CASE (3) 
    ! Binning  
    IX(1)=nint(RBIN(1)/DEL)+1 
    IY(1)=nint(RBIN(2)/DEL)+1
    IZ(1)=nint(RBIN(3)/DEL)+1 ! Note +1 because fortran starts a 1 
        
    IX(2)=IX(1)-1
    IY(2)=IY(1)-1
    IZ(2)=IZ(1)-1
     
    ! Calculate the bin weighting
    WX(2)=(DEL*IX(1)-0.5_dp*DEL-RBIN(1))/DEL  
    WX(1)=1.0_dp-WX(2)                          
    WY(2)=(DEL*IY(1)-0.5_dp*DEL-RBIN(2))/DEL  
    WY(1)=1.0_dp-WY(2)                          
    WZ(2)=(DEL*IZ(1)-0.5_dp*DEL-RBIN(3))/DEL  
    WZ(1)=1.0_dp-WZ(2)                           
CASE (4) ! Box from 0-LBOX, Bins split by boundaries
    ! Periodic BC
    RBIN(1)=RBIN(1)-floor(RBIN(1)/LBOX(1))*LBOX(1)
    RBIN(2)=RBIN(2)-floor(RBIN(2)/LBOX(2))*LBOX(2)
    RBIN(3)=RBIN(3)-floor(RBIN(3)/LBOX(3))*LBOX(3)

    ! Binning  
    IX(1)=ceiling(RBIN(1)/DEL)
    IY(1)=ceiling(RBIN(2)/DEL)
    IZ(1)=ceiling(RBIN(3)/DEL)
    
    IX(2)=IX(1)-1
    IY(2)=IY(1)-1
    IZ(2)=IZ(1)-1
    
    ! Calculate the bin weighting
    WX(2)=(DEL*IX(1)-RBIN(1))/DEL   
    WX(1)=1.0_dp-WX(2)              
    WY(2)=(DEL*IY(1)-RBIN(2))/DEL   
    WY(1)=1.0_dp-WY(2)              
    WZ(2)=(DEL*IZ(1)-RBIN(3))/DEL   
    WZ(1)=1.0_dp-WZ(2)              

    ! Periodic BC on Bins:
    IX(1)=IX(1)-floor(REAL((IX(1)-1))/REAL(NBINX(1))) * NBINX(1)
    IX(2)=IX(2)-floor(REAL((IX(2)-1))/REAL(NBINX(1))) * NBINX(1)
    IY(1)=IY(1)-floor(REAL((IY(1)-1))/REAL(NBINX(2))) * NBINX(2)
    IY(2)=IY(2)-floor(REAL((IY(2)-1))/REAL(NBINX(2))) * NBINX(2)
    IZ(1)=IZ(1)-floor(REAL((IZ(1)-1))/REAL(NBINX(3))) * NBINX(3)
    IZ(2)=IZ(2)-floor(REAL((IZ(2)-1))/REAL(NBINX(3))) * NBINX(3)
END SELECT
return
end subroutine
