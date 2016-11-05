Subroutine MC_interp(IX,IY,IZ,WX,WY,WZ,RBIN,LBOX,confineType,DEL)
IMPLICIT NONE
Double PRECISION RBIN(3)
Double PRECISION LBox
INTEGER IX(2), IY(2), IZ(2)
DOUBLE PRECISION WX(2), WY(2), WZ(2)
INTEGER NBINX
Double precision DEL
INTEGER confineType

SELECT CASE (confineType)
CASE (0) ! Box from 0-LBOX, Bins split by boundaries
    ! Periodic BC
    RBIN(1)=RBIN(1)-nint(RBIN(1)/LBOX-0.5)*LBOX
    RBIN(2)=RBIN(2)-nint(RBIN(2)/LBOX-0.5)*LBOX
    RBIN(3)=RBIN(3)-nint(RBIN(3)/LBOX-0.5)*LBOX

    ! Binning  
    IX(1)=nint(RBIN(1)/DEL+0.5)
    IY(1)=nint(RBIN(2)/DEL+0.5)
    IZ(1)=nint(RBIN(3)/DEL+0.5)
    
    IX(2)=IX(1)-1
    IY(2)=IY(1)-1
    IZ(2)=IZ(1)-1
    
    ! Calculate the bin weighting
    WX(2)=(DEL*IX(1)-RBIN(1))/DEL   ! WX(2)=(RBIN(1)-IX(1)*DEL)/(IX(2)*DEL-IX(1)*DEL)
    WX(1)=1-WX(2)                   ! WX(1)=(IX(2)*DEL-RBIN(1))/(IX(2)*DEL-IX(1)*DEL)
    WY(1)=(DEL*IY(1)-RBIN(2))/DEL   ! WY(2)=(RBIN(2)-IY(1)*DEL)/(IY(2)*DEL-IY(1)*DEL)
    WY(2)=1-WY(1)                   ! WY(1)=(IY(2)*DEL-RBIN(2))/(IY(2)*DEL-IY(1)*DEL)
    WZ(2)=(DEL*IZ(1)-RBIN(3))/DEL   ! WZ(2)=(RBIN(3)-IZ(1)*DEL)/(IZ(2)*DEL-IZ(1)*DEL)
    WZ(1)=1-WZ(2)                   ! WZ(1)=(IZ(2)*DEL-RBIN(3))/(IZ(2)*DEL-IZ(1)*DEL)

    ! Periodic BC on Bins:
    IX(1)=IX(1)-floor(REAL((IX(1)-1))/REAL(NBINX)) * NBINX
    IX(2)=IX(2)-floor(REAL((IX(2)-1))/REAL(NBINX)) * NBINX
    IY(1)=IY(1)-floor(REAL((IY(1)-1))/REAL(NBINX)) * NBINX
    IY(2)=IY(2)-floor(REAL((IY(2)-1))/REAL(NBINX)) * NBINX
    IZ(1)=IZ(1)-floor(REAL((IZ(1)-1))/REAL(NBINX)) * NBINX
    IZ(2)=IZ(2)-floor(REAL((IZ(2)-1))/REAL(NBINX)) * NBINX
CASE (1)
    ! Periodic BC
    RBIN(1)=RBIN(1)-nint(RBIN(1)/LBOX-0.5)*LBOX
    RBIN(2)=RBIN(2)-nint(RBIN(2)/LBOX-0.5)*LBOX

    ! Binning  
    IX(1)=nint(RBIN(1)/DEL+0.5)
    IY(1)=nint(RBIN(2)/DEL+0.5)
    IZ(1)=nint(RBIN(3)/DEL+1.0) ! Note 1.0 so that box centers are on half intigers 
    
    IX(2)=IX(1)-1
    IY(2)=IY(1)-1
    IZ(2)=IZ(1)-1
    
    ! Calculate the bin weighting
    WX(2)=(DEL*IX(1)-RBIN(1))/DEL   ! WX(2)=(RBIN(1)-IX(1)*DEL)/(IX(2)*DEL-IX(1)*DEL)
    WX(1)=1-WX(2)                   ! WX(1)=(IX(2)*DEL-RBIN(1))/(IX(2)*DEL-IX(1)*DEL)
    WY(1)=(DEL*IY(1)-RBIN(2))/DEL   ! WY(2)=(RBIN(2)-IY(1)*DEL)/(IY(2)*DEL-IY(1)*DEL)
    WY(2)=1-WY(1)                   ! WY(1)=(IY(2)*DEL-RBIN(2))/(IY(2)*DEL-IY(1)*DEL)
    WZ(2)=(DEL*IZ(1)-0.5*DEL-RBIN(3))/DEL   ! WZ(2)=(RBIN(3)-IZ(1)*DEL)/(IZ(2)*DEL-IZ(1)*DEL)
    WZ(1)=1-WZ(2)                   ! WZ(1)=(IZ(2)*DEL-RBIN(3))/(IZ(2)*DEL-IZ(1)*DEL)

    if ((WZ(1).lt.0).OR.(WZ(2).lt.0)) then
        print*, "negitive W"
        stop 1
    endif

    ! Periodic BC on Bins:
    IX(1)=IX(1)-floor(REAL((IX(1)-1))/REAL(NBINX)) * NBINX
    IX(2)=IX(2)-floor(REAL((IX(2)-1))/REAL(NBINX)) * NBINX
    IY(1)=IY(1)-floor(REAL((IY(1)-1))/REAL(NBINX)) * NBINX
    IY(2)=IY(2)-floor(REAL((IY(2)-1))/REAL(NBINX)) * NBINX
CASE (2) ! Box confinement
    ! Binning  
    IX(1)=nint(RBIN(1)/DEL+1.0)
    IY(1)=nint(RBIN(2)/DEL+1.0)
    IZ(1)=nint(RBIN(3)/DEL+1.0) ! Note 1.0 so that box centers are on half intigers 
    
    IX(2)=IX(1)-1
    IY(2)=IY(1)-1
    IZ(2)=IZ(1)-1
     
    ! Calculate the bin weighting
    WX(2)=(DEL*IX(1)-0.5*DEL-RBIN(1))/DEL   ! WX(2)=(RBIN(1)-IX(1)*DEL)/(IX(2)*DEL-IX(1)*DEL)
    WX(1)=1-WX(2)                   ! WX(1)=(IX(2)*DEL-RBIN(1))/(IX(2)*DEL-IX(1)*DEL)
    WY(1)=(DEL*IY(1)-0.5*DEL-RBIN(2))/DEL   ! WY(2)=(RBIN(2)-IY(1)*DEL)/(IY(2)*DEL-IY(1)*DEL)
    WY(2)=1-WY(1)                   ! WY(1)=(IY(2)*DEL-RBIN(2))/(IY(2)*DEL-IY(1)*DEL)
    WZ(2)=(DEL*IZ(1)-0.5*DEL-RBIN(3))/DEL   ! WZ(2)=(RBIN(3)-IZ(1)*DEL)/(IZ(2)*DEL-IZ(1)*DEL)
    WZ(1)=1-WZ(2)                   ! WZ(1)=(IZ(2)*DEL-RBIN(3))/(IZ(2)*DEL-IZ(1)*DEL)
CASE (3)
    print*, "DEL", DEL," RBIN(1)",RBIN(1)
    ! Binning  
    IX(2)=nint(RBIN(1)/DEL)
    IY(2)=nint(RBIN(2)/DEL)
    IZ(2)=nint(RBIN(3)/DEL) ! Note 1.0 so that box centers are on half intigers 
    !write(*,"(A,I4,A,I4,A,I4)"),"IX(1):",IX(1),"  IY(1):",IY(1),"  IZ(1)",IZ(1)
        
    IX(1)=IX(2)+1
    IY(1)=IY(2)+1
    IZ(1)=IZ(2)+1
     
    ! Calculate the bin weighting
    WX(2)=(DEL*IX(1)-0.5*DEL-RBIN(1))/DEL   ! WX(2)=(RBIN(1)-IX(1)*DEL)/(IX(2)*DEL-IX(1)*DEL)
    WX(1)=1-WX(2)                   ! WX(1)=(IX(2)*DEL-RBIN(1))/(IX(2)*DEL-IX(1)*DEL)
    WY(1)=(DEL*IY(1)-0.5*DEL-RBIN(2))/DEL   ! WY(2)=(RBIN(2)-IY(1)*DEL)/(IY(2)*DEL-IY(1)*DEL)
    WY(2)=1-WY(1)                   ! WY(1)=(IY(2)*DEL-RBIN(2))/(IY(2)*DEL-IY(1)*DEL)
    WZ(2)=(DEL*IZ(1)-0.5*DEL-RBIN(3))/DEL   ! WZ(2)=(RBIN(3)-IZ(1)*DEL)/(IZ(2)*DEL-IZ(1)*DEL)
    WZ(1)=1-WZ(2)                   ! WZ(1)=(IZ(2)*DEL-RBIN(3))/(IZ(2)*DEL-IZ(1)*DEL) 
    !write(*,"(A,f4.3,A,f4.3,A,f4.3)"),"WX(1):",WX(1),"  WY(1):",WY(1),"  WZ(1)",WZ(1)
END SELECT

end subroutine
