!---------------------------------------------------------------*
      
PROGRAM test_calcVolume
      
!     

      use mt19937, only : grnd, init_genrand, rnorm, mt, mti

  IMPLICIT NONE
!   Random number generator initiation
    integer IDUM
    character*8 datedum
    character*10 timedum
    character*5 zonedum
    integer seedvalues(8)

    integer, parameter :: NBINX=30
    integer confineType
    double precision DEL, lbox
    double precision Vol(NBINX*NBINX*NBINX)

    integer I  !for loops
    double precision Voltotal
    confineType=3 
    DEL=1.0
    lbox=30.0

    !   Seed the random number generator off the computer clock

    call date_and_time(datedum,timedum,zonedum,seedvalues)
    ! concatenate filename, time within mins, secs, millisecs to seed random number generator	
    IDUM=-seedvalues(5)*1E7-seedvalues(6)*1E5-seedvalues(7)*1E3-seedvalues(8)
    call init_genrand(IDUM)

        
    call MC_caclVolume(confineType,NBINX,DEL, LBox, &
                               Vol)
    Voltotal=0.0
    Do I=1,NBINX**3
        VolTotal=volTotal+Vol(I)
    enddo
    print*, "Vol total: ", VolTotal




      

end program
