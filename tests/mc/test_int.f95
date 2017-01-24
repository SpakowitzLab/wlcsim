!---------------------------------------------------------------*
      
PROGRAM scratch
      
!     

      use mt19937, only : grnd, init_genrand, rnorm, mt, mti

  IMPLICIT NONE
!   Random number generator initiation
    integer IDUM
    character*8 datedum
    character*10 timedum
    character*5 zonedum
    integer seedvalues(8)


      
      INTEGER, PARAMETER :: NT=1
      INTEGER, PARAMETER :: NBIN=27
      INTEGER, PARAMETER :: NBINX=3

      DOUBLE PRECISION R(NT,3)  ! Conformation of polymer chains
      DOUBLE PRECISION RP(NT,3) ! Conformation of polymer chains
      INTEGER AB(NT)            ! Chemical identity of beads    
      INTEGER ABP(NT)           ! Proposed chemical identities
      
!     Simulation input variables
      
      INTEGER confineType ! Boundary conditions
      DOUBLE PRECISION V        ! Monomer volume
      DOUBLE PRECISION HP1_Bind ! HP1_Binding energy
      DOUBLE PRECISION CHI      ! Solvent polymer chi parameter value
      DOUBLE PRECISION KAP      ! Compressibility value
      DOUBLE PRECISION LBOX     ! Simulation box size (approximate)
      DOUBLE PRECISION DEL      ! Discretization size (approximate)
      INTEGER POLY              ! Polydisperse (step-growth stat)
      
!     Variables for density calculation
      
      DOUBLE PRECISION PHIA(NBIN) ! Volume fraction of A
      DOUBLE PRECISION PHIB(NBIN) ! Volume fraction of B
      DOUBLE PRECISION DPHIA(NBIN) ! Delta Volume fraction of A
      DOUBLE PRECISION DPHIB(NBIN) ! Delta Volume fraction of B
      INTEGER INDPHI(NBIN)      ! Indices of the phi
      INTEGER NPHI              ! Number of phi values that change 
      DOUBLE PRECISION DEINT    ! Change in Self-interaction energy
      INTEGER I1                ! Test bead position 1
      INTEGER I2                ! Test bead position 2
      DOUBLE PRECISION Vol(NBIN)  ! Volume of bins 

      INTEGER ii
      INTEGER IX,IY,IZ
      I1=1
      I2=NT
      Do ii=1,NBIN
          Vol(ii)=1.0
          PHIA(ii)=0.0
          PHIB(ii)=0.0
      enddo    
      DEL=1.0
      LBOX=3.0
      
      R(1,1)= 0.0001; R(1,2)= 0.0001; R(1,3)= 0.0001
      RP(1,1)= 2.999; RP(1,2)= 2.999; RP(1,3)= 2.999
      AB(1)=1; ABP=1;
      V=1.0
      CHI=0.0
      KAP=0.0
      HP1_Bind=1
      confineType=3
      
      
      
    !   Seed the random number generator off the computer clock

    call date_and_time(datedum,timedum,zonedum,seedvalues)
    ! concatenate filename, time within mins, secs, millisecs to seed random number generator	
    IDUM=-seedvalues(5)*1E7-seedvalues(6)*1E5-seedvalues(7)*1E3-seedvalues(8)
    call init_genrand(IDUM)

    
!    print*, "hello world"
!    print*, grnd()

!    wdble=0.0
!    wInt=0.0
!    wInt2=100
!    print*, floor(REAL(wInt)/REAL(wInt2))

      
      call MC_int(DEINT,R,AB,NT,NBIN, &
                        V,CHI,KAP,LBOX,DEL,PHIA,PHIB,DPHIA,DPHIB, &
                        INDPHI,NPHI,RP,I1,I2,HP1_Bind,ABP,confineType, &
                        Vol)

!      do ii=1,9
!          print*, DPHIA(1+ii), DPHIA(2+ii), DPHIA(3+ii)
!      enddo
      print*, " "
      do IZ=1,3
         print*, "z=",IZ
         do IY=1,3
             do IX=1,3
                 write(*,"(1f10.4,$)"),  DPHIA( IX+(IY-1)*NBINX+(IZ-1)*NBINX**2 )
             enddo
             print*, " "
         enddo
      enddo 
end program
