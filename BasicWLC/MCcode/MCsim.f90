!---------------------------------------------------------------*

!     This subroutine performs a Monte Carlo simulation on the
!     polymer chain.

      SUBROUTINE MCsim(R,U,NT,N,NP,NSTEP,BROWN, &
           INTON,IDUM,PARA,MCAMP,SUCCESS,MOVEON,WINDOW,SIMTYPE)

      use mt19937, only : grnd

      PARAMETER (PI=3.141592654) ! Value of pi

      DOUBLE PRECISION R(NT,3)  ! Bead positions
      DOUBLE PRECISION U(NT,3)  ! Tangent vectors
      DOUBLE PRECISION RP(NT,3)  ! Bead positions
      DOUBLE PRECISION UP(NT,3)  ! Tangent vectors
      INTEGER N,NP,NT           ! Number of beads
      INTEGER NSTEP             ! Number of MC steps
      INTEGER BROWN            ! Turn on fluctuations
      INTEGER INTON             ! Include polymer interactions

!     Variables for the simulation

      INTEGER ISTEP             ! Current MC step index
      DOUBLE PRECISION PROB     ! Calculated test prob
      DOUBLE PRECISION TEST     ! Random test variable
      INTEGER IB                ! Test bead
      INTEGER IP                ! Test polymer
      INTEGER IB1               ! Test bead position 1
      INTEGER IT1               ! Index of test bead 1
      INTEGER IB2               ! Test bead position 2
      INTEGER IT2               ! Index of test bead 2

      INTEGER TEMP
      REAL ran1                 ! Random number generator
      INTEGER IDUM              ! Seed for the generator
      INTEGER NOW(3)            ! Time now (hr,min,sec)
      INTEGER I
      DOUBLE PRECISION R0(3)

!     Energy variables

      DOUBLE PRECISION DEELAS   ! Change in bending energy
      DOUBLE PRECISION DESELF   ! Change in self energy
      DOUBLE PRECISION DEEX     ! Change in external energy
      DOUBLE PRECISION ENERGY
      DOUBLE PRECISION DECOM    ! Change in the compression energy

!     MC adaptation variables

      DOUBLE PRECISION MCAMP(6) ! Amplitude of random change
      INTEGER MCTYPE            ! Type of MC move
      INTEGER NADAPT(6)         ! Num steps btwn adapt
      DOUBLE PRECISION PHIT     ! % hits per total steps
      DOUBLE PRECISION PDESIRE(6) ! Desired hit rate
      INTEGER SUCCESS(6)        ! Number of successes
      DOUBLE PRECISION MINAMP(6)! Minimum amp to stop
      DOUBLE PRECISION MAXAMP(6)! Minimum amp to stop
      INTEGER MOVEON(6)        ! Is the move active
      INTEGER WINDOW(6)            ! Size of window for bead selection
      INTEGER SIMTYPE           ! Simulation method (WLC=1,SSWLC=2,GC=3)

!     Variables in the simulation

      DOUBLE PRECISION EB,EPAR,EPERP
      DOUBLE PRECISION GAM,ETA
      DOUBLE PRECISION XIR,XIU
      DOUBLE PRECISION LBOX     ! Box edge length
      DOUBLE PRECISION LHC      ! Length of HC int
      DOUBLE PRECISION VHC      ! HC strength
      DOUBLE PRECISION PARA(10)
      DOUBLE PRECISION FCOM

!     Load the input parameters

      EB=PARA(1)
      EPAR=PARA(2)
      EPERP=PARA(3)
      GAM=PARA(4)
      ETA=PARA(5)
      XIR=PARA(6)
      XIU=PARA(7)
      LBOX=PARA(8)
      FCOM=PARA(9)
      VHC=PARA(10)

      MINAMP(1)=0.0*PI
      MINAMP(2)=0.0
      MINAMP(3)=0.0*PI
      MINAMP(4)=0.0*PI
      MINAMP(5)=0.1*PI
      MINAMP(6)=0.01

      MAXAMP(1)=2.*PI
      MAXAMP(2)=LBOX
      MAXAMP(3)=2.*PI
      MAXAMP(4)=2.*PI
      MAXAMP(5)=2.*PI
      MAXAMP(6)=LBOX

      NADAPT(1)=1000
      NADAPT(2)=1000
      NADAPT(3)=1000
      NADAPT(4)=1000
      NADAPT(5)=1000
      NADAPT(6)=1000
      if (NSTEP.LE.NADAPT(1)) then
         NADAPT(1)=NSTEP
      endif
      if (NSTEP.LE.NADAPT(2)) then
         NADAPT(2)=NSTEP
      endif
      if (NSTEP.LE.NADAPT(3)) then
         NADAPT(3)=NSTEP
      endif
      if (NSTEP.LE.NADAPT(4)) then
         NADAPT(4)=NSTEP
      endif
      if (NSTEP.LE.NADAPT(5)) then
         NADAPT(5)=NSTEP
      endif

      PDESIRE(1)=0.5
      PDESIRE(2)=0.5
      PDESIRE(3)=0.5
      PDESIRE(4)=0.5
      PDESIRE(5)=0.5
      PDESIRE(6)=0.5

      SUCCESS(1)=0
      SUCCESS(2)=0
      SUCCESS(3)=0
      SUCCESS(4)=0
      SUCCESS(5)=0
      SUCCESS(6)=0

      DEELAS=0.
      DESELF=0.
      DEEX=0.


!     Begin Monte Carlo simulation

      ISTEP=1

      DO WHILE (ISTEP.LE.NSTEP)

         DO 10 MCTYPE=1,6

            if (MOVEON(MCTYPE).EQ.0) then
               goto 60
            endif

            call MC_move(R,U,RP,UP,NT,N,NP,IP,IB1,IB2,IT1,IT2,IDUM,MCTYPE,MCAMP)

            if (SIMTYPE.EQ.1.AND.abs(IB2-IB1).LE.1) then
               goto 60
            endif

!     Calculate the change in polymer elastic energy using
!     SIMTYPE indicates simulation method (WLC=1,SSWLC=2,GC=3)

            call MC_eelas(DEELAS,R,U,RP,UP,NT,N,NP,IP, &
                 IB1,IB2,IT1,IT2,EB,EPAR,EPERP,GAM,ETA,SIMTYPE)

!     Calculate the change in the self-interaction energy

            if (INTON.EQ.1) then
               call MC_self(DESELF,R,U,RP,UP,NT,N,NP,IP,IB1,IB2,IT1,IT2,LHC,VHC,LBOX,GAM)
            endif

!     Calculate the change in the external force energy

            DECOM=0.
            if (IB1.EQ.1.AND.IB2.EQ.N) then
               DECOM=FCOM*(sqrt((RP(IT2,1)-RP(IT1,1))**2.+(RP(IT2,2)-RP(IT1,2))**2.+(RP(IT2,3)-RP(IT1,3))**2.) &
                    -sqrt((R(IT2,1)-R(IT1,1))**2.+(R(IT2,2)-R(IT1,2))**2.+(R(IT2,3)-R(IT1,3))**2.))
            elseif (IB1.EQ.1.AND.IB2.NE.N) then
               DECOM=FCOM*(sqrt((R(N*IP,1)-RP(IT1,1))**2.+(R(N*IP,2)-RP(IT1,2))**2.+(R(N*IP,3)-RP(IT1,3))**2.) &
                    -sqrt((R(N*IP,1)-R(IT1,1))**2.+(R(N*IP,2)-R(IT1,2))**2.+(R(N*IP,3)-R(IT1,3))**2.))
            elseif (IB1.NE.1.AND.IB2.EQ.N) then
               DECOM=FCOM*(sqrt((RP(IT2,1)-R(1+N*(IP-1),1))**2.+(RP(IT2,2)-R(1+N*(IP-1),2))**2.+(RP(IT2,3)-R(1+N*(IP-1),3))**2.) &
                    -sqrt((R(IT2,1)-R(1+N*(IP-1),1))**2.+(R(IT2,2)-R(1+N*(IP-1),2))**2.+(R(IT2,3)-R(1+N*(IP-1),3))**2.))
            endif

!     Change the position if appropriate

            ENERGY=DEELAS+DESELF+DECOM

            PROB=exp(-ENERGY)
            if (BROWN.EQ.1) then
               TEST=grnd()
            else
               TEST=1.
            endif
            if (TEST.LE.PROB) then
               DO 20 I=IT1,IT2
                  R(I,1)=RP(I,1)
                  R(I,2)=RP(I,2)
                  R(I,3)=RP(I,3)
                  U(I,1)=UP(I,1)
                  U(I,2)=UP(I,2)
                  U(I,3)=UP(I,3)
 20            CONTINUE
               SUCCESS(MCTYPE)=SUCCESS(MCTYPE)+1
            endif

!     Adapt the amplitude of step every NADAPT steps

            if (mod(ISTEP,NADAPT(MCTYPE)).EQ.0) then
               PHIT=real(SUCCESS(MCTYPE))/real(NADAPT(MCTYPE))
               if (PHIT.GT.PDESIRE(MCTYPE)) then
                  MCAMP(MCTYPE)=MCAMP(MCTYPE)*1.05
               else
                  MCAMP(MCTYPE)=MCAMP(MCTYPE)*0.95
               endif
               if (MCAMP(MCTYPE).GT.MAXAMP(MCTYPE)) then
                  MCAMP(MCTYPE)=MAXAMP(MCTYPE)
               endif
               if (MCAMP(MCTYPE).LT.MINAMP(MCTYPE)) then
                  MCAMP(MCTYPE)=MINAMP(MCTYPE)
               endif

               SUCCESS(MCTYPE)=0

               IB=1
               DO 40 I=1,NP
                  R0(1)=nint(R(IB,1)/LBOX-0.5)*LBOX
                  R0(2)=nint(R(IB,2)/LBOX-0.5)*LBOX
                  R0(3)=nint(R(IB,3)/LBOX-0.5)*LBOX
                  DO 50 J=1,N
                     R(IB,1)=R(IB,1)-R0(1)
                     R(IB,2)=R(IB,2)-R0(2)
                     R(IB,3)=R(IB,3)-R0(3)
                     IB=IB+1
 50              CONTINUE
 40           CONTINUE

            endif

 60         CONTINUE


 10      CONTINUE

         ISTEP=ISTEP+1

      ENDDO

      RETURN
      END

!---------------------------------------------------------------*
