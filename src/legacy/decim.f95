!---------------------------------------------------------------*

      subroutine decim(R,U,NT,N,NP,PARA,DT)

      PARAMETER (PI = 3.141593)

      real(dp) R(3,NT)  ! Bead positions
      real(dp) U(3,NT)  ! Tangent vectors
      real(dp) RT(NT,3)  ! Bead positions
      real(dp) UT(NT,3)  ! Tangent vectors
      real(dp) TTOT     ! Time of BD simulation
      integer N,NP,NT           ! Number of beads
      real(dp) PARA(10)
      real(dp) DEL
      real(dp) PVEC(60,8)
      integer inD,CRS
      real(dp) EB,EPAR,EPERP
      real(dp) GAM,ETA
      real(dp) XIR,XIU
      real(dp) LBOX     ! Box edge length
      real(dp) LHC      ! Length of HC int
      real(dp) VHC      ! HC strength
      real(dp) M
      real(dp) DT
      integer I,J,IP
      real(dp) L,LP

!     Reset the positions and orientations

      inD = 1
      do 10 IP = 1,NP
         I = 1
         J = 1
         do while (I <= N)
            RT(inD,1) = R(1,I + N*(IP-1))
            RT(inD,2) = R(2,I + N*(IP-1))
            RT(inD,3) = R(3,I + N*(IP-1))
            UT(inD,1) = U(1,I + N*(IP-1))
            UT(inD,2) = U(2,I + N*(IP-1))
            UT(inD,3) = U(3,I + N*(IP-1))
            I = I + 2
            J = J + 1
            inD = inD + 1
         enddo
 10   continue
      N = J-1

      do 20 I = 1,NP
         do 30 J = 1,N
            inD = J + N*(I-1)
            R(1,inD) = RT(inD,1)
            R(2,inD) = RT(inD,2)
            R(3,inD) = RT(inD,3)
            U(1,inD) = UT(inD,1)
            U(2,inD) = UT(inD,2)
            U(3,inD) = UT(inD,3)
 30      continue
 20   continue

!     Load in the parameters for the simulation

      open (unit = 5, file = 'input/input')
      read (unit = 5, fmt = '(4(/))')
      read (unit = 5, fmt = *) LP
      read (unit = 5, fmt = '(2(/))')
      read (unit = 5, fmt = *) L
      read (unit = 5, fmt = '(2(/))')
      read (unit = 5, fmt = *) LBOX
      read (unit = 5, fmt = '(2(/))')
      read (unit = 5, fmt = *) LHC
      read (unit = 5, fmt = '(2(/))')
      read (unit = 5, fmt = *) VHC
      close(5)
      L = L/LP
      DEL = L/(N-1.)

!     Load the tabulated parameters

      open (UNIT = 5,FILE = 'input/dssWLCparams',STATUS = 'OLD')
      do 40 I = 1,60
         READ(5,*) PVEC(I,1),PVEC(I,2),PVEC(I,3),PVEC(I,4),PVEC(I,5),PVEC(I,6),PVEC(I,7),PVEC(I,8)
 40   continue
      CLOSE(5)

      if (DEL < PVEC(1,1)) then
         DEL = PVEC(1,1)
      endif
      if (DEL > PVEC(60,1)) then
         DEL = PVEC(60,1)
      endif

      CRS = 0
      inD = 1
      do while (CRS == 0)
         if (DEL <= PVEC(inD,1)) then
            CRS = 1
         else
            inD = inD + 1
         endif
      enddo

      I = 2
      M = (PVEC(inD,I)-PVEC(inD-1,I))/(PVEC(inD,1)-PVEC(inD-1,1))
      EB = M*(DEL-PVEC(inD,1)) + PVEC(inD,I)

      I = 3
      M = (PVEC(inD,I)-PVEC(inD-1,I))/(PVEC(inD,1)-PVEC(inD-1,1))
      GAM = M*(DEL-PVEC(inD,1)) + PVEC(inD,I)

      I = 4
      M = (PVEC(inD,I)-PVEC(inD-1,I))/(PVEC(inD,1)-PVEC(inD-1,1))
      EPAR = M*(DEL-PVEC(inD,1)) + PVEC(inD,I)

      I = 5
      M = (PVEC(inD,I)-PVEC(inD-1,I))/(PVEC(inD,1)-PVEC(inD-1,1))
      EPERP = M*(DEL-PVEC(inD,1)) + PVEC(inD,I)

      I = 6
      M = (PVEC(inD,I)-PVEC(inD-1,I))/(PVEC(inD,1)-PVEC(inD-1,1))
      ETA = M*(DEL-PVEC(inD,1)) + PVEC(inD,I)

      I = 7
      M = (PVEC(inD,I)-PVEC(inD-1,I))/(PVEC(inD,1)-PVEC(inD-1,1))
      XIU = M*(DEL-PVEC(inD,1)) + PVEC(inD,I)

!      I = 8
!      M = (PVEC(inD,I)-PVEC(inD-1,I))/(PVEC(inD,1)-PVEC(inD-1,1))
!      DT = XIU*(M*(DEL-PVEC(inD,1)) + PVEC(inD,I))

      EB = EB/DEL
      EPAR = EPAR/DEL
      EPERP = EPERP/DEL
      GAM = DEL*GAM

      XIU = XIU*L/N
      XIR = L/N
      DT = 0.5*XIU/(EPERP*GAM**2.)

      PARA(1) = EB
      PARA(2) = EPAR
      PARA(3) = EPERP
      PARA(4) = GAM
      PARA(5) = ETA
      PARA(6) = XIR
      PARA(7) = XIU
      PARA(8) = LBOX
      PARA(9) = LHC
      PARA(10) = VHC

      RETURN
      END

!---------------------------------------------------------------*
