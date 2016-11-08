! *---------------------------------------------------------------*
!
!     The parameters provided to the simulation describe the "real" polymer
!     chain that we want to simulate and how to discretize it. This procedure
!     calculates the parameters of the corresponding discrete chain (i.e.
!     renormalized parameters) that we will actually be simulating.
!
!     The parameters are chosen to reproduce the
!     desired end-to-end distribution of that length of chain with the given
!     number of beads, and the drag coefficient is chosen such that, independent
!     of time step, the center of mass of each persistence length-long segment
!     of the chain diffuses one persistence length's distance per unit time.
!
!     For the Gaussian chain regime, this is derived in the documentation. For
!     the stretchable-shearable wormlike chain parameters are chosen according
!     to Elena's paper (Koslover and Spakowitz, "Discretizing elastic chains for
!     coarse-grained polymer models," Soft Matter, 2013).
!

      subroutine getpara(para,DT,SIMtype)

      real(dp) para(10)
      real(dp) del
      real(dp) pvec(679,8)
      integer ind,crs
      integer ring
      real(dp) eb,epar,eperp
      real(dp) gam,eta
      real(dp) xir,xiu
      real(dp) lbox     ! Box edge length
      real(dp) lhc      ! Length of HC int
      real(dp) rend     ! Fixed end-to-end distance (dimensionless)
      real(dp) vhc      ! HC strength
      real(dp) M
      real(dp) DT
      integer I,N
      real(dp) L,LP,LT
      integer Lk
      integer SIMtype           ! Simulation method (WLC=1,SSWLC=2,GC=3)

!     Load in the parameters for the simulation

      open (unit=5, file='input/input')
      read (unit=5, fmt='(4(/))')
      read (unit=5, fmt=*) LP
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) LT
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) L
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) LK
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) lbox
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) lhc
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) rend
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) vhc
      read (unit=5, fmt='(2(/))')
      read (unit=5, fmt=*) N
      close(5)
      if (N.EQ.1.0d0) then
          print*, 'Non-dimensionalization used requires at least two beads, 1 requested.'
          stop 1
      endif
      rend=rend*L/LP

      if (ring.EQ.0) then
         del=L/LP/(N-1.0_DP)
      else
         del=L/LP/(N-1.0_DP)
      endif

!     Load the tabulated parameters

      open (unit=5,file='input/dssWLCparams',status='OLD')
      do 10 I=1,679
         read(5,*) pvec(I,1),pvec(I,2),pvec(I,3),pvec(I,4),pvec(I,5),pvec(I,6),pvec(I,7),pvec(I,8)
 10   continue
      close(5)


!     Setup the parameters for WLC simulation

      if (del.LT.pvec(1,1)) then
         print*, 'It has never been known if the WLC code actually works.'
         print*, 'An entire summer student (Luis Nieves) was thrown at this'
         print*, 'problem and it is still not solved.'
         stop 1
         eb=LP/del
         gam=del
         xir=L/LP/N
         SIMtype=1
      endif

!    Setup the parameters for GC simulation

      if (del.GT.pvec(679,1)) then
         epar=1.5/del
         gam=0.0_dp
         SIMtype=3
         xir=L/N/LP
      endif

!    Setup the parameters for ssWLC simulation

      if (del.GE.pvec(1,1).AND.del.LE.pvec(679,1)) then
         SIMtype=2

        crs=0
        ind=1
        do while (crs.EQ.0)
            if (del.LE.pvec(ind,1)) then
                crs=1
            else
                ind=ind+1
            endif
        enddo

        !     Perform linear interpolations

        I=2
        M=(pvec(ind,I)-pvec(ind-1,I))/(pvec(ind,1)-pvec(ind-1,1))
        eb=M*(del-pvec(ind,1))+pvec(ind,I)

        I=3
        M=(pvec(ind,I)-pvec(ind-1,I))/(pvec(ind,1)-pvec(ind-1,1))
        gam=M*(del-pvec(ind,1))+pvec(ind,I)

        I=4
        M=(pvec(ind,I)-pvec(ind-1,I))/(pvec(ind,1)-pvec(ind-1,1))
        epar=M*(del-pvec(ind,1))+pvec(ind,I)

        I=5
        M=(pvec(ind,I)-pvec(ind-1,I))/(pvec(ind,1)-pvec(ind-1,1))
        eperp=M*(del-pvec(ind,1))+pvec(ind,I)

        I=6
        M=(pvec(ind,I)-pvec(ind-1,I))/(pvec(ind,1)-pvec(ind-1,1))
        eta=M*(del-pvec(ind,1))+pvec(ind,I)

        I=7
        M=(pvec(ind,I)-pvec(ind-1,I))/(pvec(ind,1)-pvec(ind-1,1))
        xiu=M*(del-pvec(ind,1))+pvec(ind,I)

        eb=LP*eb/(del*LP)
        epar=epar/(del*LP*LP)
        eperp=eperp/(del*LP*LP)
        gam=del*LP*gam
        eta=eta/LP
        xiu=xiu*L/N/LP
        xir=L/LP/N
        DT=0.5*xiu/(eperp*gam**2.)

      endif

      para(1)=eb
      para(2)=epar
      para(3)=eperp
      para(4)=gam
      para(5)=eta
      para(6)=xir
      para(7)=xiu
      para(8)=lbox
      para(9)=rend
      para(10)=del

      RETURN
      end

!---------------------------------------------------------------*
