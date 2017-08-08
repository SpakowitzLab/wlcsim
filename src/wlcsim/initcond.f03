!! ---------------------------------------------------------------*

!
!     This subroutine sets the initial condition for a chain within
!     the capsid
!
!     Andrew Spakowitz
!     Written 4-16-04
!
!     Updated by Quinn in 2016
!
subroutine initcond(R,U,NT,NB,NP,FRMFILE,PARA,LBOX, &
                    initCondType,rand_stat,ring, wlc_p)

!use mt19937, only : grnd, init_genrand, rnorm, mt, mti
use mersenne_twister
use params, only: dp, pi, wlcsim_params
use inputparams, only: MAXPARAMLEN

implicit none

type(wlcsim_params), intent(in) :: wlc_p

logical ring
integer NB,NP,NT           ! Number of beads
real(dp) R(3,NT)  ! Bead positions
real(dp) U(3,NT)  ! Tangent vectors
real(dp) GAM      ! Equil bead separation
real(dp) LBOX(3)  ! Box edge length
integer I,J,IB            ! Index Holders
LOGICAL FRMFILE           ! Is conformation in file?
!real(dp) RMin
real(dp) R0(3)
real(dp) PARA(10)
character(MAXPARAMLEN) initCondType           ! select what type of configurateion

!     Varibles for type 2

real(dp) Uold(3) ! save previous direction
real(dp) Rold(3) ! save previous position
real(dp) theta   ! random angle
real(dp) z       ! random z position
real(dp) rr      ! random radial position
LOGICAL search           ! for exiting while loop
real(dp) test(3) ! test position for inside confinment
real(dp) Rc      ! radius of confinement
integer ii !for testing

!   temporary variables
real(dp) mag    ! magnitude of U for reload, or of U when smoothing

!      Random number generator initiation
type(random_stat) rand_stat
real urand(3)
real nrand(3)

GAM = PARA(4)

!     Setup the choice parameters

if (FRMFILE)then
   open (UNIT = 5, FILE = 'input/r0', STATUS = 'OLD')
   Do I = 1,NT
      READ(5,*) R(1,I),R(2,I),R(3,I)
   enddo
   CLOSE(5)

   open (UNIT = 5, FILE = 'input/u0', STATUS = 'OLD')
   Do I = 1,NT
      READ(5,*) U(1,I),U(2,I),U(3,I)
       mag = sqrt(U(1,I)**2 + U(2,I)**2 + U(3,I)**2)
       U(1,I) = U(1,I)/mag
       U(2,I) = U(2,I)/mag
       U(3,I) = U(3,I)/mag
   enddo

   CLOSE(5)
   return
endif

!     Fix the initial condition
if (initCondType.eq.'lineInYFromOrigin') then
! straight line in y direction starting at origin, equilibrium bond lengths
    iB = 1
    do i = 1,nP
        do j = 1,nB
            R(1,iB) = 0.0_dp
            R(2,iB) = gam*j
            R(3,iB) = 0.0_dp
            U(1,iB) = 0.0_dp
            U(2,iB) = 1.0_dp
            U(3,iB) = 0.0_dp
            iB = iB + 1
        enddo
    enddo

else if (initCondType.eq.'lineInY') then
! staight line in y direction with random starting position

    IB = 1
    do I = 1,NP
       call random_number(urand,rand_stat)
       R0(1) = urand(1)*LBOX(1)
       R0(2) = urand(2)*LBOX(2)
       R0(3) = urand(3)*LBOX(3)

       do J = 1,NB
          R(1,IB) = R0(1)
          R(2,IB) = R0(2) + GAM*(J - NB/2.0_dp - 0.5_dp) ! center on box
          R(3,IB) = R0(3)
          U(1,IB) = 0.
          U(2,IB) = 1.
          U(3,IB) = 0.
          IB = IB + 1

       enddo
    enddo

else if (initCondType.eq.'randomLineSlitInZBoundary') then
    ! travel in radom direction
    ! rerandomize when reaching boundary
    ! slit boundary in z direction

    IB = 1
    do  I = 1,NP
        call random_number(urand,rand_stat)
        Rold(1) = urand(1)*LBOX(1)
        Rold(2) = urand(2)*LBOX(2)
        Rold(3) = urand(3)*LBOX(3)
        call random_number(urand,rand_stat)
        theta = urand(1)*2*PI
        z = urand(2)*2-1
        Uold(1) = sqrt(1-z*z)*cos(theta)
        Uold(2) = sqrt(1-z*z)*sin(theta)
        Uold(3) = z

        do J = 1,NB
           search = .TRUE.
           ii = 0
           do while(search)
                ii = ii + 1
                if (ii.gt.100) then
                    print*,'stuck in loop'
                    print*,'Rold = ',Rold(1),Rold(2),Rold(3)
                    print*,'test = ',test(1),test(2),test(3)
                    exit
                endif
                test(1) = Rold(1) + Uold(1)*GAM
                test(2) = Rold(2) + Uold(2)*GAM
                test(3) = Rold(3) + Uold(3)*GAM
                search = .FALSE.
                if (test(3).gt.LBOX(3))then
                    search = .TRUE.
                endif
                if (test(3).lt.0)then
                    search = .TRUE.
                endif
                if (search) then
                     call random_number(urand,rand_stat)
                     theta = urand(1)*2*PI
                     z = urand(2)*2-1
                     Uold(1) = sqrt(1-z*z)*cos(theta)
                     Uold(2) = sqrt(1-z*z)*sin(theta)
                     Uold(3) = z
                endif
           enddo
           R(1,IB) = test(1)
           R(2,IB) = test(2)
           R(3,IB) = test(3)
           Rold(1) = test(1)
           Rold(2) = test(2)
           Rold(3) = test(3)
           U(1,IB) = Uold(1)
           U(2,IB) = Uold(2)
           U(3,IB) = Uold(3)

           IB = IB + 1
        enddo
    enddo
else if (initCondType.eq.'randomLineCubeBoundary') then
    ! travel in radom direction
    ! rerandomize when reaching boundary
    ! square boundary

    IB = 1
    do  I = 1,NP
       call random_number(urand,rand_stat)
       Rold(1) = urand(1)*LBOX(1)
       Rold(2) = urand(2)*LBOX(2)
       Rold(3) = urand(3)*LBOX(3)
       call random_number(urand,rand_stat)
       theta = urand(1)*2*PI
       z = urand(2)*2-1
       Uold(1) = sqrt(1-z*z)*cos(theta)
       Uold(2) = sqrt(1-z*z)*sin(theta)
       Uold(3) = z

       do J = 1,NB
          search = .TRUE.
          ii = 0
          do while(search)
               ii = ii + 1
               if (ii.gt.100) then
                   print*,'stuck in loop'
                   print*,'Rold = ',Rold(1),Rold(2),Rold(3)
                   print*,'test = ',test(1),test(2),test(3)
                   exit
               endif
               test(1) = Rold(1) + Uold(1)*GAM
               test(2) = Rold(2) + Uold(2)*GAM
               test(3) = Rold(3) + Uold(3)*GAM
               search = .FALSE.
               if (test(1).gt.LBOX(1))then
                   search = .TRUE.
               endif
               if (test(1).lt.0)then
                   search = .TRUE.
               endif
               if (test(2).gt.LBOX(2))then
                   search = .TRUE.
               endif
               if (test(2).lt.0)then
                   search = .TRUE.
               endif
               if (test(3).gt.LBOX(3))then
                   search = .TRUE.
               endif
               if (test(3).lt.0)then
                   search = .TRUE.
               endif
               if (search) then
                    call random_number(urand,rand_stat)
                    theta = urand(1)*2_dp*PI
                    z = urand(2)*2.0_dp-1.0_dp
                    Uold(1) = sqrt(1-z*z)*cos(theta)
                    Uold(2) = sqrt(1-z*z)*sin(theta)
                    Uold(3) = z
               endif
          enddo
          R(1,IB) = test(1)
          R(2,IB) = test(2)
          R(3,IB) = test(3)
          Rold(1) = test(1)
          Rold(2) = test(2)
          Rold(3) = test(3)
          U(1,IB) = Uold(1)
          U(2,IB) = Uold(2)
          U(3,IB) = Uold(3)

          IB = IB + 1
       enddo
    enddo

else if (initCondType.eq.'randomLineSphereBoundary') then
    ! travel in radom direction
    ! rerandomize when reaching boundary
    ! shpere boundary
    ! radius of LBox/2 centered at LBox/2
    Rc = LBOX(1)/2.0_dp ! use LBOX as radius
    IB = 1
    do  I = 1,NP
       call random_number(urand,rand_stat)
       theta = urand(1)*2.0_dp*PI
       z = urand(2)*2.0_dp-1.0_dp
       rr = Rc*urand(3)  ! should have an r**2 from jacobian
       Rold(1) = sqrt(1.0_dp-z*z)*cos(theta)*rr + Rc
       Rold(2) = sqrt(1.0_dp-z*z)*sin(theta)*rr + Rc
       Rold(3) = z*rr + Rc
       call random_number(urand,rand_stat)
       theta = urand(1)*2_dp*PI
       z = urand(2)*2.0_dp-1.0_dp
       Uold(1) = sqrt(1-z*z)*cos(theta)
       Uold(2) = sqrt(1-z*z)*sin(theta)
       Uold(3) = z
       do J = 1,NB
           search = .TRUE.
           do while(search)
               test(1) = Rold(1) + Uold(1)*GAM
               test(2) = Rold(2) + Uold(2)*GAM
               test(3) = Rold(3) + Uold(3)*GAM
               search = .FALSE.
               if ((test(1)-Rc)**2 + &
                   (test(2)-Rc)**2 + &
                   (test(3)-Rc)**2.gt.Rc**2)then
                   search = .TRUE.
               endif
               if (search) then
                    call random_number(urand,rand_stat)
                    theta = urand(1)*2.0_dp*PI
                    z = urand(2)*2.0_dp-1.0_dp
                    Uold(1) = sqrt(1.0_dp-z*z)*cos(theta)
                    Uold(2) = sqrt(1.0_dp-z*z)*sin(theta)
                    Uold(3) = z
               endif
           enddo
           R(1,IB) = test(1)
           R(2,IB) = test(2)
           R(3,IB) = test(3)
           Rold(1) = test(1)
           Rold(2) = test(2)
           Rold(3) = test(3)
           U(1,IB) = Uold(1)
           U(2,IB) = Uold(2)
           U(3,IB) = Uold(3)
           IB = IB + 1
       enddo ! loop to N
    enddo ! loop to np
else if (initCondType.eq.'randomlyDistributeBeadsInSphere') then
    ! randomly distribute beads in shereical confinement
    do IB = 1,NT
        search = .true.
        do while(search)
             call random_number(urand,rand_stat)
             test(1) = urand(1)*LBox(1)
             test(2) = urand(2)*LBox(2)
             test(3) = urand(3)*LBox(3)
             if (((test(1)-LBox(1)/2.0_dp)**2+ &
                 (test(2)-LBox(1)/2.0_dp)**2+ &
                 (test(3)-LBox(1)/2.0_dp)**2).lt.(LBox(1)*LBox(1)*0.25_dp)) then
                 search = .false.
             endif
        enddo
        R(1,IB) = test(1)
        R(2,IB) = test(2)
        R(3,IB) = test(3)
        U(1,IB) = 0.00_dp
        U(2,IB) = 0.00_dp
        U(3,IB) = 0.00_dp
    enddo
else if (initCondType == 'ring') then
    IB = 1
    do  I = 1,NP
        call random_number(urand,rand_stat)
        R0(1) = urand(1)*LBOX(1)
        call random_number(urand,rand_stat)
        R0(2) = urand(1)*LBOX(1)
        call random_number(urand,rand_stat)
        R0(3) = urand(1)*LBOX(1)
        do  J = 1,NB
             R(1,IB) = R0(1) + ((GAM*NB)/(2*PI))*Cos(J*2.0_dp*PI/NB)
             R(2,IB) = R0(2) + ((GAM*NB)/(2*PI))*Sin(J*2.0_dp*PI/NB)
             R(3,IB) = 0.0_dp
             U(1,IB) = -Sin(J*2.0_dp*PI/NB)
             U(2,IB) = Cos(J*2.0_dp*PI/NB)
             U(3,IB) = 0.0_dp;
             IB = IB + 1
        ENDdo
    ENDdo

else if (initCondType == 'randomWalkWithBoundary') then
    ! initialize as if it were a gaussian chain, first bead at zero
    IB = 1
    do I = 1,NP
        R(1,IB) = 0.0_dp
        R(2,IB) = 0.0_dp
        R(3,IB) = 0.0_dp
        IB = IB + 1
        do J = 2,NB
            call random_gauss(nrand, rand_stat)
            R(1,IB) = R(1,IB-1) + wlc_p%sigma*nrand(1)
            R(2,IB) = R(2,IB-1) + wlc_p%sigma*nrand(2)
            R(3,IB) = R(3,IB-1) + wlc_p%sigma*nrand(3)
            IB = IB + 1
        enddo
    enddo
    ! now initialize the orientation vectors by smoothing
    IB = 1
    do I = 1,NP
        U(1,IB) = R(1,IB + 1) - R(1,IB)
        U(2,IB) = R(2,IB + 1) - R(2,IB)
        U(3,IB) = R(3,IB + 1) - R(3,IB)
        mag = sqrt(U(1,IB)*U(1,IB) + U(2,IB)*U(2,IB) + U(3,IB)*U(3,IB))
        U(1,IB) = U(1,IB)/mag
        U(2,IB) = U(2,IB)/mag
        U(3,IB) = U(3,IB)/mag
        IB = IB + 1
        do J = 2,NB-1
            U(1,IB) = R(1,IB + 1) - R(1,IB-1)
            U(2,IB) = R(2,IB + 1) - R(2,IB-1)
            U(3,IB) = R(3,IB + 1) - R(3,IB-1)
            mag = sqrt(U(1,IB)*U(1,IB) + U(2,IB)*U(2,IB) + U(3,IB)*U(3,IB))
            U(1,IB) = U(1,IB)/mag
            U(2,IB) = U(2,IB)/mag
            U(3,IB) = U(3,IB)/mag
            IB = IB + 1
        enddo
        U(1,IB) = R(1,IB) - R(1,IB-1)
        U(2,IB) = R(2,IB) - R(2,IB-1)
        U(3,IB) = R(3,IB) - R(3,IB-1)
        mag = sqrt(U(1,IB)*U(1,IB) + U(2,IB)*U(2,IB) + U(3,IB)*U(3,IB))
        U(1,IB) = U(1,IB)/mag
        U(2,IB) = U(2,IB)/mag
        U(3,IB) = U(3,IB)/mag
        IB = IB + 1
    enddo
else
    print*, "Unknown version of chain initialization initCondType....."
    stop 1
endif


RETURN
END

!---------------------------------------------------------------*
