!! ---------------------------------------------------------------*
!
!     This subroutine sets the initial positions (and orientations,
!     if applicable) of the chain.
!
subroutine initcond(R,U,AB,NT,N,NP,FRMfile,para,lbox, &
                    settype,rand_stat)

!use mt19937, only : grnd, init_genrand, rnorm, mt, mti
use mersenne_twister
use setPrecision

IMPLICIT NONE

real(dp) PI
paraMETER (PI=3.141593_dp)

real(dp) R(NT,3)  ! Bead positions
real(dp) U(NT,3)  ! Tangent vectors
integer AB(NT)            ! Chemical identity of beads
integer N,NP,NT           ! Number of beads
real(dp) gam      ! Equil bead separation
real(dp) lbox(3)  ! Box edge length
integer I,J,IB            ! index Holders
logical FRMfile           ! Is conformation in file?
!real(dp) RMIN
real(dp) R0(3)
real(dp) para(10)
integer settype           ! select what type of configurateion

!     Varibles for type 2

real(dp) Uold(3) ! save previous direction
real(dp) Rold(3) ! save previous position
real(dp) theta   ! random angle
real(dp) z       ! random z position
real(dp) rr      ! random radial position
logical search           ! for exiting while loop
real(dp) test(3) ! test position for inside confinment
real(dp) Rc      ! radius of confinement
integer ii !for testing


!      Random number generator initiation
type(random_stat) rand_stat
real urand(3)
real(dp) mag  ! magnitude of U for reload

!     Setup the choice parameters

if(FRMfile)then
   open (unit = 5, file = 'input/r0', status = 'OLD')
   do I=1,NT
      read(5,*) R(I,1),R(I,2),R(I,3),AB(I)
   enddo
   close(5)

   open (unit = 5, file = 'input/u0', status = 'OLD')
   do I=1,NT
      read(5,*) U(I,1),U(I,2),U(I,3)
       mag=sqrt(U(I,1)**2+U(I,2)**2+U(I,3)**2)
       U(I,1)=U(I,1)/mag
       U(I,2)=U(I,2)/mag
       U(I,3)=U(I,3)/mag
   enddo

   close(5)
   return
endif


!     Fix the initial condition
if(settype.eq.1) then
! staight line in y direction with random starting position
    gam=para(4)

    IB=1
    do I=1,NP
       call random_number(urand,rand_stat)
       R0(1)=urand(1)*lbox(1)
       R0(2)=urand(2)*lbox(2)
       R0(3)=urand(3)*lbox(3)
       do J=1,N
          R(IB,1)=R0(1)
          R(IB,2)=R0(2)+gam*(J-N/2.0_dp-0.5_dp) ! center on box
          R(IB,3)=R0(3)
          U(IB,1)=0.
          U(IB,2)=1.
          U(IB,3)=0.
          IB=IB+1
       enddo
    enddo
else if(settype.eq.2) then
    ! travel in radom direction
    ! rerandomize when reaching boundary
    ! slit boundary in z direction
    gam=para(4)

    IB=1
    do  I=1,NP
        call random_number(urand,rand_stat)
        Rold(1)=urand(1)*lbox(1)
        Rold(2)=urand(2)*lbox(2)
        Rold(3)=urand(3)*lbox(3)
        call random_number(urand,rand_stat)
        theta=urand(1)*2*PI
        z=urand(2)*2-1
        Uold(1)=sqrt(1-z*z)*cos(theta)
        Uold(2)=sqrt(1-z*z)*sin(theta)
        Uold(3)=z

        do J=1,N
           search=.TRUE.
           ii=0
           do while(search)
                ii=ii+1
                if(ii.gt.100) then
                    print*,'stuck in loop'
                    print*,'Rold=',Rold(1),Rold(2),Rold(3)
                    print*,'test=',test(1),test(2),test(3)
                    exit
                endif
                test(1)=Rold(1)+Uold(1)*gam
                test(2)=Rold(2)+Uold(2)*gam
                test(3)=Rold(3)+Uold(3)*gam
                search=.FALSE.
                if(test(3).gt.lbox(3))then
                    search=.TRUE.
                endif
                if(test(3).lt.0)then
                    search=.TRUE.
                endif
                if(search) then
                     call random_number(urand,rand_stat)
                     theta=urand(1)*2*PI
                     z=urand(2)*2-1
                     Uold(1)=sqrt(1-z*z)*cos(theta)
                     Uold(2)=sqrt(1-z*z)*sin(theta)
                     Uold(3)=z
                endif
           enddo
           R(IB,1)=test(1)
           R(IB,2)=test(2)
           R(IB,3)=test(3)
           Rold(1)=test(1)
           Rold(2)=test(2)
           Rold(3)=test(3)
           U(IB,1)=Uold(1)
           U(IB,2)=Uold(2)
           U(IB,3)=Uold(3)

           IB=IB+1
        enddo
    enddo
else if(settype.eq.3) then
    ! travel in radom direction
    ! rerandomize when reaching boundary
    ! square boundary
    gam=para(4)

    IB=1
    do  I=1,NP
       call random_number(urand,rand_stat)
       Rold(1)=urand(1)*lbox(1)
       Rold(2)=urand(2)*lbox(2)
       Rold(3)=urand(3)*lbox(3)
       call random_number(urand,rand_stat)
       theta=urand(1)*2*PI
       z=urand(2)*2-1
       Uold(1)=sqrt(1-z*z)*cos(theta)
       Uold(2)=sqrt(1-z*z)*sin(theta)
       Uold(3)=z

       do J=1,N
          search=.TRUE.
          ii=0
          do while(search)
               ii=ii+1
               if(ii.gt.100) then
                   print*,'stuck in loop'
                   print*,'Rold=',Rold(1),Rold(2),Rold(3)
                   print*,'test=',test(1),test(2),test(3)
                   exit
               endif
               test(1)=Rold(1)+Uold(1)*gam
               test(2)=Rold(2)+Uold(2)*gam
               test(3)=Rold(3)+Uold(3)*gam
               search=.FALSE.
               if(test(1).gt.lbox(1))then
                   search=.TRUE.
               endif
               if(test(1).lt.0)then
                   search=.TRUE.
               endif
               if(test(2).gt.lbox(2))then
                   search=.TRUE.
               endif
               if(test(2).lt.0)then
                   search=.TRUE.
               endif
               if(test(3).gt.lbox(3))then
                   search=.TRUE.
               endif
               if(test(3).lt.0)then
                   search=.TRUE.
               endif
               if(search) then
                    call random_number(urand,rand_stat)
                    theta=urand(1)*2_dp*PI
                    z=urand(2)*2.0_dp-1.0_dp
                    Uold(1)=sqrt(1-z*z)*cos(theta)
                    Uold(2)=sqrt(1-z*z)*sin(theta)
                    Uold(3)=z
               endif
          enddo
          R(IB,1)=test(1)
          R(IB,2)=test(2)
          R(IB,3)=test(3)
          Rold(1)=test(1)
          Rold(2)=test(2)
          Rold(3)=test(3)
          U(IB,1)=Uold(1)
          U(IB,2)=Uold(2)
          U(IB,3)=Uold(3)

          IB=IB+1
       enddo
    enddo

else if(settype.eq.4) then
    ! travel in radom direction
    ! rerandomize when reaching boundary
    ! shpere boundary
    ! radius of lbox/2 centered at lbox/2
    Rc=lbox(1)/2.0_dp ! use lbox as radius
       gam=para(4)
    IB=1
    do  I=1,NP
       call random_number(urand,rand_stat)
       theta=urand(1)*2.0_dp*PI
       z=urand(2)*2.0_dp-1.0_dp
       rr=Rc*urand(3)  ! should have an r**2 from jacobian
       Rold(1)=sqrt(1.0_dp-z*z)*cos(theta)*rr + Rc
       Rold(2)=sqrt(1.0_dp-z*z)*sin(theta)*rr + Rc
       Rold(3)=z*rr + Rc
       call random_number(urand,rand_stat)
       theta=urand(1)*2_dp*PI
       z=urand(2)*2.0_dp-1.0_dp
       Uold(1)=sqrt(1-z*z)*cos(theta)
       Uold(2)=sqrt(1-z*z)*sin(theta)
       Uold(3)=z
       do J=1,N
           search=.TRUE.
           do while(search)
               test(1)=Rold(1)+Uold(1)*gam
               test(2)=Rold(2)+Uold(2)*gam
               test(3)=Rold(3)+Uold(3)*gam
               search=.FALSE.
               if((test(1)-Rc)**2+&
                   (test(2)-Rc)**2+&
                   (test(3)-Rc)**2.gt.Rc**2)then
                   search=.TRUE.
               endif
               if(search) then
                    call random_number(urand,rand_stat)
                    theta=urand(1)*2.0_dp*PI
                    z=urand(2)*2.0_dp-1.0_dp
                    Uold(1)=sqrt(1.0_dp-z*z)*cos(theta)
                    Uold(2)=sqrt(1.0_dp-z*z)*sin(theta)
                    Uold(3)=z
               endif
           enddo
           R(IB,1)=test(1)
           R(IB,2)=test(2)
           R(IB,3)=test(3)
           Rold(1)=test(1)
           Rold(2)=test(2)
           Rold(3)=test(3)
           U(IB,1)=Uold(1)
           U(IB,2)=Uold(2)
           U(IB,3)=Uold(3)
           IB=IB+1
       enddo ! loop to N
    enddo ! loop to np
else if(settype.eq.5) then
    ! randomly distribute beads in shereical confinement
    do IB=1,NT
        search=.true.
        do while(search)
             call random_number(urand,rand_stat)
             test(1)=urand(1)*lbox(1)
             test(2)=urand(2)*lbox(2)
             test(3)=urand(3)*lbox(3)
             if(((test(1)-lbox(1)/2.0_dp)**2+ &
                 (test(2)-lbox(1)/2.0_dp)**2+ &
                 (test(3)-lbox(1)/2.0_dp)**2).lt.(lbox(1)*lbox(1)*0.25_dp)) then
                 search=.false.
             endif
        enddo
        R(IB,1)=test(1)
        R(IB,2)=test(2)
        R(IB,3)=test(3)
        U(IB,1)=0.00_dp
        U(IB,2)=0.00_dp
        U(IB,3)=0.00_dp
    enddo
else if (settype == 6) then
    IB=1
    do  I=1,NP
        call random_number(urand,rand_stat)
        R0(1)=urand(1)*lbox(1)
        call random_number(urand,rand_stat)
        R0(2)=urand(1)*lbox(1)
        call random_number(urand,rand_stat)
        R0(3)=urand(1)*lbox(1)
        do  J=1,NB
            if (ring.EQ.0) then
                R(IB,1)=R0(1)
                R(IB,2)=R0(2)+gam*(J-NB/2.0_dp-0.5_dp)
                R(IB,3)=R0(3)
                U(IB,1)=0.0_dp
                U(IB,2)=1.0_dp
                U(IB,3)=0.0_dp
            else
                R(IB,1)=R0(1)+((gam*NB)/(2*PI))*Cos(J*2.0_dp*PI/NB)
                R(IB,2)=R0(2)+((gam*NB)/(2*PI))*Sin(J*2.0_dp*PI/NB)
                R(IB,3)=0.0_dp
                U(IB,1)=-Sin(J*2.0_dp*PI/NB)
                U(IB,2)=Cos(J*2.0_dp*PI/NB)
                U(IB,3)=0.0_dp;
            endif
        IB=IB+1
        enddo
    enddo

endif


RETURN
end

!---------------------------------------------------------------*
