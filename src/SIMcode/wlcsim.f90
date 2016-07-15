!---------------------------------------------------------------*

program wlcsim

!
!     wlc simulation package:
!     simulation package for brownian dynamics and
!     monte carlo simulation
!
!     andrew spakowitz
!     version 1.0
!     8/17/2015
!

!     variables within the simulation

    use mt19937, only : grnd, sgrnd, rnorm, mt, mti
    use sim_params, only : n, np, dp, tf, indmax, dt, frmfile, brown, &
                           inton, logtime, ninit, nstep, fpt_dist, col_type
    implicit none

    real(dp), parameter :: pi=3.141592654 ! value of pi

    integer, parameter :: nt = n*np ! total # beads among all polymers in sim
    real(dp) :: r(nt,3)     ! conformation of polymer chains
    real(dp) :: u(nt,3)     ! conformation of polymer chains
    real(dp) :: r0(nt,3)    ! conformation of polymer chains
    real(dp) :: u0(nt,3)    ! conformation of polymer chains

    real(dp) l0             ! equilibrium segment length
    real(dp) energy         ! total energy
    real(dp) time           ! current time
    real(dp) tsave          ! time of save point
    integer i,j,ib          ! index
    integer ind             ! ind in series
    integer tens            ! decimal of index
    character*4 fileind     ! index of output
    character*16 snapnm     ! file for output

!     simulation input variables

    real(dp) :: dt0 = dt      ! initial time step size

!     monte carlo variables

    real(dp) mcamp(6)   ! amplitude of random change
    integer moveon(6)   ! is the move active
    integer window(6)   ! size of window for bead selection
    integer success(6)  ! number of successes

!     energy variables

    real(dp) eelas(3)   ! elastic energy
    real(dp) eponp      ! poly-poly energy

!     structure analysis

    real(dp) rcom(3)    ! center of mass
    real(dp) delr(3)    ! mag of gyration tensor
    real(dp) rcom0(3)   ! init val rcom
    real(dp) delr0(3)   ! init val delr
    real(dp) drcom      ! change in rcom
    real(dp) sig(3,3)
    real(dp) cor

!     variables in the simulation

    real(dp) para(10)
    integer simtype     ! simulation method (wlc=1,sswlc=2,gc=3)

!     variables for the random number generators

    integer idum        ! seed for the generator
    real(dp) mom(6)

!     variable to hold time of first collisions between each bead
    real(dp), allocatable, dimension(:,:):: has_collided

    call getpara(para,dt0,simtype)

    if (col_type.ne.0) then
        allocate(has_collided(nt,nt))
        has_collided = -1.0d+0
    endif

!     setup the initial condition

    call initcond(r,u,nt,n,np,idum,frmfile,para)

!     turn on moves for each simulation type

    if (simtype.eq.1) then
        mcamp(1)=1.
        mcamp(2)=1.
        mcamp(3)=1.
        mcamp(4)=1.
        mcamp(5)=1.
        mcamp(6)=1.
        moveon(1)=1
        moveon(2)=0
        moveon(3)=1
        moveon(4)=0
    elseif (simtype.eq.2) then
        mcamp(1)=1.
        mcamp(2)=1.
        mcamp(3)=1.
        mcamp(4)=1.
        mcamp(5)=1.
        mcamp(6)=1.
        moveon(1)=1
        moveon(2)=1
        moveon(3)=1
        moveon(4)=1
    elseif (simtype.eq.3) then
        mcamp(1)=1.
        mcamp(2)=1.
        mcamp(3)=1.
        mcamp(4)=1.
        mcamp(5)=1.
        mcamp(6)=1.
        moveon(1)=1
        moveon(2)=1
        moveon(3)=1
        moveon(4)=0
    endif

!     turn off whole chain rotation and translation if interactions are off

    if (inton.eq.1) then
        moveon(5)=1
        moveon(6)=1
    else
        moveon(5)=0
        moveon(6)=0
    endif

!     perform an initialization mc simulation

    call mcsim(r,u,nt,n,np,ninit,brown,inton,idum,para,mcamp, &
        success,moveon,window,simtype)

!     save the conformation and psi angles

    open (unit = 1, file = 'data/r0', status = 'new')
    ib=1
    do i=1,np
        do j=1,n
            r0(ib,1)=r(ib,1)
            r0(ib,2)=r(ib,2)
            r0(ib,3)=r(ib,3)
            u0(ib,1)=u(ib,1)
            u0(ib,2)=u(ib,2)
            u0(ib,3)=u(ib,3)
            write(1,*) r(ib,1),r(ib,2),r(ib,3)
            ib=ib+1
        end do
    end do
    close(1)

    open (unit = 1, file = 'data/u0', status = 'new')
    ib=1
    do i=1,np
        do j=1,n
        write(1,*) u(ib,1),u(ib,2),u(ib,3)
        ib=ib+1
        end do
    end do
    close(1)

!     begin simulation

    ind=1
    time=0.

!     open the output files

    open (unit = 2, file = 'data/out1', status = 'new')
    open (unit = 3, file = 'data/out2', status = 'new')
    open (unit = 4, file = 'data/out3', status = 'new')

    call stress(sig,r,u,nt,n,np,para,inton)

    write(3,*) real(sig(1,1)),real(sig(1,2)),real(sig(1,3)),real(sig(2,1)),real(sig(2,2))
    write(4,*) real(sig(2,3)),real(sig(3,1)),real(sig(3,2)),real(sig(3,3))

    do while (ind.le.indmax)

!     perform a mc simulation, only if nstep.ne.0

        call mcsim(r,u,nt,n,np,nstep,brown,inton,idum,para,mcamp, &
            success,moveon,window,simtype)

!     perform a brownian dynamics simulation over time step

        if (logtime.eq.0) then
        tsave = tf*ind/indmax
        else
        tsave = dt0*exp((ind-1.)/(indmax-1.)*log(tf/dt0))
        endif
        if (nstep.eq.0) then
        call bdsim(r,u,nt,n,np,time,tsave,dt0,brown,inton,idum, &
                    para,simtype,has_collided,fpt_dist,col_type)
        endif

!     save the conformation and the metrics

        tens=nint(log10(1.*ind)-0.4999)+1
        write (fileind,'(i4)'), ind
        snapnm= 'data/r'//fileind((4-tens+1):4)
        open (unit = 1, file = snapnm, status = 'new')
        ib=1
        do i=1,np
            do j=1,n
                write(1,*) r(ib,1),r(ib,2),r(ib,3)
                ib=ib+1
            end do
        end do
        close(1)

        snapnm= 'data/u'//fileind((4-tens+1):4)
        open (unit = 1, file = snapnm, status = 'new')
        ib=1
        do i=1,np
            do j=1,n
                write(1,*) u(ib,1),u(ib,2),u(ib,3)
                ib=ib+1
            end do
        end do
        close(1)

        snapnm='data/coltimes'
        open (unit=1, file=snapnm, status='replace')
        do, i=1,nt
            write(1,*) ( has_collided(i,j), j=1,nt )
        enddo
        close(1)

        call stress(sig,r,u,nt,n,np,para,inton,simtype)
        call stressp(cor,r,u,r0,u0,nt,n,np,para,inton,simtype)

        call energy_elas(eelas,r,u,nt,n,np,para)
        eponp=0.
        if (inton.eq.1) then
        call energy_ponp(eponp,r,nt,n,np,para)
        endif
        write(2,*) real(time),real(eelas(1)),real(eelas(2)),real(eelas(3)),real(eponp),real(cor)


        write(3,*) real(sig(1,1)),real(sig(1,2)),real(sig(1,3)),real(sig(2,1)),real(sig(2,2))
        write(4,*) real(sig(2,3)),real(sig(3,1)),real(sig(3,2)),real(sig(3,3))

        print*, '________________________________________'
        print*, 'time point ',ind, ' out of', indmax
        print*, 'current time ',time
        print*, 'bending energy ', eelas(1)
        print*, 'par compression energy ', eelas(2)
        print*, 'perp compression energy ', eelas(3)
        print*, 'polymer-polymer energy ', eponp
        print*, 'current number of beads ', n
        print*, 'time step ', dt
        print*, 'end-to-end distance poly 1 ', norm2(r(n,:) - r(1,:))
        print*, 'simulation type ', simtype

        ind=ind+1

    enddo
    end


!---------------------------------------------------------------*

