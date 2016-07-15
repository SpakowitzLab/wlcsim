! *---------------------------------------------------------------*
!
!     subroutine getpara.f90
!     Setup the parameters for the simulation
!
!     1. Determine the simulation type
!     2. Evaluate the polymer elastic parameters
!     3. Determine the parameters for Brownian dynamics simulation
!
!     Andrew Spakowitz
!     8/17/15
!

subroutine getpara(para, dt, simtype)

    use, intrinsic :: iso_fortran_env, dp=>real64
    use sim_params, only : p_lp=>lp, p_l=>l, lbox, p_rend=>rend, vhc, n

    implicit none

    real(dp) para(10)
    real(dp) del
    real(dp) pvec(679,8)
    integer ind,crs
    real(dp) eb,epar,eperp
    real(dp) gam,eta
    real(dp) xir,xiu
    real(dp) m
    real(dp) dt
    integer i
    integer simtype           ! simulation method (wlc=1,sswlc=2,gc=3)
    real(dp) :: l, lp, rend


!     load in the parameters for the simulation

    l=p_l/p_lp
    lp=1.0d0
    del=l/(n-1.0d0)
    rend=p_rend*l

!     load the tabulated parameters

    open (unit=5,file='input/dssWLCparams',status='old')
    do i=1,679
        read(5,*) pvec(i,1),pvec(i,2),pvec(i,3),pvec(i,4),pvec(i,5),pvec(i,6),pvec(i,7),pvec(i,8)
    end do
    close(5)


!     setup the parameters for wlc simulation

    if (del.lt.pvec(1,1)) then
        eb=lp/del
        gam=del
        xir=l/n
        simtype=1
    endif

!    setup the parameters for gc simulation

    if (del.gt.pvec(679,1)) then
        epar=1.5/(del*lp**2.)
        gam=0.
        simtype=3
        xir=l/n
    endif

!    setup the parameters for sswlc simulation

    if (del.ge.pvec(1,1).and.del.le.pvec(679,1)) then
        simtype=2

        crs=0
        ind=1
        do while (crs.eq.0)
        if (del.le.pvec(ind,1)) then
            crs=1
        else
            ind=ind+1
        endif
        enddo

        i=2
        m=(pvec(ind,i)-pvec(ind-1,i))/(pvec(ind,1)-pvec(ind-1,1))
        eb=m*(del-pvec(ind,1))+pvec(ind,i)

        i=3
        m=(pvec(ind,i)-pvec(ind-1,i))/(pvec(ind,1)-pvec(ind-1,1))
        gam=m*(del-pvec(ind,1))+pvec(ind,i)

        i=4
        m=(pvec(ind,i)-pvec(ind-1,i))/(pvec(ind,1)-pvec(ind-1,1))
        epar=m*(del-pvec(ind,1))+pvec(ind,i)

        i=5
        m=(pvec(ind,i)-pvec(ind-1,i))/(pvec(ind,1)-pvec(ind-1,1))
        eperp=m*(del-pvec(ind,1))+pvec(ind,i)

        i=6
        m=(pvec(ind,i)-pvec(ind-1,i))/(pvec(ind,1)-pvec(ind-1,1))
        eta=m*(del-pvec(ind,1))+pvec(ind,i)

        i=7
        m=(pvec(ind,i)-pvec(ind-1,i))/(pvec(ind,1)-pvec(ind-1,1))
        xiu=m*(del-pvec(ind,1))+pvec(ind,i)

!         i=8
!         m=(pvec(ind,i)-pvec(ind-1,i))/(pvec(ind,1)-pvec(ind-1,1))
!         dt=xiu*(m*(del-pvec(ind,1))+pvec(ind,i))

        eb=eb/del
        epar=epar/del
        eperp=eperp/del
        gam=del*gam

        xiu=xiu*l/n
        xir=l/n
        dt=0.5*xiu/(eperp*gam**2.)

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

    return
end

!---------------------------------------------------------------*
