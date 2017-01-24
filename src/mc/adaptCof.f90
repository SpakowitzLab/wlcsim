Subroutine adaptCof(downSuccess,nPTReplicas,cof,N_average,&
                     lowerRepExe,upperRepExe,&
                     lowerCofRail,upperCofRail,&
                     repAnnealSpeed,replicaBounds)
use params, only: dp
implicit none

! inputs
integer, intent(in) :: nPTReplicas
integer, intent(in) :: downSuccess(nPTReplicas)
integer, intent(in) :: N_average
double precision, intent(in) :: lowerRepExe
double precision, intent(in) :: upperRepExe
double precision, intent(in) :: lowerCofRail
double precision, intent(in) :: upperCofRail
double precision, intent(in) :: repAnnealSpeed
logical, intent(in) :: replicaBounds ! output between 0 and 1
! input/output
double precision, intent(inout) :: cof(nPTReplicas)

! internal variables
double precision, allocatable :: newCof(:)
integer rep
double precision successRate

allocate( newCof(1:nPTReplicas) )

if(.false.) then
    if ((downSuccess(2).lt.0.99_dp).and.(newCof(1).lt.1.0_dp)) then
        newCof(1)=Cof(1)+repAnnealSpeed
    else
        newCof(1)=Cof(1)
    endif
else
    newCof(1)=Cof(1)
endif
do rep=2,nPTReplicas
    successRate=dble(downSuccess(rep))/dble(N_average)
    if (Cof(rep).lt.Cof(rep-1)) then
        print*, "Error in adaptCof!"
        stop 1
    endif
    ! addapt spacing
    if (successRate.lt.lowerRepExe) then
        newCof(rep)=newCof(rep-1)+(Cof(rep)-Cof(rep-1))*0.95_dp
    elseif (successRate.gt.upperRepExe) then
        newCof(rep)=newCof(rep-1)+(Cof(rep)-Cof(rep-1))*1.05_dp
    else
        newCof(rep)=newCof(rep-1)+Cof(rep)-Cof(rep-1)
    endif
    ! don't addapt too fast.  a.k.a. anneal
    if (newCof(rep).gt.repAnnealSpeed+Cof(rep)) then
        newCof(rep)=Cof(rep)+repAnnealSpeed
    elseif(newCof(rep).lt.Cof(rep)-repAnnealSpeed) then
        newCof(rep)=Cof(rep)-repAnnealSpeed
    endif

    ! inforce limits on addaption
    if ((newCof(rep)-newCof(rep-1)).lt.lowerCofRail) then
        newCof(rep)=newCof(rep-1)+lowerCofRail
    endif
    if ((newCof(rep)-newCof(rep-1)).gt.upperCofRail) then
        newCof(rep)=newCof(rep-1)+upperCofRail
    endif
    if (newCof(rep).lt.newCof(rep-1)) then
        print*, "Error in adaptCof!!!"
        stop 1
    endif

    if (replicaBounds) then
        if(newCof(rep).lt.0.0_dp) newCof(rep)=0.0_dp
        if(newCof(rep).gt.1.0_dp) newCof(rep)=1.0_dp
    endif
enddo
do rep=1,nPTReplicas
    if (newcof(rep).lt.0.00001) then
        print*, "Error in adaptCof! how did that happen?"
        print*, newCof
        print*, Cof
        stop 1
    endif
    cof(rep)=newcof(rep)
enddo
deallocate (newCof)
end subroutine
