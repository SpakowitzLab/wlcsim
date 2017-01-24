subroutine save_repHistory(upSuccess,downSuccess,nPTReplicas, &
                           cofMtrx,xMtrx,nodeNumber,N_average,nExchange,ind,nTerms,s)
    use params, only: dp, MAXFILENAMELEN, outFileUnit
! print Energy data
    IMPLICIT NONE
    integer, intent(in) :: nTerms
    integer, intent(in) :: nPTReplicas
    integer, intent(in) :: upSuccess(nPTReplicas)
    integer, intent(in) :: downSuccess(nPTReplicas)
    integer, intent(in) :: nodeNumber(nPTReplicas)
    integer, intent(in) :: N_average
    integer, intent(in) :: nExchange
    integer, intent(in) :: ind
    real(dp), intent(in) :: cofMtrx(nPTReplicas,nTerms)
    real(dp), intent(in) :: xMtrx(nPTReplicas,nTerms)
    real(dp), intent(in) :: s(nPTReplicas)
    integer rep
    logical isfile
    character(MAXFILENAMELEN) fullName
    fullName=  'data/repHistory'
    inquire(file = fullName, exist=isfile)
    if (isfile) then
        open (unit=outFileUnit, file = fullName, status ='OLD', POSITION="append")
    else
        open (unit=outFileUnit, file = fullName, status = 'new')
    endif

    write(outFileUnit,*) "~~~~~~~~~~~exchange: ",nExchange,", ind:",ind,"~~~~~~~~~~~~~~~~~~~~"
    write(outFileUnit,*) " rep | node|  up  | down |",&
               " chi  |  x_chi |",&
               " h_A  |  x_h_A |",&
               " mu   |  x_mu  |  s   |"
    do rep=1,nPTReplicas
        write(outFileUnit,"(2I6,2f7.4,f7.4,f9.1,f7.4,f9.1,f7.4,f9.1,f7.4)")&
                 rep, nodeNumber(rep), &
                 real(upSuccess(rep))/real(N_average), &
                 real(downSuccess(rep))/real(N_average),&
                 cofMtrx(rep,1), xMtrx(rep,1),&
                 cofMtrx(rep,3), xMtrx(rep,3),&
                 cofMtrx(rep,2), xMtrx(rep,2),s(rep)
    enddo
    close(outFileUnit)
    ! save which node each replica is running on
    fullName=  'data/nodeNumber'
    inquire(file = fullName, exist=isfile)
    if (isfile) then
        open (unit=outFileUnit, file = fullName, status ='OLD', POSITION="append")
    else
        open (unit=outFileUnit, file = fullName, status = 'new')
    endif
    write(outFileUnit,*) ind, nodeNumber
    close(outFileUnit)
    ! save chi values
    fullName=  'data/chi'
    inquire(file = fullName, exist=isfile)
    if (isfile) then
        open (unit=outFileUnit, file = fullName, status ='OLD', POSITION="append")
    else
        open (unit=outFileUnit, file = fullName, status = 'new')
    endif
    write(outFileUnit,*) ind, cofMtrx(:,1)
    close(outFileUnit)
    ! save field strength
    fullName=  'data/h_A'
    inquire(file = fullName, exist=isfile)
    if (isfile) then
        open (unit=outFileUnit, file = fullName, status ='OLD', POSITION="append")
    else
        open (unit=outFileUnit, file = fullName, status = 'new')
    endif
    write(outFileUnit,*) ind, cofMtrx(:,3)
    close(outFileUnit)
    ! save kap
    fullName=  'data/kap'
    inquire(file = fullName, exist=isfile)
    if (isfile) then
        open (unit=outFileUnit, file = fullName, status ='OLD', POSITION="append")
    else
        open (unit=outFileUnit, file = fullName, status = 'new')
    endif
    write(outFileUnit,*) ind, cofMtrx(:,5)
    close(outFileUnit)
    ! save mu
    fullName=  'data/mu'
    inquire(file = fullName, exist=isfile)
    if (isfile) then
        open (unit=outFileUnit, file = fullName, status ='OLD', POSITION="append")
    else
        open (unit=outFileUnit, file = fullName, status = 'new')
    endif
    write(outFileUnit,*) ind, cofMtrx(:,2)
    close(outFileUnit)
    ! save s
    fullName=  'data/s'
    inquire(file = fullName, exist=isfile)
    if (isfile) then
        open (unit=outFileUnit, file = fullName, status ='OLD', POSITION="append")
    else
        open (unit=outFileUnit, file = fullName, status = 'new')
    endif
    write(outFileUnit,*) ind, s
    close(outFileUnit)
end subroutine
