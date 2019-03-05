subroutine save_repHistory(upSuccess,downSuccess,nPTReplicas, &
                           cofMtrx,xMtrx,nodeNumber,N_average,nExchange,ind,nTerms,s)
    use params, only: dp, MAXFILENAMELEN, outFileUnit, print_11char_float
    use energies, only: energyOf, NUMBER_OF_ENERGY_TYPES
! print Energy data
    implicit none
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
    integer ii
    fullName=  'data/repHistory'
    inquire(file = fullName, exist = isfile)
    if (isfile) then
        open (unit = outFileUnit, file = fullName, status ='OLD', POSITION = "append")
    else
        open (unit = outFileUnit, file = fullName, status = 'new')
    endif

    write(outFileUnit,*) "~~~~~~~~~~~exchange: ",nExchange,", ind:",ind,"~~~~~~~~~~~~~~~~~~~~"
    write(outFileUnit,"(35A)", advance="no") " rep | node|  up   | down  |  s   |"
    do ii = 1, NUMBER_OF_ENERGY_TYPES
        if (energyOf(ii)%parallel_temper) then
            write(outFileUnit,"(12A)",advance="no")  " c-",energyOf(ii)%name_str, " "
            write(outFileUnit,"(12A)",advance="no")  " x-",energyOf(ii)%name_str, " "
        endif
    enddo
    write(outFileUnit,*) " "

    do rep = 1,nPTReplicas
        write(outFileUnit,"(2I6,2f8.5,f7.4)", advance="no")&
                 rep, nodeNumber(rep), &
                 real(upSuccess(rep),dp)/real(N_average,dp), &
                 real(downSuccess(rep),dp)/real(N_average,dp),s(rep)
        do ii = 1, NUMBER_OF_ENERGY_TYPES
            if (energyOf(ii)%parallel_temper) then
                call print_11char_float(outFileUnit,cofMtrx(rep,ii))
                write(outFileUnit,"(A)",advance="no") " "
                call print_11char_float(outFileUnit,xMtrx(rep,ii))
                write(outFileUnit,"(A)",advance="no") " "
            endif
        enddo
        write(outFileUnit,*) " "
    enddo
    close(outFileUnit)

    ! save which node each replica is running on
    fullName=  'data/nodeNumber'
    inquire(file = fullName, exist = isfile)
    if (isfile) then
        open (unit = outFileUnit, file = fullName, status ='OLD', POSITION = "append")
    else
        open (unit = outFileUnit, file = fullName, status = 'new')
    endif
    write(outFileUnit,*) ind, nodeNumber
    close(outFileUnit)

    ! save cofs that are subject to pt
    do ii= 1, NUMBER_OF_ENERGY_TYPES
        if (.not. energyOf(ii)%parallel_temper) cycle
        fullName = 'data/' // trim(energyOf(ii)%name_str)
        inquire(file = fullName, exist = isfile)
        if (isfile) then
            open (unit = outFileUnit, file = fullName, status ='OLD', POSITION = "append")
        else
            open (unit = outFileUnit, file = fullName, status = 'new')
        endif
        write(outFileUnit,*) ind, cofMtrx(:,ii)
        close(outFileUnit)
    enddo
    ! save s
    fullName=  'data/s'
    inquire(file = fullName, exist = isfile)
    if (isfile) then
        open (unit = outFileUnit, file = fullName, status ='OLD', POSITION = "append")
    else
        open (unit = outFileUnit, file = fullName, status = 'new')
    endif
    write(outFileUnit,*) ind, s
    close(outFileUnit)
end subroutine
