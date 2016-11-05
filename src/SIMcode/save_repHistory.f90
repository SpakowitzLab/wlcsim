Subroutine save_repHistory(upSuccess,downSuccess,nPTReplicas, &
                           cofMtrx,xMtrx,nodeNumber,N_average,nExchange,IND,nTerms,s)
    use setPrecision
! Print Energy data
    IMPLICIT NONE
    integer, intent(in) :: nTerms
    integer, intent(in) :: nPTReplicas
    integer, intent(in) :: upSuccess(nPTReplicas)
    integer, intent(in) :: downSuccess(nPTReplicas)
    integer, intent(in) :: nodeNumber(nPTReplicas)
    integer, intent(in) :: N_average
    integer, intent(in) :: nExchange
    integer, intent(in) :: IND
    double precision, intent(in) :: cofMtrx(nPTReplicas,nTerms)
    double precision, intent(in) :: xMtrx(nPTReplicas,nTerms)
    double precision, intent(in) :: s(nPTReplicas)
    integer rep
    LOGICAL isfile
    character*32 fullName
    fullName=  'data/repHistory'
    inquire(file = fullName, exist=isfile)
    if (isfile) then
        OPEN (UNIT = 1, FILE = fullName, STATUS ='OLD', POSITION="append")
    else 
        OPEN (UNIT = 1, FILE = fullName, STATUS = 'new')
    endif

    write(1,*) "~~~~~~~~~~~exchange: ",nExchange,", IND:",IND,"~~~~~~~~~~~~~~~~~~~~"
    write(1,*) " rep | node|  up  | down |",&
               " chi  |  x_chi |",&
               " h_A  |  x_h_A |",&
               " mu   |  x_mu  |  s   |"
    do rep=1,nPTReplicas
        write(1,"(2I6,2f7.4,f7.4,f9.1,f7.4,f9.1,f7.4,f9.1,f7.4)"),&
                 rep, nodeNumber(rep), &
                 real(upSuccess(rep))/real(N_average), &
                 real(downSuccess(rep))/real(N_average),& 
                 cofMtrx(rep,1), xMtrx(rep,1),&
                 cofMtrx(rep,3), xMtrx(rep,3),&
                 cofMtrx(rep,2), xMtrx(rep,2),s(rep)
    enddo  
    Close(1)
    ! save which node each replica is running on
    fullName=  'data/nodeNumber'
    inquire(file = fullName, exist=isfile)
    if (isfile) then
        OPEN (UNIT = 1, FILE = fullName, STATUS ='OLD', POSITION="append")
    else 
        OPEN (UNIT = 1, FILE = fullName, STATUS = 'new')
    endif
    write(1,*), IND, nodeNumber
    Close(1)
    ! save chi values
    fullName=  'data/chi'
    inquire(file = fullName, exist=isfile)
    if (isfile) then
        OPEN (UNIT = 1, FILE = fullName, STATUS ='OLD', POSITION="append")
    else 
        OPEN (UNIT = 1, FILE = fullName, STATUS = 'new')
    endif
    write(1,*), IND, cofMtrx(:,1)
    Close(1)
    ! save field strength
    fullName=  'data/h_A'
    inquire(file = fullName, exist=isfile)
    if (isfile) then
        OPEN (UNIT = 1, FILE = fullName, STATUS ='OLD', POSITION="append")
    else 
        OPEN (UNIT = 1, FILE = fullName, STATUS = 'new')
    endif
    write(1,*), IND, cofMtrx(:,3)
    Close(1)
    ! save kap
    fullName=  'data/kap'
    inquire(file = fullName, exist=isfile)
    if (isfile) then
        OPEN (UNIT = 1, FILE = fullName, STATUS ='OLD', POSITION="append")
    else 
        OPEN (UNIT = 1, FILE = fullName, STATUS = 'new')
    endif
    write(1,*), IND, cofMtrx(:,5)
    Close(1)
    ! save mu
    fullName=  'data/mu'
    inquire(file = fullName, exist=isfile)
    if (isfile) then
        OPEN (UNIT = 1, FILE = fullName, STATUS ='OLD', POSITION="append")
    else 
        OPEN (UNIT = 1, FILE = fullName, STATUS = 'new')
    endif
    write(1,*), IND, cofMtrx(:,2)
    Close(1)
    ! save s
    fullName=  'data/s'
    inquire(file = fullName, exist=isfile)
    if (isfile) then
        OPEN (UNIT = 1, FILE = fullName, STATUS ='OLD', POSITION="append")
    else 
        OPEN (UNIT = 1, FILE = fullName, STATUS = 'new')
    endif
    write(1,*), IND, s
    Close(1)
end subroutine
