! ------------------------------------------------------
!
!   This module does binning for fast neighbor retrival
!
!        By Quinn MacPherson ~ summer 2017
!
!  Example Usage Pseudo-code:
!
!  use binning
!  type(binType) bin
!  double precision R(NT,3) ! all bead locations
!  
!  !  Set up binning object
!  setBinSize = [sideLength, sideLength, sideLength] ! size of bin
!  setMinXYZ = [0.0,0.0,0.0]  ! location of corner of bin
!  setBinShape = [10,10,10]   ! Specify first level of binning
!  call constructBin(bin,setBinShape,setMinXYZ,setBinSize) 
!  
!  ! Add 13th bead to bin
!  call addBead(bin,R,NT,13)
!
!  ! Remove 13th bead 
!  call removeBead(bin,R(13,:),13)
!
!  ! Get neighbors within radius of coordinate
!  double precision coordinate(3)  ! Position about which to search
!  double precision distances(maxN) ! Returned distances
!  integer neighbors(maxN) ! ID of neighboring beads
!  integer nn ! number of neighbors
!  call findNeighbors(bin,coordinate,radius,R,NT,maxN,neighbors,distances,nn)
!  
! ------------------------------------------------------
Module binning
    Implicit none
   
    ! The following are optimization constants.
    ! Their optimal values will, in general, depend on the application
    Integer, Parameter :: maxBeadsPerBin = 20
    Integer, Parameter :: subdivisionRatio = 4
    Integer, Parameter :: nBeadsToTrigerBinMerge = 3

    Type binType
        !-----------
        ! General variables
        !-----------
        logical isSuperBin
        integer numberOfBeads

        ! ----------
        ! simple bin only variables
        ! ----------
        integer beads(maxBeadsPerBin)

        ! ----------
        ! super-bin variables
        ! ----------
        integer numberOfBins   ! Number of sub-bins
        double precision minXYZ(3) ! xmin, ymin, zmin
        double precision binSize(3) ! size of system in floating point units
        integer binsShape(3)   ! bins in x,y and z dirrections
        type(binType), pointer :: bins(:) ! bins which are of type binType
        !type(binType), Allocatable, dimension(:):: bins ! (bin,element)


    end Type

contains

    subroutine constructBin(bin,setBinShape,setMinXYZ,setBinSize)
        Type(binType), intent(inout) :: bin
        integer, intent(in) :: setBinShape(3)
        double precision, intent(in) :: setBinSize(3) 
        double precision, intent(in) :: setMinXYZ(3) 
        bin%isSuperBin=.False.
        bin%numberOfBeads=0
        bin%binsShape=setBinShape
        bin%binSize=setBinSize
        bin%minXYZ=setMinXYZ
        bin%numberOfBins = setBinShape(1)*setBinShape(2)*setBinShape(3)
    end subroutine

    subroutine promoteBin(bin,R,NT)
        Type(binType), intent(inout) :: bin
        integer, intent(in) :: NT ! total number of beads
        double precision, intent(in) :: R(NT,3)  ! locations of all beads
        integer ii,ix,iy,iz,dd
        integer numberToTransfer
        double precision setBinSize(3) ! size of sub-bins
        double precision setMinXYZ(3) ! corners of sub-bins
        integer, parameter :: setBinShape(3) = [subdivisionRatio,&
                                                subdivisionRatio,&
                                                subdivisionRatio] ! How to sub-sub-divide bins
        ! calculate sub-bin size
        do dd=1,3
            setBinSize(dd)=bin%binSize(dd)/float(bin%binsShape(dd))
        enddo    
        ! Allocate bins
        Allocate(bin%Bins(bin%numberOfBins))
        ! Define bin locations
        ii=1
        setMinXYZ(1)=bin%minXYZ(1)
        do ix=1,bin%binsShape(1)
            setMinXYZ(2)=bin%minXYZ(2)
            do iy=1,bin%binsShape(2)
                setMinXYZ(3)=bin%minXYZ(3)
                do iz=1,bin%binsShape(3)
                    call constructBin(bin%bins(ii),setBinShape,setMinXYZ,setBinSize)
                    ii=ii+1
                    setMinXYZ(3)=setMinXYZ(3)+setBinSize(3)
                enddo
                setMinXYZ(2)=setMinXYZ(2)+setBinSize(2)
            enddo
            setMinXYZ(1)=setMinXYZ(1)+setBinSize(1)
        enddo
        ! Transfer beads from list to bins
        bin%isSuperBin=.True. ! now it is a supper bin
        numberToTransfer=bin%numberOfBeads ! number of beads to place in bins
        bin%numberOfBeads=0  ! add bead sill increment this
        do ii =1,numberToTransfer
            call addBead(bin,R,NT,bin%beads(ii))
        enddo
        !if (bin%numberOfBeads .ne. numberToTransfer) then
        !    print*, "Error in promote.  Values don't add up"
        !    print*, bin%numberOfBeads, numberToTransfer
        !    stop
        !endif
    end subroutine

    recursive subroutine demote(bin)
        Type(binType), intent(inout) :: bin
        integer ii
        integer numberOfMovedBeads,jj
        if (.not.bin%isSuperBin) then
            print*, "you tried to demote a simble bin"
            stop
        endif
        if (bin%numberOfBeads .ne. 0) then
            numberOfMovedBeads=0
            do ii=1,bin%numberOfBins
                if (bin%bins(ii)%isSuperBin) then
                    call demote(bin%bins(ii))
                endif
                do jj=1,bin%bins(ii)%numberOfBeads
                    numberOfMovedBeads=numberOfMovedBeads+1
                    bin%beads(numberOfMovedBeads)=bin%bins(ii)%beads(jj)
                enddo
            enddo
            !if (numberOfMovedBeads.ne.bin%numberOfBeads) then
            !    print*, "Error durring demote"
            !    print*, "Number of beads doesn't add up"
            !    stop
            !endif
            bin%isSuperBin=.False.
        else
            do ii=1,bin%numberOfBins
                if (bin%bins(ii)%isSuperBin) then
                    call demote(bin%bins(ii))
                endif
            enddo
            deallocate(bin%bins)
            bin%isSuperBin=.False.
        endif
    end subroutine

    recursive subroutine addBead(bin,R,NT,beadID)
        Type(binType), intent(inout) :: bin 
        integer, intent(in) :: NT ! total number of beads
        double precision, intent(in) :: R(NT,3)  ! locations of all beads
        integer, intent(in) :: beadID
        integer binIndex

        double precision, parameter :: eps=0.00000001
        integer dd
        do dd=1,3
            if (R(beadID,dd)+eps.lt.bin%minXYZ(dd)) then
                print*, "In Add"
                print*, "Error, location out of bin.  Under."
                print*, "location", R(beadID,:)
                print*, "minXYZ", bin%minXYZ
                print*, "maxXYZ", bin%minXYZ+bin%binSize
                stop
            endif
            if (R(beadID,dd)-eps-bin%binSize(dd)>bin%minXYZ(dd)) then
                print*, "In Add"
                print*, "Error, location out of bin. Over"
                print*, "location", R(beadID,:)
                print*, "minXYZ", bin%minXYZ
                print*, "maxXYZ", bin%minXYZ+bin%binSize
                stop
            endif
        enddo            
            
        if (bin%isSuperBin) then
            !print*, 'adding Bead to superbin minXYZ',bin%minXYZ
            !print*, 'getting index'
            !print*, 'beadID',beadID
            !print*, 'R:',R(beadID,:)
            call getBinIndex(bin,R(beadID,:),binIndex)
            !print*, 'index', binIndex
            call addBead(bin%bins(binIndex),R,NT,beadID)
            bin%numberOfBeads=bin%numberOfBeads+1
        else
            if (bin%numberOfBeads == maxBeadsPerBin) then
                call promoteBin(bin,R,NT)
                call addBead(bin,R,NT,beadID)
            else
                bin%numberOfBeads=bin%numberOfBeads+1
                bin%beads(bin%numberOfBeads)=beadID
            endif
        endif
    end subroutine

    recursive subroutine countBeads(bin,numberInBin)
        Type(binType), intent(in) :: bin 
        integer, intent(out) :: numberInBin
        integer ii, temp
        if (bin%isSuperBin) then
            numberInBin=0
            do ii=1,bin%numberOfBins
                call countBeads(bin%bins(ii),temp)
                numberInBin=numberInBin+temp
            enddo
            if (numberInBin.ne.bin%numberOfBeads) then
                print*, "Error in countBeads"
                print*, "Inconsistant number of beads"
                print*, "numberInBin", numberInBin,"numberOfBeads",bin%numberOfBeads
                stop
            endif
        else
            numberInBin=bin%numberOfBeads
        endif
        return
    end

    recursive subroutine removeBead(bin,location,beadID)
        Type(binType), intent(inout) :: bin 
        double precision, intent(in) :: location(3)
        integer, intent(in) :: beadID
        integer indexOfBeadToBeRemoved, ii
        integer binIndex ! index of sub-bin
        if (bin%numberOfBeads.lt.1) then
            print*, "Error: No beads to remove"
            stop
        endif
        indexOfBeadToBeRemoved=-1 ! -1 means not found (yet)
        if (bin%isSuperBin) then
            call getBinIndex(bin,location,binIndex)
            call removeBead(bin%bins(binIndex),location,beadID)
            bin%numberOfBeads=bin%numberOfBeads-1
            ! Possibly demote unneeded levels of binning
            if (bin%numberOfBeads .le. nBeadsToTrigerBinMerge) then
                call demote(bin)
            endif
        else
            ! find
            do ii=1,bin%numberOfBeads
                if (bin%beads(ii)==beadID) then
                    indexOfBeadToBeRemoved=ii
                    exit
                endif
            enddo
            ! Not found?
            if (indexOfBeadToBeRemoved==-1) then
                print*, "Error, bead not found"
                stop
            endif
            ! shift rest of list
            do ii=indexOfBeadToBeRemoved,bin%numberOfBeads-1
                bin%beads(ii)=bin%beads(ii+1)
            enddo
            bin%numberOfBeads=bin%numberOfBeads-1
        endif
    end subroutine

    subroutine getBinIndex(bin,location,binIndex)
        Type(binType), intent(in) :: bin
        double precision, intent(in) :: location(3)
        integer, intent(out) :: binIndex
        integer XYZ(3)
        integer dd
        
        !double precision, parameter :: eps=0.000001
        !do dd=1,3
        !    if (location(dd)+eps.lt.bin%minXYZ(dd)) then
        !        print*, "Error, location out of bin.  Under."
        !        print*, "location", location
        !        print*, "minXYZ", bin%minXYZ
        !        print*, "maxXYZ", bin%minXYZ+bin%binSize
        !        stop
        !    endif
        !    if (location(dd)-eps-bin%binSize(dd)>bin%minXYZ(dd)) then
        !        print*, "Error, location out of bin. Over"
        !        print*, "location", location
        !        print*, "minXYZ", bin%minXYZ
        !        print*, "maxXYZ", bin%minXYZ+bin%binSize
        !        stop
        !    endif
        !enddo            

        do dd=1,3
            XYZ(dd)=ceiling((location(dd)-bin%minXYZ(dd))*&
                bin%binsShape(dd)/bin%binSize(dd))
            if (XYZ(dd) == 0) then  ! corner case where x=0.000000
                XYZ(dd) = 1
            endif
            if (XYZ(dd) .gt. bin%binsShape(dd)) then ! in case of rounding error
                XYZ(dd) = bin%binsShape(dd)
            endif
        enddo
        binIndex=XYZ(3)+bin%binsShape(3)*((XYZ(2)-1)+bin%binsShape(2)*(XYZ(1)-1))
    end subroutine

    recursive subroutine findNeighbors(bin,location,radius,R,NT,&
                  maxNeighbors,neighbors,distances,nNeighbors)
        Type(binType), intent(in) :: bin
        double precision, intent(in) :: location(3) ! 
        double precision, intent(in) :: radius ! cuttoff on interaction distance
        integer, intent(in) :: NT ! total number of beads
        double precision, intent(in) :: R(NT,3)  ! locations of all beads
        integer, intent(in) :: maxNeighbors ! equal length of lists 
        integer, intent(inout) :: nNeighbors !number of neighbors found (so far)
        integer, intent(inout):: neighbors(maxNeighbors) ! list of bead ID's
        double precision, intent(inout) :: distances(maxNeighbors) ! list of |r-r| values
        integer dd,xx,yy,zz,ii,binIndex
        integer lower(3),upper(3)
        double precision distance


        if (bin%isSuperBin) then
            do dd = 1,3
                ! bound which sub-bins need to be searched
                lower(dd)=ceiling( (location(dd)-radius - bin%minXYZ(dd))*&
                                     float(bin%binsShape(dd))/bin%binSize(dd) )

                if (lower(dd)<1) then
                    lower(dd)=1
                endif
                upper(dd)=ceiling( (location(dd)+radius - bin%minXYZ(dd))*&
                                     float(bin%binsShape(dd))/bin%binSize(dd) )
                if (upper(dd)>bin%binsShape(dd)) then
                    upper(dd)=bin%binsShape(dd)
                endif
            enddo
            ! search relivant sub-bins
            do xx=lower(1),upper(1)
                do yy=lower(2),upper(2)
                    do zz=lower(3),upper(3)
                        binIndex=zz+bin%binsShape(3)*((yy-1)+bin%binsShape(2)*(xx-1))
                        call findNeighbors(bin%bins(binIndex),location,radius,R,NT,&
                            maxNeighbors,neighbors,distances,nNeighbors)
                    enddo
                enddo
            enddo
        else
            do ii=1,bin%numberOfBeads
                distance= (R(bin%beads(ii),1)-location(1))**2&
                         +(R(bin%beads(ii),2)-location(2))**2&
                         +(R(bin%beads(ii),3)-location(3))**2
                if (distance.le.(radius**2)) then
                    !if (nNeighbors .ge. maxNeighbors) then
                    !    print*, "Error: Too many neighbors"
                    !    stop
                    !endif
                    nNeighbors=nNeighbors+1
                    distances(nNeighbors)=distance
                    neighbors(nNeighbors)=bin%beads(ii)
                endif
            enddo
        endif
    end subroutine
end module
