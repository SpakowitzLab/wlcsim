program test_binning
use binning
implicit none
integer NT
double precision, allocatable, dimension(:,:):: R
double precision sideLength
integer ii,jj, seed
type(binType) bin
double precision radius

double precision setBinSize(3)
double precision setMinXYZ(3)
integer setBinShape(3)

integer, parameter :: maxNeighbors = 1000
double precision distances(maxNeighbors)
integer neighbors(maxNeighbors)
integer nNeighbors
integer temp

double precision distance
integer totalNumberOfNeighbors

real :: start, finish

!---------------------
!
!  Set of fake simulation
!
!---------------------
NT=100000
sideLength=12.234
radius=0.4233
seed=86456
call srand(seed)
allocate(R(NT,3))
do ii=1,NT
    R(ii,1)=rand()*sideLength
    R(ii,2)=rand()*sideLength
    R(ii,3)=rand()*sideLength
enddo

!-------------------------
!
!  Set up binning object
!
!--------------------------
setBinSize = [sideLength, sideLength, sideLength]
setMinXYZ = [0.0,0.0,0.0]
setBinShape = [6,6,6]
call constructBin(bin,setBinShape,setMinXYZ,setBinSize)


!------------------------
!
!   Add beads to bin object
!
!-------------------------
print*, 'adding',NT,' beads to bin...'
call cpu_time(start)
do ii=1,NT
    call addBead(bin,R,NT,ii)
enddo
call cpu_time(finish)
print*, "Average time per addition", (finish-start)*(10**6)/NT, " microseconds"

!-------------------
!
! Move some beads Beads
!
!------------------
do ii=1,NT/2 
    call removeBead(bin,R(ii,:),ii)
    R(ii,1)=rand()*sideLength
    R(ii,2)=rand()*sideLength
    R(ii,3)=rand()*sideLength
enddo
do ii=1,NT/2
    call addBead(bin,R,NT,ii)
enddo


!--------------------
! Run double check
!------------------
print*, "double checking..."
call countBeads(bin,temp)
print*, "found total of ",temp,"beads"


!-----------------------
!
!  Find neighbors of all beads
!
!-------------------------
print*, 'now Find neighbors...'
totalNumberOfNeighbors=0
call cpu_time(start)
do ii=1,NT
    nNeighbors=0
    call findNeighbors(bin,R(ii,:),radius,R,NT,&
                          maxNeighbors,neighbors,distances,nNeighbors)
    totalNumberOfNeighbors=totalNumberOfNeighbors+nNeighbors
enddo
print*, "found", totalNumberOfNeighbors, "neighbors. "
print*, "Including self there are ~",totalNumberOfNeighbors/float(NT)," neighbors ber bead" 
call cpu_time(finish)
print*, "Average time per search point", (finish-start)*(10**6)/NT, " microseconds"


!-------------------------------
!
!  Check against brute force algorithm
!
!------------------------------------
print*, "now brute force"
totalNumberOfNeighbors=NT ! self interactions
do ii=1,NT
    do jj=ii+1,NT
        distance=((R(ii,1)-R(jj,1))**2+&
                  (R(ii,2)-R(jj,2))**2+&
                  (R(ii,3)-R(jj,3))**2)
        if (distance<(radius**2)) then
            totalNumberOfNeighbors=totalNumberOfNeighbors+2
        endif
    enddo
enddo
print*, "found", totalNumberOfNeighbors, "neighbors"

!-----------------------------------
!
!   Remove all the beads
!
!---------------------------------
print*, 'removing beads...'
do ii=1,NT
    call removeBead(bin,R(ii,:),ii)
enddo


end program
