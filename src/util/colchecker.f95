
subroutine insertion_sort(n,a)
    implicit none
    integer n,i,j
    double precision a(n),x
    do 30 i=2,n
        x=a(i)
        j=i
10      j=j-1
        if (j.eq.0 .or. a(j).le.x) go to 20
        a(j+1)=a(j)
        go to 10
20      a(j+1)=x
30  continue
end

subroutine check_collisions(r, nt, has_collided, fpt_dist, time, col_type)
    implicit none
    integer, intent(in) :: nt, col_type
    double precision, intent(in) :: fpt_dist, time
    double precision, intent(in) :: r(nt,3)
    double precision, intent(inout) :: has_collided(nt, nt)
    if (col_type.eq.0) then
        return
    else if (col_type.eq.1) then
        call check_collisions_brute(r, nt, has_collided, fpt_dist, time)
    else if (col_type.eq.2) then
        call check_collisions_kd(r, nt, has_collided, fpt_dist, time)
    else if (col_type.eq.3) then
        call check_collisions_bb(r, nt, has_collided, fpt_dist, time)
    else if (col_type.eq.4) then
        print*, 'Unimplemented: col_type == 4: check_collisions_bin.'
        stop 1
        call check_collisions_bin(r, nt, has_collided, fpt_dist, time)
    end if
end

subroutine check_collisions_brute(r, nt, has_collided, fpt_dist, &
        time)
    integer nt
    double precision fpt_dist, time
    double precision r(nt,3), has_collided(nt, nt)
    !     check if the particles have collided
    do k1 = 1, nt
        do k2 = 1, nt
            if (has_collided(k1,k2).lt.0.0d0 .and. k1.ne.k2 &
                .and. abs(r(k1,1) - r(k2,1)) < fpt_dist &
                .and. abs(r(k1,2) - r(k2,2)) < fpt_dist &
                .and. abs(r(k1,3) - r(k2,3)) < fpt_dist) then
                has_collided(k1,k2) = time
            end if
        end do
    end do
end

subroutine check_collisions_kd(r, nt, has_collided, fpt_dist, time)
    use kdtree2_module, only : kdtree2, kdtree2_result, kdtree2_create, &
                                kdtree2_r_nearest_around_point

    integer nt, nfound, nalloc, k1, k2, i
    double precision fpt_dist, time
    double precision r(nt,3), has_collided(nt, nt)
    type(kdtree2), pointer :: col_tree
    type(kdtree2_result), allocatable :: kd_results(:)

    col_tree => kdtree2_create(r, rearrange = .true., sort = .false.)
    do k1 = 1,nt
        call kdtree2_r_nearest_around_point(col_tree, idxin = k1, &
            correltime = 1, r2 = fpt_dist, nfound = nfound, nalloc = nalloc, &
            results = kd_results)
        do i = 1,nfound
            k2 = kd_results(i)%idx
            if (has_collided(k1,k2) .lt. 0) then
                has_collided(k1,k2) = time
            endif
        enddo
    enddo
end

subroutine check_collisions_bb(r, nt, has_collided, fpt_dist, time)
! at each time point, we want to have 2 "pointer arrays", ind & indi
! r(ind(:,k),k) is in order for k in 1,2,3    i.e. [~,ind(:,1)] = sort(r(:,1)
! ind(indi(i,k),k) == i for k in 1,2,3
! at all time points
!
! acf means "after collision found" with ith bead in dimension d
    use params, only : dp

    implicit none

    integer, intent(in) :: nt
    double precision, intent(in) :: fpt_dist, time, r(nt,3)
    double precision, intent(inout) :: has_collided(nt, nt)

    integer :: neighbors(nt,nt)  ! most of array won't be used, probably
        ! neighbors(1:num_neighbors(i),i) holds neighbors of bead i for each i
        ! acf: neighbors(?,i) = j iff found in all three (d == 3, neighbor_triplet_keeper(j) == 2/3)
    integer :: num_neighbors(nt) ! to prevent O(nt^2) access to neighbors
        ! acf: num_neighbors(i)++ iff found in all three (d == 3, neighbor_triplet_keeper(j) == 2/3)
    integer :: neighbor_triplet_keeper(nt) ! O(nt)-space "hash table"
        ! every time we find that bead j is a neighbor in one of the three
        ! dimensions, then we increment neighbor_triplet_keeper(j), until we
        ! realize it's not a neighbor in one of the dimensions, or that it is
        ! in all three
        ! acf: neighbor_triplet_keeper(j)++
    integer :: neighbor_zeroer(nt), num_zeros ! to zero out "hash table" quickly
        ! reports that neighbor_triplet_keeper(neighbor_zeroer(1:num_zeros))
        ! should be zeroed after checking for all the neighbors of a particular
        ! bead
        ! acf: neighbor_zeroer(++num_zeros) = j if neighbor_triplet_keeper == 0/1
        ! better: "" "" if d == 1 (i.e. we're adding j to triplet array)
    integer, save, allocatable, dimension(:,:) :: ind, indi
    integer, save :: is_allocated = 0 ! "static" variable, allow initial setup
    integer :: curr_indi, curr_ind, i, j, d, rd0
    double precision :: rneighbor
    ! initialize ind and indi on first pass, requires O(n log n) sort
    if (is_allocated == 0) then
        is_allocated = 1
        allocate(ind(nt,3))
        allocate(indi(nt,3))
        do d = 1, 3
            do i = 1, nt
                ind(i,d) = i
                indi(i,d) = i
            enddo
            call qcolsort(nt, indi(:,d), ind(:,d), r(:,d))
        enddo
        ! ind and indi should satisfy desired property
    else
        do d = 1, 3
            call icolsort(nt, indi(:,d), ind(:,d), r(:,d))
        enddo
    endif
    ! initialize loop variables
    num_zeros = 0
    do i = 1, nt
        num_neighbors(i) = 0
        neighbor_triplet_keeper(i) = 0
        neighbor_zeroer(i) = 0
    enddo
    ! fills neighbors(:,i) with num_neighbors(i) indices of particles that the
    ! ith particle has collided with
    do i = 1, nt
        ! look at three dimensions one-by-one
        do d = 1, 3
            curr_indi = indi(i,d)
            curr_ind = ind(curr_indi,d)
            rd0 = r(curr_ind,d)
            ! first we're going to look for particles to the "right"
            j = 1
            if (curr_indi == nt) exit
            curr_ind = ind(curr_indi + j,d)
            rneighbor = r(curr_ind,d)
            do while (rneighbor < rd0 + fpt_dist)
! on the first pass, just mark that a collision happened in this coord
! then mark that index in the "hash table" as "needs zeroing"
                if (d == 1) then
                    num_zeros = num_zeros + 1
                    neighbor_zeroer(num_zeros) = curr_ind
                    ! basically do: neighbor_triplet_keeper(curr_ind) = neighbor_triplet_keeper(curr_ind) + 1
                    neighbor_triplet_keeper(curr_ind) = 1
! on the second pass, only mark if it also happened on the first pass
                elseif (d == 2) then
                    if (neighbor_triplet_keeper(curr_ind) == 1) then
                        ! basically do: neighbor_triplet_keeper(curr_ind) = neighbor_triplet_keeper(curr_ind) + 1
                        neighbor_triplet_keeper(curr_ind) = 2
                    endif
! on the third pass, make the ones that have hit the official colliders
                else ! (d == 3)
                    if (neighbor_triplet_keeper(curr_ind) == 2) then
                        num_neighbors(i) = num_neighbors(i) + 1
                        neighbors(num_neighbors(i),i) = curr_ind
                    endif
                endif
                j = j + 1
                if (curr_indi + j > nt) exit
                curr_ind = ind(curr_indi + j,d)
                rneighbor = r(curr_ind,d)
            enddo
            ! now we look for particles to the "left"
            j = 1
            if (curr_indi == 1) exit
            curr_ind = ind(curr_indi - j,d)
            rneighbor = r(curr_ind,d)
            do while (rneighbor > rd0 - fpt_dist)
! on the first pass, just mark that a collision happened in this coord
! then mark that index in the "hash table" as "needs zeroing"
                if (d == 1) then
                    num_zeros = num_zeros + 1
                    neighbor_zeroer(num_zeros) = curr_ind
                    ! basically do: neighbor_triplet_keeper(curr_ind) = neighbor_triplet_keeper(curr_ind) + 1
                    neighbor_triplet_keeper(curr_ind) = 1
! on the second pass, only mark if it also happened on the first pass
                elseif (d == 2) then
                    if (neighbor_triplet_keeper(curr_ind) == 1) then
                        ! basically do: neighbor_triplet_keeper(curr_ind) = neighbor_triplet_keeper(curr_ind) + 1
                        neighbor_triplet_keeper(curr_ind) = 2
                    endif
! on the third pass, make the ones that have hit the official colliders
                else ! (d == 3)
                    if (neighbor_triplet_keeper(curr_ind) == 2) then
                        num_neighbors(i) = num_neighbors(i) + 1
                        neighbors(num_neighbors(i),i) = curr_ind
                    endif
                endif
                j = j + 1
                if (curr_indi - j < 1) exit
                curr_ind = ind(curr_indi - j,d)
                rneighbor = r(curr_ind,d)
            enddo
        enddo
        ! zero out the "hash table"
        ! might be faster to sort neighbors to zero or zero entire thing for small nt?
        ! should probably check this at some point
        do j = 1, num_zeros
            neighbor_triplet_keeper(neighbor_zeroer(j)) = 0
        enddo
        num_zeros = 0
    enddo
    ! from the neighbors, num_neighbors arrays, we can rapidly extract
    ! exactly those elements that have collided
    ! neighbors(neighborj, beadi), num_neighbors(beadi)
    do i = 1, nt
        do j = 1, num_neighbors(i)
            if (has_collided(neighbors(j,i),i) < 0.0_dp) then
                has_collided(neighbors(j,i),i) = time
            endif
        enddo
    enddo
end

subroutine check_collisions_bin(r, nt, has_collided, fpt_dist, &
        time)
!TODO actually implement this function. copy/pasta of _brute for now
    integer nt
    double precision fpt_dist, time
    double precision r(nt,3), has_collided(nt, nt)
    !     check if the particles have collided
    do k1 = 1, nt
        do k2 = 1, nt
            if (has_collided(k1,k2).lt.0.0d0 .and. k1.ne.k2 &
                .and. abs(r(k1,1) - r(k2,1)) < fpt_dist &
                .and. abs(r(k1,2) - r(k2,2)) < fpt_dist &
                .and. abs(r(k1,3) - r(k2,3)) < fpt_dist) then
                has_collided(k1,k2) = time
            end if
        end do
    end do
end

