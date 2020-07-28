#include "../defines.inc"
!--------------------------------------------------------------
!
!     subroutine MC_sterics
!
!     Determine if moved beads intersect old beads.
!     Written by Nicole Pagane, Mar 2020
!     Inspired from Quinn's cylinder collision code
!     GJK algorithm used for collision detection
!--------------------------------------------------------------

subroutine MC_sterics(collisions, netSterics, MCTYPE)
! checks for steric collisons using the GJK algorithm, check under the 
! utils/sterics file for actual implentation of GJK. 
   use params, only: dp, NAN, wlc_RP, wlc_UP, wlc_VP, wlc_R, wlc_U, wlc_V, wlc_GJK, &
                     wlc_R_GJK, wlc_nucleosomeWrap, wlc_nPointsMoved, wlc_bin, wlc_basepairs, &
                     wlc_basepairs_prop, wlc_pointsMoved
   use GJKAlgorithm, only: constructPolygonPrism
   use polydispersity, only: get_IP, first_bead_of_chain, last_bead_of_chain
   use binning, only: addBead, removeBead, find_neighbors
   implicit none

   integer, intent(out) :: collisions ! number of steric collisions
   logical, intent(in) :: netSterics ! whether to check the sterics of the OG confiration
   integer, intent(in) :: MCTYPE ! current monte carlo move
   real(dp) RALL(3, WLC_P__NT) ! all bead R
   real(dp) UALL(3, WLC_P__NT) ! all bead U
   real(dp) VALL(3, WLC_P__NT) ! all bead V
   real(dp) SGJK(WLC_P__GJK_POLYGON, 3, WLC_P__NT) ! all vertices for GJK
   real(dp) RGJK(3, WLC_P__NT) ! all centers for GJK
   real(dp) basepairs(WLC_P__NT) ! basepairs
   real(dp) distances(WLC_P__NT) ! Returned distances
   integer neighbors(WLC_P__NT) ! ID of neighboring beads
   integer nn ! number of neighbors
   integer i, j, k, minr
   integer IP ! chain index 
   real(dp) poly(WLC_P__GJK_POLYGON, 3)
   integer ignore(WLC_P__NT), ignore_bin(WLC_P__NT)

! only if the MC move moved a bead
   if (wlc_nPointsMoved > 0) then
      ! set up for collision searching
      RGJK = wlc_R_GJK
      SGJK = wlc_GJK
      RALL = wlc_R
      UALL = wlc_U
      VALL = wlc_V
      collisions = 0
      if (netSterics) then 
         k = 1
         ignore = NAN
         ! check for neighbors on old beads
         do j = 1, wlc_nPointsMoved
            i = wlc_pointsMoved(j)
            IP = get_IP(i) 
            ! check segment preceding bead
            if ((ANY(ignore == i - 1) .eqv. .false.) .AND. (i > first_bead_of_chain(IP))) then 
               ignore(k) = i - 1
               k = k + 1
               call findNeighbors(RGJK(:, i - 1), 2*WLC_P__GJK_RADIUS, &
                                  RGJK, WLC_P__NT, WLC_P__NT, neighbors, distances, nn)
               ! check for collisions
               call sterics_check(collisions, RALL, UALL, VALL, SGJK, wlc_basepairs, ignore, i - 1, &
                                  nn, neighbors(1:nn), distances(1:nn), .true.)
            endif
            ! check succeeding bead
            if ((ANY(ignore == i) .eqv. .false.) .AND. (i < last_bead_of_chain(IP))) then 
               ignore(k) = i
               k = k + 1
               call findNeighbors(RGJK(:, i), 2*WLC_P__GJK_RADIUS, &
                                  RGJK, WLC_P__NT, WLC_P__NT, neighbors, distances, nn)
               ! check for collisions
               call sterics_check(collisions, RALL, UALL, VALL, SGJK, wlc_basepairs, ignore, i, &
                                  nn, neighbors(1:nn), distances(1:nn), .true.)
            endif
         enddo
      endif
      ! set basepairs vector
      if (WLC_P__MOVEON_NUCLEOSOMESLIDE == 1 .AND. MCTYPE == 13) then 
         basepairs = wlc_basepairs_prop
      else
         basepairs = wlc_basepairs
      endif
      ! replace old beads with new moved beads
      k = 1
      ignore = NAN
      do j = 1, wlc_nPointsMoved
         i = wlc_pointsMoved(j)
         IP = get_IP(i) 
         ! update real bead locations
         RALL(:, i) = wlc_RP(:, i)
         UALL(:, i) = wlc_UP(:, i)
         VALL(:, i) = wlc_VP(:, i)/norm2(wlc_VP(:, i))
         ! check segment preceding bead
         if ((i > first_bead_of_chain(IP)) .AND. (ANY(ignore == i - 1) .eqv. .false.)) then
            ignore(k) = i - 1
            k = k + 1
            if (isnan(wlc_RP(1, i - 1))) then 
               poly = constructPolygonPrism(wlc_R(:, i - 1), wlc_RP(:, i), &
                                            wlc_nucleosomeWrap(i - 1), wlc_U(:, i - 1), &
                                            wlc_V(:, i - 1), WLC_P__GJK_POLYGON)
            else
               poly = constructPolygonPrism(wlc_RP(:, i - 1), wlc_RP(:, i), &
                                            wlc_nucleosomeWrap(i - 1), wlc_UP(:, i - 1), &
                                            wlc_VP(:, i - 1)/norm2(wlc_VP(:, i - 1)), WLC_P__GJK_POLYGON)
            endif
            SGJK(:, :, i - 1) = poly
            RGJK(1, i - 1) = sum(poly(:, 1)/WLC_P__GJK_POLYGON)
            RGJK(2, i - 1) = sum(poly(:, 2)/WLC_P__GJK_POLYGON)
            RGJK(3, i - 1) = sum(poly(:, 3)/WLC_P__GJK_POLYGON)
         endif
         ! check succeeding bead
         if ((i < last_bead_of_chain(IP)) .AND. (ANY(ignore == i) .eqv. .false.)) then 
            ignore(k) = i
            k = k + 1
            if (isnan(wlc_RP(1, i + 1))) then 
               poly = constructPolygonPrism(wlc_RP(:, i), wlc_R(:, i + 1), &
                                            wlc_nucleosomeWrap(i), wlc_UP(:, i), &
                                            wlc_VP(:, i)/norm2(wlc_VP(:, i)), WLC_P__GJK_POLYGON)
            else 
               poly = constructPolygonPrism(wlc_RP(:, i), wlc_RP(:, i + 1), &
                                            wlc_nucleosomeWrap(i), wlc_UP(:, i), &
                                            wlc_VP(:, i)/norm2(wlc_VP(:, i)), WLC_P__GJK_POLYGON)
            endif
            SGJK(:, :, i) = poly
            RGJK(1, i) = sum(poly(:, 1)/WLC_P__GJK_POLYGON)
            RGJK(2, i) = sum(poly(:, 2)/WLC_P__GJK_POLYGON)
            RGJK(3, i) = sum(poly(:, 3)/WLC_P__GJK_POLYGON)
         endif
      enddo
      collisions = -collisions
      ! set up ignore beads
      if (WLC_P__NEIGHBOR_BINS .AND. (MCTYPE /= 13) .AND. (netSterics .eqv. .false.)) then
         ignore_bin = NAN
         ignore_bin(1:wlc_nPointsMoved) = wlc_pointsMoved(1:wlc_nPointsMoved)
      endif
      k = 1
      ignore = NAN
      ! check for neighbors on new beads
      do j = 1, wlc_nPointsMoved
         i = wlc_pointsMoved(j)
         IP = get_IP(i) 
         ! check segment preceding bead
         if ((ANY(ignore == i - 1) .eqv. .false.) .AND. (i > first_bead_of_chain(IP))) then 
            ignore(k) = i - 1
            k = k + 1
            if (WLC_P__NEIGHBOR_BINS .AND. (MCTYPE /= 13) .AND. (netSterics .eqv. .false.)) then
               nn = 0
               call find_neighbors(wlc_bin, RGJK(:, i - 1), 2*WLC_P__GJK_RADIUS, &
                                   wlc_R_GJK, WLC_P__NT, WLC_P__NT, neighbors, distances, nn)
               ! check for collisions
               call sterics_check(collisions, RALL, UALL, VALL, SGJK, basepairs, ignore_bin, i - 1, &
                                  nn, neighbors(1:nn), distances(1:nn), netSterics)
            else
               call findNeighbors(RGJK(:, i - 1), 2*WLC_P__GJK_RADIUS, &
                                  RGJK, WLC_P__NT, WLC_P__NT, neighbors, distances, nn)
               ! check for collisions
               call sterics_check(collisions, RALL, UALL, VALL, SGJK, basepairs, ignore, i - 1, &
                                  nn, neighbors(1:nn), distances(1:nn), netSterics)
            endif
         endif
         ! check succeeding bead
         if ((ANY(ignore == i) .eqv. .false.) .AND. (i < last_bead_of_chain(IP))) then 
            ignore(k) = i
            k = k + 1
            if (WLC_P__NEIGHBOR_BINS .AND. (MCTYPE /= 13) .AND. (netSterics .eqv. .false.)) then
               nn = 0
               call find_neighbors(wlc_bin, RGJK(:, i), 2*WLC_P__GJK_RADIUS, &
                                   wlc_R_GJK, WLC_P__NT, WLC_P__NT, neighbors, distances, nn)
               ! check for collisions
               call sterics_check(collisions, RALL, UALL, VALL, SGJK, basepairs, ignore_bin, i, &
                                  nn, neighbors(1:nn), distances(1:nn), netSterics)
            else
               call findNeighbors(RGJK(:, i), 2*WLC_P__GJK_RADIUS, &
                                  RGJK, WLC_P__NT, WLC_P__NT, neighbors, distances, nn)
               ! check for collisions
               call sterics_check(collisions, RALL, UALL, VALL, SGJK, basepairs, ignore, i, &
                                  nn, neighbors(1:nn), distances(1:nn), netSterics)
            endif
         endif
      enddo
   endif
END subroutine MC_sterics

subroutine sterics_check(collisions, RALL, UALL, VALL, SGJK, basepairs, ignore, ii, nn, neighbors, distances, netSterics)
! sterics check subroutine to check for different types of collisions (i.e. nuc vs nuc, nuc vs dna, dna vs dna)
! iterates through the nearest neighbors of the specified bead to check for any collisions
   use GJKAlgorithm, only: GJK, constructPolygonPrism
   use params, only: dp, wlc_basepairs, wlc_nucleosomeWrap, wlc_pointsMoved, wlc_nPointsMoved
   use nucleosome, only: nucleosome_prop
   use polydispersity, only: get_IP, first_bead_of_chain, last_bead_of_chain
   implicit none

   integer, intent(inout) :: collisions
   real(dp), intent(in) :: RALL(3, WLC_P__NT) ! all bead R
   real(dp), intent(in) :: UALL(3, WLC_P__NT) ! all bead U
   real(dp), intent(in) :: VALL(3, WLC_P__NT) ! all bead V
   real(dp), intent(in) :: SGJK(WLC_P__GJK_POLYGON, 3, WLC_P__NT) ! all vertices for GJK
   real(dp), intent(in) :: basepairs(WLC_P__NT) ! basepair discretization
   integer, intent(in) :: ignore(WLC_P__NT) ! beads that have already been checked, i.e. ignore
   integer, intent(in) :: ii ! index of moved bead
   integer, intent(in) :: nn ! number of neighbors
   integer, intent(in) :: neighbors(nn) ! ID of neighboring beads
   real(dp), intent(in) :: distances(nn) ! ID of neighboring beads
   logical, intent(in) :: netSterics ! whether or not to check all collisions, i.e. previous and new
   logical :: iiIsNucleosome, jjIsNucleosome ! whether or not the moved bead is a nucleosome
   logical :: jjGreaterThanii ! as it says
   integer, parameter :: s = WLC_P__GJK_POLYGON ! num sides of desired polygon
   real(dp), dimension(3) :: tempR, tempU, tempV
   real(dp), dimension(s, 3) :: poly1Plus, poly2Plus, poly1ExitDNA, poly2ExitDNA
   integer jj

   if (ii == last_bead_of_chain(get_IP(ii))) return

   ! determine identity of moving bead
   if (wlc_nucleosomeWrap(ii) /= 1) then ! is nucleosome
      iiIsNucleosome = .TRUE.
      call nucleosome_prop(UALL(:, ii), VALL(:, ii), RALL(:, ii), basepairs(ii), wlc_nucleosomeWrap(ii), &
                          tempU, tempV, tempR)
      poly1ExitDNA = constructPolygonPrism(tempR, RALL(:, ii + 1), 1, tempU, tempV, s)
   else
      iiIsNucleosome = .FALSE.
   endif

   ! construct polygon for i to i+1 bead (should be on virtual bead, so i+1 is safe near end of chain)
   poly1Plus = SGJK(:, :, ii)

   ! iterate through all possible collisions 
   do jj = 1, nn
      if (ANY(ignore == neighbors(jj)) .or. (neighbors(jj) == last_bead_of_chain(get_IP(neighbors(jj))))) cycle
      ! determine identity of potentially collided bead
      if (wlc_nucleosomeWrap(neighbors(jj)) /= 1) then ! is nucleosome
         jjIsNucleosome = .TRUE.
         call nucleosome_prop(UALL(:, neighbors(jj)), VALL(:, neighbors(jj)), RALL(:, neighbors(jj)), &
                             basepairs(neighbors(jj)), wlc_nucleosomeWrap(neighbors(jj)), tempU, tempV, tempR)
         poly2ExitDNA = constructPolygonPrism(tempR, RALL(:, neighbors(jj) + 1), 1, tempU, tempV, s)
      else
         jjIsNucleosome = .FALSE.
      endif
      ! construct forward polygon
      poly2Plus = SGJK(:, :, neighbors(jj))
      ! figure out ordering, THIS IS IMPORTANT FOR CONSISTENT GJK CONVERSION
      if (neighbors(jj) > ii) then 
         jjGreaterThanii = .true.
      else
         jjGreaterThanii = .false.
      endif
      ! check identity of all other beads in chain 
      if (iiIsNucleosome .AND. jjIsNucleosome) then ! nuc i + nuc j 
         ! check for collision of nucleosomes
         if (jjGreaterThanii) then ! ORDER MATTERS
            collisions = collisions + GJK(poly1Plus, poly2Plus, s)
         else
            collisions = collisions + GJK(poly2Plus, poly1Plus, s)
         endif
         ! check for exit DNA collisions within relevant distance cutoff
         if (distances(jj) < (2*basepairs(neighbors(jj))*WLC_P__LENGTH_PER_BP) + WLC_P__GJK_RADIUS .OR. &
             distances(jj) < (2*basepairs(ii)*WLC_P__LENGTH_PER_BP) + WLC_P__GJK_RADIUS ) then 
            ! check for collision of nucleosome i with exit DNA j ,
            ! collision of exit DNA i with nuclesome j, and
            ! colllision of exit DNA i and j
            if ((ii + 1 < neighbors(jj)) .OR. (neighbors(jj) + 1 < ii) ) then
               if (jjGreaterThanii) then ! ORDER MATTERS
                  collisions = collisions + GJK(poly1Plus, poly2ExitDNA, s)
                  collisions = collisions + GJK(poly1ExitDNA, poly2Plus, s)
                  collisions = collisions + GJK(poly1ExitDNA, poly2ExitDNA, s)
               else
                  collisions = collisions + GJK(poly2ExitDNA, poly1Plus, s)
                  collisions = collisions + GJK(poly2Plus, poly1ExitDNA, s)
                  collisions = collisions + GJK(poly2ExitDNA, poly1ExitDNA, s)
               endif
            endif
         endif
      ! check for moved bead nuc + DNA within relevant distance cutoff
      else if (iiIsNucleosome .AND. (jjIsNucleosome .EQV. .FALSE.) .AND. & 
               distances(jj) < (2*basepairs(neighbors(jj))*WLC_P__LENGTH_PER_BP) + WLC_P__GJK_RADIUS) then ! nuc i + DNA j 
         ! ignore 10bp nearest nuc, this is inspired from what elena koslover did in fibermodel
         if ((neighbors(jj) + 1 <= ii - 1 .AND. sum(basepairs(neighbors(jj) + 1:ii - 1)) > 10) .OR. &
             (neighbors(jj) - 1 >= ii + 1 .AND. sum(basepairs(ii:neighbors(jj) - 1)) > 10)) then 
            ! check for collision of nucleosome i with DNA j,
            ! and for collision of exit DNA i with DNA j
            if (jjGreaterThanii) then ! ORDER MATTERS
               collisions = collisions + GJK(poly1Plus, poly2Plus, s)
               collisions = collisions + GJK(poly1ExitDNA, poly2Plus, s)
            else
               collisions = collisions + GJK(poly2Plus, poly1Plus, s)
               collisions = collisions + GJK(poly2Plus, poly1ExitDNA, s)
            endif
         endif 
      ! check for moved bead DNA + nuc within relevant cutoff distance
      else if ((iiIsNucleosome .EQV. .FALSE.) .AND. jjIsNucleosome .AND. &
                distances(jj) < (2*basepairs(ii)*WLC_P__LENGTH_PER_BP) + WLC_P__GJK_RADIUS) then ! DNA i  + nuc j
         ! ignore 10bp nearest nuc, this is inspired from what elena koslover did in fibermodel
         if ((ii + 1 <= neighbors(jj) - 1 .AND. sum(basepairs(ii + 1:neighbors(jj) - 1)) > 10) .OR. &
             (ii - 1 >= neighbors(jj) + 1 .AND. sum(basepairs(neighbors(jj):ii - 1)) > 10)) then 
            ! check for collision of DNA i with nucleosome j,
            ! and for collision DNA i with exit DNA j
            if (jjGreaterThanii) then ! ORDER MATTERS
               collisions = collisions + GJK(poly1Plus, poly2Plus, s)
               collisions = collisions + GJK(poly1Plus, poly2ExitDNA, s)
            else
               collisions = collisions + GJK(poly2Plus, poly1Plus, s)
               collisions = collisions + GJK(poly2ExitDNA, poly1Plus, s)
            endif
         endif
      ! check for moved dna and dna within relevant cutoff distance
      else if (distances(jj) < 1.5*max(basepairs(neighbors(jj)), basepairs(ii))*WLC_P__LENGTH_PER_BP) then ! DNA i  + DNA j
         if ((ii + 1 < neighbors(jj)) .OR. (neighbors(jj) + 1 < ii)) then
            ! check for collision
            if (jjGreaterThanii) then ! ORDER MATTERS
               collisions = collisions + GJK(poly1Plus, poly2Plus, s)
            else
               collisions = collisions + GJK(poly2Plus, poly1Plus, s)
            endif
         endif
      endif
      ! return early if collision found and there no preexisting clashes
      if (collisions > 0 .and. (netSterics .eqv. .false.)) return
   enddo
END subroutine sterics_check 

subroutine findNeighbors(pos, radius, beads, nBeads, neighboringMax, neighbors, distances, nn)
! default to this instead of quinn's binning find_neighbors if binning is turned off
! rather than using recursive binning, this functon checks pairwise distances between beads.
! this will be SLOW for very many beads, but if you are under 1000-ish then it shouldn't be 
! terrible; otherwise, turn WLC_P__NEIGHBOR_BINS on
   use params, only: dp
   use vector_utils, only: distance
   implicit none

   real(dp), intent(in) :: pos(3) ! position of interest
   real(dp), intent(in) :: radius ! cutoff distance for relevant collisions
   integer, intent(in) :: nBeads ! total number of beads in simulation
   real(dp), intent(in) :: beads(3, nBeads) ! positions of all other beads
   integer, intent(in) :: neighboringMax ! max number of neighbors that could be found
   integer, intent(out) :: neighbors(neighboringMax) ! list of neighbors within cutoff distance
   real(dp), intent(out) :: distances(neighboringMax) ! list of distances of found neighbors
   integer, intent(out) :: nn ! total number of neighbors found within cutoff distance
   integer i 
   real(dp) dist

   nn = 0
   do i = 1, nBeads
      dist = distance(pos, beads(:, i))
      if (dist <= radius) then 
         nn = nn + 1
         neighbors(nn) = i
         distances(nn) = dist
      endif
   enddo

end subroutine findNeighbors
! ---------------------------------------------------------------*