#include "../defines.inc"
subroutine mc_internucleosome()
!Calculate the change in the polymer internucleosome attraction energy
! the actual energy function is in the nucleosome module file and is 
! based off of work from the de pablo group

! values from wlcsim_data
   use params, only: dp, NAN, wlc_U, wlc_nucleosomeWrap, wlc_VP, wlc_V&
                     , wlc_R, wlc_UP, wlc_RP, wlc_nPointsMoved, wlc_pointsMoved, &
                     wlc_bin, wlc_R_GJK
   use GJKAlgorithm, only: constructPolygonPrism
   use polydispersity, only: get_IP, first_bead_of_chain, last_bead_of_chain
   use nucleosome, only: internucleosome_energy
   use energies, only: energyOf, internucleosome_
   use binning, only: find_neighbors

   implicit none
   real(dp) delInt     ! change in internucleosome attraction energy
   real(dp), parameter :: cutoff = 25.0 ! nm
   integer i, j, k, IP
   integer ignore(WLC_P__NT) 
   real(dp) RGJK(3, WLC_P__NT) ! all centers for GJK
   real(dp) RALL(3, WLC_P__NT) ! all bead R
   real(dp) UALL(3, WLC_P__NT) ! all bead U
   real(dp) VALL(3, WLC_P__NT) ! all bead V
   real(dp) distances(WLC_P__NT) ! Returned distances
   real(dp) poly(WLC_P__GJK_POLYGON, 3)
   integer neighbors(WLC_P__NT) ! ID of neighboring beads
   integer nn ! number of neighbors

   ! set up
   delInt = 0
   RGJK = wlc_R_GJK
   RALL = wlc_R
   UALL = wlc_U
   VALL = wlc_V

   ! check energetics of old beads
   ignore = NAN
   do k = 1, wlc_nPointsMoved
      i = wlc_pointsMoved(k)
      if (wlc_nucleosomeWrap(i) == 1) cycle
      ignore(k) = i
      ! only check beads within cutoff distance, since n choose k grows quick
      if (WLC_P__NEIGHBOR_BINS) then
         nn = 0
         call find_neighbors(wlc_bin, RGJK(:, i), cutoff, RGJK, WLC_P__NT, &
                              WLC_P__NT, neighbors, distances, nn)
      else
         call findNeighbors(RGJK(:, i), cutoff, RGJK, WLC_P__NT, &
                              WLC_P__NT, neighbors, distances, nn)
      endif
      do j = 1, nn
         if (wlc_nucleosomeWrap(neighbors(j)) == 1 .or. ANY(ignore == neighbors(j))) cycle
         ! old config
         delInt = delInt - internucleosome_energy(RALL(:, i), RALL(:, neighbors(j)), &
                                                  UALL(:, i), UALL(:, neighbors(j)), &
                                                  VALL(:, i), VALL(:, neighbors(j)))
      enddo
   enddo
   ! replace old beads with new moved beads
   do j = 1, wlc_nPointsMoved
      i = wlc_pointsMoved(j)
      if (wlc_nucleosomeWrap(i) == 1) cycle ! only nucs
      IP = get_IP(i) 
      ! update real bead locations
      RALL(:, i) = wlc_RP(:, i)
      UALL(:, i) = wlc_UP(:, i)
      VALL(:, i) = wlc_VP(:, i)/norm2(wlc_VP(:, i))
      ! set center of nuc bead
      if ((i < last_bead_of_chain(IP))) then 
         if (isnan(wlc_RP(1, i + 1))) then 
            poly = constructPolygonPrism(wlc_RP(:, i), wlc_R(:, i + 1), &
                                          wlc_nucleosomeWrap(i), wlc_UP(:, i), &
                                          wlc_VP(:, i)/norm2(wlc_VP(:, i)), WLC_P__GJK_POLYGON)
         else 
            poly = constructPolygonPrism(wlc_RP(:, i), wlc_RP(:, i + 1), &
                                          wlc_nucleosomeWrap(i), wlc_UP(:, i), &
                                          wlc_VP(:, i)/norm2(wlc_VP(:, i)), WLC_P__GJK_POLYGON)
         endif
         RGJK(1, i) = sum(poly(:, 1)/WLC_P__GJK_POLYGON)
         RGJK(2, i) = sum(poly(:, 2)/WLC_P__GJK_POLYGON)
         RGJK(3, i) = sum(poly(:, 3)/WLC_P__GJK_POLYGON)
      endif
   enddo
   ! check energetics of new beads
   ignore = NAN
   do k = 1, wlc_nPointsMoved
      i = wlc_pointsMoved(k)
      if (wlc_nucleosomeWrap(i) == 1) cycle
      ignore(k) = i
      ! only check beads within cutoff distance, since n choose k grows quick
      if (WLC_P__NEIGHBOR_BINS) then
         nn = 0
         call find_neighbors(wlc_bin, RGJK(:, i), cutoff, RGJK, WLC_P__NT, &
                              WLC_P__NT, neighbors, distances, nn)
      else
         call findNeighbors(RGJK(:, i), cutoff, RGJK, WLC_P__NT, &
                              WLC_P__NT, neighbors, distances, nn)
      endif
      do j = 1, nn
         if (wlc_nucleosomeWrap(neighbors(j)) == 1 .or. ANY(ignore == neighbors(j))) cycle
         ! new config
         delInt = delInt + internucleosome_energy(RALL(:, i), RALL(:, neighbors(j)), &
                                                  UALL(:, i), UALL(:, neighbors(j)), &
                                                  VALL(:, i), VALL(:, neighbors(j)))
      enddo
   enddo
   energyOf(internucleosome_)%dx = delInt

END subroutine 