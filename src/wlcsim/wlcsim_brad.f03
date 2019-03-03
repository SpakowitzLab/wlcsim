#include "../defines.inc"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!These subroutines allow for a Monte Carlo simulation of 
!a polymer chain to be performed
!Cases for parallel tempering (replica exchange) and without are included
!Without parallel tempering, nSTEPS are performed.
!With parallel tempering, nSTEP are performed between nREPLICAexchange iterations 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!wlcsim subroutine for performing monte carlo simulation
subroutine wlcsim_brad(wlc_p)
! values from wlcsim_data
use params, only: wlc_eelasREPLICAS, wlc_id, wlc_R, wlc_Wrs, wlc_nodeNUMBER&
    , wlc_LKs, wlc_Wr, wlc_nLKs
  use params
  use mersenne_twister
  use energies, only:energyOf, bend_, stretch_, shear_, twist_ 
#if MPI_VERSION
  use mpi
#endif

  implicit none

  !Structures for simulation
  type(wlcsim_params) :: wlc_p
  
  !MPI status variables
#if MPI_VERSION
  integer (kind = 4)  status(MPI_STATUS_SIZE)
#endif

  integer (kind = 4) error
  !Integers for commuinication among replicas
  integer rep
  integer dest
  integer source
  integer LK
  real(dp) Wr
  real(dp) eelas(4)
  integer repSTART,repEND
  !Counting variables
  integer repinD

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Run simulation
  !Cases:
  !parallel tempering on
  !parallel tempering off
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  !If parallel tempering is on
  if (WLC_P__PT_TWIST) then

#if MPI_VERSION
     !Cases:
     !Head node (id = 0)
     !Worker node (id /= 0)

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !Head node
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     if (wlc_id.eq.0) then
       
        !During replica exchange, replicas alternate between switching with
        !replica above and replica below. Initially, set the replica looking above
        !to the be first one
        repSTART = 1
        repEND = wlc_nLKs - 1

        !Begin replica exchange loop
        do repinD = 1,WLC_P__NREPLICAEXCHANGEPERSAVEPOINT
           
           !Send out Lk to the worker nodes to begin their simulations
           do rep = 1,wlc_nLKs
              LK = wlc_LKs(rep)
              dest = wlc_nodeNUMBER(rep)
              call MPI_Send (LK,1, MPI_integer, dest,   0, &
                   MPI_COMM_WORLD,error )
           enddo

           !LK, Wr, and EELAS from the worker nodes to perform
           !replica exchange

           do rep = 1,wlc_nLKs
              LK = wlc_LKs(rep)
              source = wlc_nodeNUMBER(rep)
              call MPI_Recv(Wr,1, MPI_doUBLE_PRECISION, source, 0, &
                   MPI_COMM_WORLD,status,error )
              call MPI_Recv(eelas,4, MPI_doUBLE_PRECISION, source, 0, &
                   MPI_COMM_WORLD,status,error )
              wlc_Wrs(rep) = Wr
              wlc_eelasREPLICAS(rep,:) = eelas
           enddo

           !peform replica exchange
           call replicaEXCHANGE_brad(wlc_p)
           

        enddo

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !Worker nodes
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     else
        !Begin replica exchange loop
        do repinD = 1,WLC_P__NREPLICAEXCHANGEPERSAVEPOINT
           !Receive LK from head node to begin simulation
           source = 0
           call MPI_Recv(LK,1, MPI_integer, source,   0, &
                MPI_COMM_WORLD,status,error )
           wlc_p%LK = LK

           !Run a monte carlo simulation for NSTEPS
           call MCsim(wlc_p,WLC_P__STEPSPEREXCHANGE)
         
           !Recalculate structural quantities and energies
           call writhe(wlc_R,WLC_P__NB, wlc_Wr)
           call energy_elas(eelas,wlc_p)
           energyOf(bend_)%E    =eelas(1)
           energyOf(shear_)%E   =eelas(2)
           energyOf(stretch_)%E =eelas(3)
           energyOf(twist_)%E   =eelas(4)

           !Communicate with the head node for replica exchange
         
           !Send back writhe and elastic energy back to the head node
           call MPI_SEND(wlc_Wr,1, MPI_doUBLE_PRECISION, source, 0, &
                MPI_COMM_WORLD,error )
           call MPI_SEND(eelas,4, MPI_doUBLE_PRECISION, source, 0, &
                MPI_COMM_WORLD,error )


           !Receive new LKs from the head node after replica exchange
           call MPI_Recv(LK,1, MPI_integer, source,   0, &
                MPI_COMM_WORLD,status,error )
           wlc_p%LK = LK
           
        enddo
     endif

#endif
  !No parallel tempering
  else
     !Run a MC simulation for nstepsPerExchange

     call MCsim(wlc_p,WLC_P__STEPSPEREXCHANGE)

     !Recalculate structural quantities and energies
     call writhe(wlc_R,WLC_P__NB, wlc_Wr)
     call energy_elas(eelas,wlc_p)
     energyOf(bend_)%E    =eelas(1)
     energyOf(shear_)%E   =eelas(2)
     energyOf(stretch_)%E =eelas(3)
     energyOf(twist_)%E   =eelas(4)

  end if
end subroutine wlcsim_brad


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Replica exchange subroutine
!after performing a monte carlo simulation
!this subroutine loops over pairs of replicas
!with adjacent linking numbers and tests for 
!exchange of their linking numbers 

#if MPI_VERSION
subroutine replicaEXCHANGE_brad(wlc_p)
! values from wlcsim_data
use params, only: wlc_replicaSTART, wlc_eelasREPLICAS, wlc_Wrs, wlc_rand_stat, wlc_nTRIALdown&
    , wlc_replicaEND, wlc_LKs, wlc_nTRIALup, wlc_nodeNUMBER, wlc_nSWAPdown, wlc_nLKs&
    , wlc_nSWAPup
  use params
  use mersenne_twister
  use mpi
  
  implicit none

  !Global structures for wlc params and data
  type(wlcsim_params) :: wlc_p

  !Variables for replica exchange
  integer rep             !replica index
  real(dp) deEXCHANGE     !change in energy for exchange for linking number
  real(dp) urand(1)       !random number for exchange
  real(dp) prob               !probability of exchange
  integer tempLK
  integer dest
  integer (kind = 4) error
  integer LK
  !Loop over every other replica and check for it to exchange with the one above it
  !For the replica exchange, each replica keeps its polymer conformation, and exchanges
  !linking number

  do rep = wlc_replicaSTART, wlc_replicaEND,2

     !Update number of trial exchanges
     wlc_nTRIALup(rep) = wlc_nTRIALup(rep) + 1 
     wlc_nTRIALdown(rep + 1) = wlc_nTRIALdown(rep + 1) + 1 

     deEXCHANGE = 0.0_dp

     !Include energy from original twist energy of each replica
     deEXCHANGE = -wlc_eelasREPLICAS(rep,4) - wlc_eelasREPLICAS(rep + 1,4)

     !Add change in twist energy due to rep taking LK from above
     deEXCHANGE = deEXCHANGE +((2.0_dp*pi*(wlc_LKs(rep + 1)-wlc_Wrs(rep)))**2.0_dp)*WLC_P__LT/(2.0_dp*WLC_P__L)

     !Add change in twist energy due to rep + 1 takign LK from below
     deEXCHANGE = deEXCHANGE +((2.0_dp*pi*(wlc_LKs(rep)-wlc_Wrs(rep + 1)))**2.0_dp)*WLC_P__LT/(2.0_dp*WLC_P__L)

     !Generate a random number and test for exchange
     call random_number(urand,wlc_rand_stat)
     prob = exp(-deEXCHANGE)

     !Exchange linking numbers according to Metropolis-Hastings criterion
     !This means that the node number associated with each LK is swapped
     if (urand(1).le.prob) then

        tempLK = wlc_nodeNUMBER(rep)
        wlc_nodeNUMBER (rep) = wlc_nodeNUMBER (rep + 1)
        wlc_nodeNUMBER(rep + 1) = tempLK

        !update number of swaps
        wlc_nSWAPup(rep) = wlc_nSWAPup(rep) + 1
        wlc_nSWAPdown(rep + 1) = wlc_nSWAPdown(rep + 1) + 1

     endif
  enddo

  !Alternatve value of replica start and replica end
  if (wlc_replicaSTART.eq.1) then
     wlc_replicaSTART = 2
  else
     wlc_replicaSTART = 1
  endif

  if (wlc_replicaEND.eq.wlc_nLKs-1) then
     wlc_replicaEND = wlc_nLKs-2
  else
     wlc_replicaEND = wlc_nLKs -1
  endif

  !Tell nodes their new LKs

  do rep = 1,wlc_nLKs
     LK = wlc_LKs(rep)
     dest = wlc_nodeNUMBER(rep)
     call MPI_Send (LK,1, MPI_integer, dest,   0, &
          MPI_COMM_WORLD,error )
  enddo


end subroutine replicaEXCHANGE_brad

#endif
        


  
  


