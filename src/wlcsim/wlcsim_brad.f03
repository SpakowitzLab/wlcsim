!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!These subroutines allow for a Monte Carlo simulation of
!a polymer chain to be performed
!Cases for parallel tempering (replica exchange) and without are included
!Without parallel tempering, nSTEPS are performed.
!With parallel tempering, nSTEP are performed between nREPLICAexchange iterations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Replica exchange subroutine
!after performing a monte carlo simulation
!this subroutine loops over pairs of replicas
!with adjacent linking numbers and tests for
!exchange of their linking numbers


subroutine replicaEXCHANGE(wlc_p,wlc_d)
  use params
  use mersenne_twister

  implicit none

  !Global structures for wlc params and data
  type(wlcsim_params) :: wlc_p
  type(wlcsim_data) :: wlc_d

  !Variables for replica exchange
  integer rep             !replica index
  real(dp) deEXCHANGE     !change in energy for exchange for linking number
  real urand(1)       !random number for exchange
  real(dp) prob               !probability of exchange
  integer tempLK


  !Loop over every other replica and check for it to exchange with the one above it
  !For the replica exchange, each replica keeps its polymer conformation, and exchanges
  !linking number

  do rep = wlc_d%replicaSTART, wlc_d%replicaEND,2

     !Update number of trial exchanges
     wlc_d%nTRIALup(rep) = wlc_d%nTRIALup(rep) + 1
     wlc_d%nTRIALdown(rep+1) = wlc_d%nTRIALdown(rep+1) + 1

     deEXCHANGE = 0.0_dp

     !Include energy from original twist energy of each replica
     deEXCHANGE = -wlc_d%eelasREPLICAS(rep,4) - wlc_d%eelasREPLICAS(rep+1,4)

     !Add change in twist energy due to rep taking LK from above
     deEXCHANGE = deEXCHANGE +((2.0_dp*pi*(wlc_d%LKs(rep + 1)-wlc_d%Wrs(rep)))**2.0_dp)*wlc_p%lt/(2.0_dp*wlc_p%l)

     !Add change in twist energy due to rep + 1 takign LK from below
     deEXCHANGE = deEXCHANGE +((2.0_dp*pi*(wlc_d%LKs(rep)-wlc_d%Wrs(rep+1)))**2.0_dp)*wlc_p%lt/(2.0_dp*wlc_p%l)

     !Generate a random number and test for exchange
     call random_number(urand,wlc_d%rand_stat)
     prob = exp(-deEXCHANGE)

     !Exchange linking numbers according to Metropolis-Hastings criterion
     if (urand(1).le.prob) then
        tempLK = wlc_d%LKs(rep)
        wlc_d%LKs(rep) = wlc_d%LKs(rep + 1)
        wlc_d%LKs(rep + 1) = tempLK
        tempLK = wlc_d%nodeNUMBER(rep)
        wlc_d%nodeNUMBER (rep) = wlc_d%nodeNUMBER (rep+1)
        wlc_d%nodeNUMBER(rep+1) = tempLK

        !update number of swaps
        wlc_d%nSWAPup(rep) = wlc_d%nSWAPup(rep) + 1
        wlc_d%nSWAPdown(rep+1) = wlc_d%nSWAPdown(rep+1) + 1

     endif
  enddo

  !Alternatve value of replica start and replica end
  if (wlc_d%replicaSTART.eq.1) then
     wlc_d%replicaSTART = 2
  else
     wlc_d%replicaSTART = 1
  endif

  if (wlc_d%replicaEND.eq.wlc_d%nLKs-1) then
     wlc_d%replicaEND = wlc_d%nLKs-2
  else
     wlc_d%replicaEND = wlc_d%nLKs -1
  endif

end subroutine replicaEXCHANGE


!wlcsim subroutine for performing monte carlo simulation
subroutine wlcsim_brad(save_ind,wlc_p,wlc_d)

  use params
  use mersenne_twister
  use mpi

  implicit none

  !Structures for simulation
  type(wlcsim_params) :: wlc_p
  type(wlcsim_data) :: wlc_d
  integer, intent(in) :: save_ind

  !MPI status variables
  integer (kind = 4) status(MPI_STATUS_SIZE)
  integer (kind = 4) error
  !Integers for commuinication among replicas
  integer rep
  integer dest
  integer source
  integer LK
  real(dp) Wr
  real(dp) eelas(4)
  integer repSTART,repEND
  real(dp) deEXCHANGE
  !Counting variables
  integer i,repIND

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Run simulation
  !Cases:
  !parallel tempering on
  !parallel tempering off
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !If parallel tempering is on
  if (wlc_p%ptON) then

     !Cases:
     !Head node (id = 0)
     !Worker node (id /= 0)

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !Head node
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     if (wlc_d%id.eq.0) then

        !During replica exchange, replicas alternate between switching with
        !replica above and replica below. Initially, set the replica looking above
        !to the be first one
        repSTART = 1
        repEND = wlc_d%nLKs - 1

        !Begin replica exchange loop
        do repIND = 1,wlc_p%nReplicaExchangePerSavePoint

           !Send out Lk to the worker nodes to begin their simulations
           do rep = 1,wlc_d%nLKs
              LK = wlc_d%LKs(rep)
              dest = wlc_d%nodeNUMBER(rep)
              call MPI_Send (LK,1, MPI_INTEGER, dest,   0, &
                   MPI_COMM_WORLD,error )
           enddo

           !LK, Wr, and EELAS from the worker nodes to perform
           !replica exchange

           do rep = 1,wlc_d%nLKs
              LK = wlc_d%LKs(rep)
              source = wlc_d%nodeNUMBER(rep)
              call MPI_Recv(Wr,1, MPI_DOUBLE_PRECISION, source, 0, &
                   MPI_COMM_WORLD,status,error )
              call MPI_Recv(eelas,4, MPI_DOUBLE_PRECISION, source, 0, &
                   MPI_COMM_WORLD,status,error )
              wlc_d%Wrs(rep) = Wr
              wlc_d%eelasREPLICAS(rep,:) = eelas
           enddo

           !peform replica exchange
           call replicaEXCHANGE(wlc_p,wlc_d)


        enddo

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !Worker nodes
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     else
        !Begin replica exchange loop
        do repIND = 1,wlc_p%nREPLICAexchangePerSavePoint
           !Receive LK from head node to begin simulation
           source = 0
           call MPI_Recv(LK,1, MPI_INTEGER, source,   0, &
                MPI_COMM_WORLD,status,error )
           wlc_p%LK = LK

           !Run a monte carlo simulation for NPT
           call MCsim(wlc_d%R,wlc_d%U,wlc_p%nT,wlc_p%nB,wlc_p%nP,wlc_p%NPT,wlc_p%INTERP_BEAD_LENNARD_JONES,wlc_p%eb,&
                wlc_p%epar,wlc_p%eperp,wlc_p%gam,wlc_p%eta, wlc_p%lhc,wlc_p%vhc,wlc_p%lbox, &
                wlc_d%MCAMP,wlc_d%SUCCESS,wlc_d%successTOTAL,wlc_p%MOVEON,&
                wlc_d%WINDOW,wlc_p%RING,wlc_p%TWIST,wlc_p%Lk,wlc_p%LT,wlc_p%LP,wlc_p%L,wlc_d%rand_stat)


           !Recalculate structural quantities and energies
           call writhe(wlc_d%R,wlc_p%nB, wlc_d%Wr)
           call energy_elas(wlc_d%eelas,wlc_d%R,wlc_d%U,wlc_p%nT,wlc_p%nB,wlc_p%nP,wlc_p%eb,wlc_p%epar, &
                wlc_p%eperp,wlc_p%gam,wlc_p%eta,wlc_p%RING,wlc_p%twist,wlc_p%LK,wlc_p%lt,wlc_p%lp,wlc_p%l)


           !Communicate with the head node for replica exchange

           !Send back writhe and elastic energy back to the head node
           call MPI_SEND(wlc_d%Wr,1, MPI_DOUBLE_PRECISION, source, 0, &
                MPI_COMM_WORLD,error )
           call MPI_SEND(wlc_d%eelas,4, MPI_DOUBLE_PRECISION, source, 0, &
                MPI_COMM_WORLD,error )

        enddo
     endif

  else


     !Run a MC simulation for NSTEPS
     call MCsim(wlc_d%R,wlc_d%U,wlc_p%nT,wlc_p%nB,wlc_p%nP,wlc_p%NPT,wlc_p%INTERP_BEAD_LENNARD_JONES,wlc_p%eb,&
          wlc_p%epar,wlc_p%eperp,wlc_p%gam,wlc_p%eta, wlc_p%lhc,wlc_p%vhc,wlc_p%lbox, &
          wlc_d%MCAMP,wlc_d%SUCCESS,wlc_d%successTOTAL,wlc_p%MOVEON,&
          wlc_d%WINDOW,wlc_p%RING,wlc_p%TWIST,wlc_p%Lk,wlc_p%LT,wlc_p%LP,wlc_p%L,wlc_d%rand_stat)
  end if
end subroutine wlcsim_brad
