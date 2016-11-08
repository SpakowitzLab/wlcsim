!---------------------------------------------------------------*
      
subroutine wlcsim(rand_stat)
      
!     
!     This simulation tracks the dynamics of a single polymer
!     chain modeled as a discrete wormlike chain with bending
!     and stretching energy.
!     
!     Andrew Spakowitz
!     Written 9-2-13
!     
!     Edited by Shifan  prior to 2016
!     Edited heavily by Quinn in spring of 2016
!
!     Variables within the simulation

  use setPrecision
  use simMod
  use mersenne_twister  ! so that we know the size of rand_stat

  IMPLICIT NONE
  type(random_stat), intent(inout) :: rand_stat ! state of random number generator

  ! miscellaneous
  integer I            
  character*4 fileind       ! index of output
  character*16 iostr       ! file for output

!     Simulation input variables
  
  integer INTON             ! Include polymer interactions

! simulation data strucutres
  type(MCvar) mc
  type(MCData) md

!-------------------------------------------------------
!
!    Set simulation parameers 
!
!--------------------------------------------------------
  iostr='input/params'
  print*, "setting parameters from: ", iostr
  call MCvar_setparams(mc,iostr)
  call MCvar_allocate(mc,md)
 
!-----------------------------------------------------
!
!   Initial condition / Restart
!
!----------------------------------------------------
  !  Calculate volume of bins
  if (mc%confinetype.eq.3) then 
      call MC_caclVolume(mc%confinetype,mc%NBINX,mc%dbin, mc%lbox(1), &
                         md%Vol,rand_stat)  ! calculate partial volumes
  else
      do I=1,mc%NBIN
           md%Vol(I)=mc%dbin**3
      enddo
  endif

  if (mc%restart) then
      call pt_restart(mc,md)
      !print*, '-----load simulation-----'
      !iostr='BinaryfileName'
      !stop 1
      !call MCvar_readBindary(mc,md,iostr)
  else
    
      
    !   Setup the initial condition
      call initcond(md%R,md%U,md%AB,mc%NT,mc%NB,mc%NP,mc%FRMfile,mc%para,mc%lbox, &
                    mc%settype,rand_stat)
    
    !   Load in AB sequence
      if (mc%FRMCHEM) then
          iostr='input/ab'
          call MCvar_loadAB(mc,md,iostr)
      else
          call initchem(md%AB,mc%NT,mc%N,mc%G,mc%NP,mc%FA,mc%LAM,rand_stat)
      endif
    
      
    !   Load methalation sequence
      if (mc%FRMMETH) then
          ! more to come here ...
          print*, "wlcsim: FRMMETH not fininshed"
          stop 1
      else
          call initchem(md%METH,mc%NT,mc%N,mc%G,mc%NP,mc%F_METH,mc%LAM_METH,rand_stat)        
      endif
    
    ! Load External field
      if (mc%FRMField) then
          iostr='input/h_A'
          Call MCvar_LoadField(mc,md,iostr)
      else
          Call MCvar_MakeField(mc,md)
      endif

      ! Get assignement from head node
      call PT_override(mc,md)

      iostr='data/r0'
      I=0;
      call MCvar_saveR(mc,md,iostr,0)
      
      iostr='data/params'
      call MCvar_saveparameters(mc,iostr)
    
      iostr='data/u0'
      call MCvar_saveU(mc,md,iostr)

      mc%ind=1
  endif

 ! call MCvar_printDescription(mc)  ! output simulation parameters

!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!              Begin simulation
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  print*, 'Beginning simulation: rep', mc%rep, " id", mc%id

  do WHILE ((mc%ind).LE.mc%indMAX) 

     ! for changing constants durring run
     if (mc%useSchedule) then
         call strength_schedule(mc,INTON)
     endif
     

!   * Perform a MC simulation *
    call MCsim(mc,md,mc%NSTEP,INTON,rand_stat)

!    Save the conformation and the metrics
    write (fileind,'(I4)'), mc%ind

    !Save various energy contiributions to file 
    iostr='data/out1'
    call MCvar_appendEnergyData(mc,iostr)
     
    !part 2.5 - adaptations
    iostr='data/out3'
    call MCvar_appendAdaptData(mc,iostr)

    if (mc%savePhi) then
        write(iostr,"(I6)"), mc%ind
        iostr='data/phi' // trim(adjustL(iostr))
        call MCVar_savePHI(mc,md,iostr)    
    endif

    write(iostr,"(I6)"), mc%ind
    iostr='data/r' // trim(adjustL(iostr))
    call MCvar_saveR(mc,md,iostr,0)
   
    if (mc%saveU) then
        write(iostr,"(I6)"), mc%ind
        iostr='data/u' // trim(adjustL(iostr))
        call MCvar_saveU(mc,md,iostr)
    endif


    print*, '________________________________________'
    print*, 'Time point ',mc%ind, ' out of', mc%indMAX
    call MCvar_printEnergies(mc)
    call MCvar_printWindowStats(mc)
    !call MCvar_printPhi(mc,md)
    mc%ind=mc%ind+1    
  enddo
  
end
!---------------------------------------------------------------*
