!---------------------------------------------------------------*
      
Subroutine wlcsim(rand_stat)
      
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
  INTEGER I            
  character*4 fileind       ! Index of output
  character*16 iostr       ! File for output

!     Simulation input variables
  
  INTEGER INTON             ! Include polymer interactions

! simulation data strucutres
  TYPE(MCvar) mc
  TYPE(MCData) md

!-------------------------------------------------------
!
!    Set simulation parameers 
!
!--------------------------------------------------------
  iostr='input/params'
  print*, "setting parameters from: ", iostr
  call MCvar_setParams(mc,iostr)
  call MCvar_allocate(mc,md)
 
!-----------------------------------------------------
!
!   Initial condition / Restart
!
!----------------------------------------------------
  !  Calculate volume of bins
  if (mc%confineType.eq.3) then 
      call MC_caclVolume(mc%confineType,mc%NBINX,mc%DEL, mc%LBox(1), &
                         md%Vol,rand_stat)  ! calculate partial volumes
  else
      do I=1,mc%NBIN
           md%Vol(I)=mc%del**3
      enddo
  endif

  if (mc%restart) then
      call pt_restart(mc,md)
      !PRINT*, '-----load simulation-----'
      !iostr='BinaryFileName'
      !stop 1
      !call MCvar_readBindary(mc,md,iostr)
  else
    
      
    !   Setup the initial condition
      call initcond(md%R,md%U,md%AB,mc%NT,mc%NB,mc%NP,mc%FRMFILE,mc%PARA,mc%LBOX, &
                    mc%setType,rand_stat)
    
    !   Load in AB sequence
      IF (mc%FRMCHEM) THEN
          iostr='input/ab'
          call MCvar_loadAB(mc,md,iostr)
      ELSE
          call initchem(md%AB,mc%NT,mc%N,mc%G,mc%NP,mc%FA,mc%LAM,rand_stat)
      ENDIF
    
      
    !   Load methalation sequence
      IF (mc%FRMMETH) THEN
          ! more to come here ...
          print*, "wlcsim: FRMMETH not fininshed"
          stop 1
      ELSE
          call initchem(md%METH,mc%NT,mc%N,mc%G,mc%NP,mc%F_METH,mc%LAM_METH,rand_stat)        
      ENDIF
    
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
      call MCvar_saveParameters(mc,iostr)
    
      iostr='data/u0'
      call MCvar_saveU(mc,md,iostr)

      mc%IND=1
  endif

 ! call MCvar_printDescription(mc)  ! output simulation parameters

!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!              Begin simulation
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  print*, 'Beginning simulation: rep', mc%rep, " id", mc%id

  DO WHILE ((mc%IND).LE.mc%INDMAX) 

     ! for changing constants durring run
     if (mc%UseSchedule) then
         call strength_schedule(mc,INTON)
     endif
     

!   * Perform a MC simulation *
    call MCsim(mc,md,mc%NSTEP,INTON,rand_stat)

!    Save the conformation and the metrics
    write (fileind,'(I4)'), mc%IND

    !Save various energy contiributions to file 
    iostr='data/out1'
    call MCvar_appendEnergyData(mc,iostr)
     
    !part 2.5 - adaptations
    iostr='data/out3'
    call MCvar_appendAdaptData(mc,iostr)

    if (mc%savePhi) then
        write(iostr,"(I6)"), mc%IND
        iostr='data/phi' // trim(adjustL(iostr))
        call MCVar_savePHI(mc,md,iostr)    
    endif

    write(iostr,"(I6)"), mc%IND
    iostr='data/r' // trim(adjustL(iostr))
    call MCvar_saveR(mc,md,iostr,0)
   
    if (mc%saveU) then
        write(iostr,"(I6)"), mc%IND
        iostr='data/u' // trim(adjustL(iostr))
        call MCvar_saveU(mc,md,iostr)
    endif


    PRINT*, '________________________________________'
    PRINT*, 'Time point ',mc%IND, ' out of', mc%INDMAX
    call MCvar_printEnergies(mc)
    call MCvar_printWindowStats(mc)
    !call MCvar_printPhi(mc,md)
    mc%IND=mc%IND+1    
  ENDDO
  
END
!---------------------------------------------------------------*
