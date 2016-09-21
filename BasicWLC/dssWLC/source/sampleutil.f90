MODULE SAMPLEUTIL
  ! Utilities for sampling from various distributions, 
  ! including rejection sampling for segment joints of the dssWLC
  USE MT19937, ONLY : GRND, RNORM,MVNORM
  use QUATUTIL, ONLY : PI
  
  IMPLICIT NONE

CONTAINS
  SUBROUTINE GETEQUILCHAIN(CHAINP,TYPESAMPLE,LASTCOORDS,STARTCOORDS)
    ! sample a configuration of the chain from the equilibrium free distribution
    ! uses rejection sampling for the coupled coordinates
    ! should be exact for any segment length (no gaussian approximations
    ! TYPESAMPLE is the type of sampling to do:
    ! 1) old version of rejection sampling (better for very flexible chains)
    ! 2) multivariate normal version of rejection sampling (better for very stiff/short segments)
    ! 3) monte carlo sampling
    ! LASTCOORDS: if doing monte carlo returns last coordinate and range values
    ! if STARTCOORDS is supplied (and using MC) then use the given starting 
    ! coordinates and step ranges, without an initialization period

    USE CHAINUTIL, ONLY : CHAIN
    USE GENUTIL, ONLY : CROSS_PRODUCT, NORMALIZE
    USE KEYS, ONLY : MCSTATSTEPS, MCINITSTEPS, GAUSSIANCHAIN, LOGRTERM, NEDGESEG

    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    INTEGER, INTENT(IN) :: TYPESAMPLE
    DOUBLE PRECISION, INTENT(OUT) :: LASTCOORDS(6)
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: STARTCOORDS(6)
    DOUBLE PRECISION :: COORDS(4,CHAINP%NPT-1), Z, R, RHO, PHI, ST, PHIPERP,PHIU
    DOUBLE PRECISION :: DEL(CHAINP%NPT), ETA(CHAINP%NPT)
    INTEGER :: BC
    DOUBLE PRECISION :: XAX(3), YAX(3), FINALRANGES(3)

    IF (TYPESAMPLE.LT.1.OR.TYPESAMPLE.GT.3) THEN
       PRINT*, 'INVALID VALUE OF TYPESAMPLE', TYPESAMPLE
       STOP 1
    ENDIF

    CHAINP%POS(:,1) = (/0D0,0D0,0D0/)
    ! place first orientation vector uniformly
    RHO = GRND()*2-1; ST = SQRT(1-RHO**2)
    PHI = GRND()*2*PI
    CHAINP%UVEC(:,1) = (/ST*COS(PHI),ST*SIN(PHI),RHO/)

    ! Get relative coordinates for each bead (relative to previous one)
    ! coordinates are (in order):
    ! z = translation along uvec
    ! r = magnitude translation perpendicular to uvec
    ! rho = cos(angle between subsequent uvecs)
    ! phi = angle between the perp translation and the subsequent uvec

    DEL = CHAINP%LS
    ETA = -CHAINP%EC/CHAINP%LP
    IF (TYPESAMPLE.EQ.3) THEN
       IF (NEDGESEG.GT.0) THEN
          PRINT*, 'MONTE CARLO INITIAL SAMPLING NOT YET SET UP WITH DIFFERENT EDGE SEGMENTS!'
          STOP 1
       ENDIF
       !monte carlo sampling
       IF (PRESENT(STARTCOORDS)) THEN
          CALL  SAMPLERELCOORDSMC(DEL(1),CHAINP%LP(1),CHAINP%GAM(1),&
               & CHAINP%EPAR(1),CHAINP%EPERP(1),ETA(1),CHAINP%NPT-1,&
               & MCSTATSTEPS,0,COORDS,FINALRANGES,STARTCOORDS)
       ELSE
          CALL  SAMPLERELCOORDSMC(DEL(1),CHAINP%LP(1),CHAINP%GAM(1),&
               & CHAINP%EPAR(1),CHAINP%EPERP(1),ETA(1),CHAINP%NPT-1,&
               & MCSTATSTEPS,MCINITSTEPS,COORDS,FINALRANGES)
       ENDIF
       LASTCOORDS(1:3) = COORDS(2:4,CHAINP%NPT-1)
       LASTCOORDS(4:6) = FINALRANGES
    ELSE
       IF (GAUSSIANCHAIN) THEN          
          ! plain old normal sampling for gaussian chain
          IF (NEDGESEG.GT.0) THEN
             CALL SAMPLERELCOORDSGAUSSIAN(DEL(1),CHAINP%EPAR(1),CHAINP%EPERP(1),NEDGESEG,COORDS(:,1:NEDGESEG))
             CALL SAMPLERELCOORDSGAUSSIAN(DEL(CHAINP%NPT-1),&
                  & CHAINP%EPAR(CHAINP%NPT-1),CHAINP%EPERP(CHAINP%NPT-1),&
                  & NEDGESEG,COORDS(:,CHAINP%NPT-NEDGESEG:CHAINP%NPT-1))
             CALL SAMPLERELCOORDSGAUSSIAN(DEL(NEDGESEG+1),&
                  & CHAINP%EPAR(NEDGESEG+1),CHAINP%EPERP(NEDGESEG+1),&
                  & CHAINP%NPT-1-2*NEDGESEG,COORDS(:,NEDGESEG+1:CHAINP%NPT-NEDGESEG-1))
          ELSE
             CALL SAMPLERELCOORDSGAUSSIAN(DEL(1),CHAINP%EPAR(1),CHAINP%EPERP(1),CHAINP%NPT-1,COORDS)
          ENDIF

       ELSE IF (CHAINP%SHEARABLE) THEN
          ! use rejection sampling
          IF (TYPESAMPLE.EQ.1) THEN
             IF (NEDGESEG.GT.0) THEN
                CALL SAMPLERELCOORDS(DEL(1),CHAINP%LP(1),CHAINP%GAM(1),&
                     & CHAINP%EPAR(1),CHAINP%EPERP(1),ETA(1),NEDGESEG,&
                     & COORDS(:,1:NEDGESEG),CHAINP%SHEARABLE)
                CALL SAMPLERELCOORDS(DEL(CHAINP%NPT-1),CHAINP%LP(CHAINP%NPT-1),CHAINP%GAM(CHAINP%NPT-1),&
                     & CHAINP%EPAR(CHAINP%NPT-1),CHAINP%EPERP(CHAINP%NPT-1),ETA(CHAINP%NPT-1),NEDGESEG,&
                     & COORDS(:,CHAINP%NPT-NEDGESEG:CHAINP%NPT-1),CHAINP%SHEARABLE)
                CALL SAMPLERELCOORDS(DEL(NEDGESEG+1),CHAINP%LP(NEDGESEG+1),CHAINP%GAM(NEDGESEG+1),&
                     & CHAINP%EPAR(NEDGESEG+1),CHAINP%EPERP(NEDGESEG+1),ETA(NEDGESEG+1),CHAINP%NPT-1-2*NEDGESEG,&
                     & COORDS(:,NEDGESEG+1:CHAINP%NPT-NEDGESEG-1),CHAINP%SHEARABLE)
             ELSE
                CALL SAMPLERELCOORDS(DEL(1),CHAINP%LP(1),CHAINP%GAM(1),&
                     & CHAINP%EPAR(1),CHAINP%EPERP(1),ETA(1),CHAINP%NPT-1,COORDS,CHAINP%SHEARABLE)
             ENDIF

          ELSEIF (TYPESAMPLE.EQ.2) THEN
             ! use multivariate normal rejection sampling
              IF (NEDGESEG.GT.0) THEN
                CALL SAMPLERELCOORDSMVN(DEL(1),CHAINP%LP(1),CHAINP%GAM(1),&
                     & CHAINP%EPAR(1),CHAINP%EPERP(1),ETA(1),NEDGESEG,&
                     & COORDS(:,1:NEDGESEG))
                CALL SAMPLERELCOORDSMVN(DEL(CHAINP%NPT-1),CHAINP%LP(CHAINP%NPT-1),CHAINP%GAM(CHAINP%NPT-1),&
                     & CHAINP%EPAR(CHAINP%NPT-1),CHAINP%EPERP(CHAINP%NPT-1),ETA(CHAINP%NPT-1),NEDGESEG,&
                     & COORDS(:,CHAINP%NPT-NEDGESEG:CHAINP%NPT-1))
                CALL SAMPLERELCOORDSMVN(DEL(NEDGESEG+1),CHAINP%LP(NEDGESEG+1),CHAINP%GAM(NEDGESEG+1),&
                     & CHAINP%EPAR(NEDGESEG+1),CHAINP%EPERP(NEDGESEG+1),ETA(NEDGESEG+1),CHAINP%NPT-1-2*NEDGESEG,&
                     & COORDS(:,NEDGESEG+1:CHAINP%NPT-NEDGESEG-1))
             ELSE
                CALL SAMPLERELCOORDSMVN(DEL(1),CHAINP%LP(1),CHAINP%GAM(1),&
                     & CHAINP%EPAR(1),CHAINP%EPERP(1),ETA(1),CHAINP%NPT-1,COORDS)
             ENDIF
          ENDIF
       ELSE
          IF (LOGRTERM) THEN
             IF (NEDGESEG.GT.0) THEN
                CALL SAMPLERELCOORDS(DEL(1),CHAINP%LP(1),CHAINP%GAM(1),&
                     & CHAINP%EPAR(1),CHAINP%EPERP(1),ETA(1),NEDGESEG,&
                     & COORDS(:,1:NEDGESEG),CHAINP%SHEARABLE)
                CALL SAMPLERELCOORDS(DEL(CHAINP%NPT-1),CHAINP%LP(CHAINP%NPT-1),CHAINP%GAM(CHAINP%NPT-1),&
                     & CHAINP%EPAR(CHAINP%NPT-1),CHAINP%EPERP(CHAINP%NPT-1),ETA(CHAINP%NPT-1),NEDGESEG,&
                     & COORDS(:,CHAINP%NPT-NEDGESEG:CHAINP%NPT-1),CHAINP%SHEARABLE)
                CALL SAMPLERELCOORDS(DEL(NEDGESEG+1),CHAINP%LP(NEDGESEG+1),CHAINP%GAM(NEDGESEG+1),&
                     & CHAINP%EPAR(NEDGESEG+1),CHAINP%EPERP(NEDGESEG+1),ETA(NEDGESEG+1),CHAINP%NPT-1-2*NEDGESEG,&
                     & COORDS(:,NEDGESEG+1:CHAINP%NPT-NEDGESEG-1),CHAINP%SHEARABLE)
             ELSE
                CALL SAMPLERELCOORDS(DEL(1),CHAINP%LP(1),CHAINP%GAM(1),&
                     & CHAINP%EPAR(1),CHAINP%EPERP(1),ETA(1),CHAINP%NPT-1,COORDS,CHAINP%SHEARABLE)
             ENDIF            
          ELSE
             IF (NEDGESEG.GT.0) THEN
                CALL SAMPLERELCOORDSNOSHEAR(DEL(1),CHAINP%LP(1),CHAINP%GAM(1),&
                     & CHAINP%EPAR(1),NEDGESEG,COORDS(:,1:NEDGESEG),CHAINP%STRETCHABLE)
                CALL SAMPLERELCOORDSNOSHEAR(DEL(CHAINP%NPT-1),CHAINP%LP(CHAINP%NPT-1),CHAINP%GAM(CHAINP%NPT-1),&
                     & CHAINP%EPAR(CHAINP%NPT-1),NEDGESEG,COORDS(:,CHAINP%NPT-NEDGESEG:CHAINP%NPT-1),CHAINP%STRETCHABLE)
                CALL SAMPLERELCOORDSNOSHEAR(DEL(NEDGESEG+1),CHAINP%LP(NEDGESEG+1),CHAINP%GAM(NEDGESEG+1),&
                     & CHAINP%EPAR(NEDGESEG+1),CHAINP%NPT-1-2*NEDGESEG,&
                     & COORDS(:,NEDGESEG+1:CHAINP%NPT-1-NEDGESEG),CHAINP%STRETCHABLE)
             ELSE
                CALL SAMPLERELCOORDSNOSHEAR(DEL(1),CHAINP%LP(1),CHAINP%GAM(1),&
                     & CHAINP%EPAR(1),CHAINP%NPT-1,COORDS,CHAINP%STRETCHABLE)
             ENDIF
          ENDIF
       ENDIF
    ENDIF

    DO BC = 2,CHAINP%NPT
       Z = COORDS(1,BC-1); R = COORDS(2,BC-1); 
       RHO = COORDS(3,BC-1); PHI = COORDS(4,BC-1)
       ST = SQRT(1-RHO**2)

       ! set up orthonormal triad for previous bead
       IF (CHAINP%UVEC(2,BC-1)**2 + CHAINP%UVEC(3,BC-1)**2.EQ.0) THEN
          YAX = (/0D0,0D0,1D0/)
          XAX = (/0D0,1D0,0D0/)
       ELSE
          CALL CROSS_PRODUCT(CHAINP%UVEC(:,BC-1),(/1D0,0D0,0D0/),YAX)
          CALL NORMALIZE(YAX)
          CALL CROSS_PRODUCT(YAX,CHAINP%UVEC(:,BC-1),XAX)
          CALL NORMALIZE(XAX)
       ENDIF
       
       IF (CHAINP%SHEARABLE) THEN
          ! pick orientation of perpendicular displacement randomly around uvec
          PHIPERP = GRND()*2*PI
          ! place bead
          !       print*, 'testx1:', bc, z, r, rho, phi, sqrt(sum(chainp%uvec(:,bc-1)**2))
          CHAINP%POS(:,BC) = CHAINP%POS(:,BC-1) + Z*CHAINP%UVEC(:,BC-1) &
               & + R*COS(PHIPERP)*XAX + R*SIN(PHIPERP)*YAX
       ELSE
          PHIPERP = 0D0
          CHAINP%POS(:,BC) = CHAINP%POS(:,BC-1) + Z*CHAINP%UVEC(:,BC-1)
       ENDIF

       ! place next orientation vector
       PHIU = PHIPERP+PHI

       CHAINP%UVEC(:,BC) = ST*COS(PHIU)*XAX + ST*SIN(PHIU)*YAX + RHO*CHAINP%UVEC(:,BC-1)
    ENDDO
   
  END SUBROUTINE GETEQUILCHAIN
    
  SUBROUTINE A0BOUNDFUNC(R,PARAM,FU,DU)
    ! function used to find enveloping Lorentzian for cylindrical gaussian
    ! PARAM is array of parameters (A,B,R0,M0)
    ! FU returns function values, DU returns derivative
    ! sign flipped upside down to do minimization rather than maximization
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: R, PARAM(:)
    DOUBLE PRECISION, INTENT(OUT) :: FU, DU
    DOUBLE PRECISION :: A, B, R0, M0

    A = PARAM(1); B = PARAM(2); R0 = PARAM(3); M0 = PARAM(4)
     

    IF (ABS(R-R0).LT.1D-6) THEN
       ! use asymptotic form near r0       
       !DU =  (-B+SQRT(2/A+B**2))/(6+3*A*B**2)   
       !FU = 1 + B*SQRT(A/(2+A*B**2))/(2*A) + DU*(R-R0)
       DU = (-2*A*B-A**2*B**3+2*SQRT(A*(2+A*B**2))+B**2*SQRT(A**3*(2+A*B**2)))/(3*A*(2+A*B**2)**2)
       FU = (2*A*B+A**2*B**3+2*SQRT(A*(2+A*B**2))+B**2*SQRT(A**3*(2+A*B**2)))/2/SQRT((A*(2+A*B**2))**3) + DU*(R-R0)
       !PRINT*, 'TESTX3', DU
    ELSE
       DU = (R-R0)*(-2*R**2-EXP(A*(B-R)**2)*M0*(R*(-3-2*A*(B-R)*(R-R0))+R0)) / &
            & (R-M0*EXP(A*(B-R)**2))**2
       FU = R*(R-R0)**2/(M0*EXP(A*(R-B)**2) - R)
    END IF
    
    FU = -FU; DU = -DU;

  !  PRINT*, 'TESTX1:', R, PARAM, FU, DU
  END SUBROUTINE A0BOUNDFUNC

  SUBROUTINE GETLORENTZENVELOPE(A,B,LORPARAM)
    ! find an enveloping Lorenzian distribution
    ! for the cylindrical normal P~r*exp(-a*(r-b)^2)
    ! returns parameters: r0,m0,a0
    ! enveloping distribution g ~ m0/(1+(r-r0)^2/a0)
    ! see 4/24/2013 notes
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: A,B
    DOUBLE PRECISION, INTENT(OUT) :: LORPARAM(3)
    DOUBLE PRECISION :: R0,M0,A0
    DOUBLE PRECISION :: FU, DU, PARAM(4)
    DOUBLE PRECISION :: AX, BX, STEPSIZE, CX, FMIN, XMIN
    ! minimization function
    DOUBLE PRECISION :: DBRENT
    EXTERNAL DBRENT
    ! max number of steps to try when bracketing
    INTEGER, PARAMETER :: MAXNSTEP = 100; 
    INTEGER :: C


    ! position of maximum likelihood
    R0 = B/2 + 0.5D0*SQRT(2/A+B**2)
    ! maximum value (at r0)
    M0 = R0*EXP(-A*(R0-B)**2)
    
    ! find an appropriate width of the Lorentzian such that it's always greater than our cylindrical Gaussian
    ! WARNING: technically this finds a local rather than global minimum, though for this particular functional form it doesn't seem to be a problem

!    PRINT*, R0, M0

    ! bracket the minimum
    AX = 0D0
    BX = R0
    STEPSIZE = R0
    DO C = 1,MAXNSTEP
       CX = R0+C*STEPSIZE
       CALL A0BOUNDFUNC(CX,(/A,B,R0,M0/),FU,DU)
       IF (FU.GT.-M0) THEN
          EXIT
       ENDIF
    ENDDO
    IF (C.GE.MAXNSTEP) THEN
       PRINT*, 'ERROR IN SAMPLECYLNORMAL: failed to bracket extremum'
       STOP 1
    ENDIF

    PARAM = (/A,B,R0,M0/)

    ! PRINT*, 'BRACKETS:'
    ! CALL A0BOUNDFUNC(AX,PARAM,FU,DU)
    ! PRINT*, AX, FU, DU
    ! CALL A0BOUNDFUNC(BX,PARAM,FU,DU)
    ! PRINT*, BX, FU, DU
    ! CALL A0BOUNDFUNC(CX,PARAM,FU,DU)
    ! PRINT*, CX, FU, DU
    
    FMIN = DBRENT(AX,BX,CX,A0BOUNDFUNC,4,PARAM,SQRT(EPSILON(1D0)),XMIN)

   ! PRINT*, 'MINIMUM:', XMIN, FMIN
    A0 = -FMIN
    
    LORPARAM = (/R0,M0,A0/)
  END SUBROUTINE GETLORENTZENVELOPE

  SUBROUTINE SAMPLECYLNORMAL(A,B,LORPARAM,RVAL,NTRY)
    ! sample from a cylindrical normal distribution P ~ r exp(-a*(r-b)^2) for r>0
    ! uses rejection sampling with a Lorentzian (Cauchy) distribution envelope
    ! defined by g ~ M0/(1+(r-r0)^2/a0) where lorparam = (r0,m0,a0)
    ! get the enveloping distribution parameters using GETLORENTZENVELOPE
    ! see 4/24/2013 notes
    ! NTRY is the number of tries required to get an accepted sample

    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: A,B, LORPARAM(3)
    DOUBLE PRECISION, INTENT(OUT) :: RVAL
    INTEGER, INTENT(OUT) :: NTRY
    DOUBLE PRECISION :: R0,M0,A0, R, U, U2, RATIO
    INTEGER :: TRY
    ! maximum number of rejection trials before giving up
    INTEGER, PARAMETER :: MAXNTRY = 1000 

    R0 = LORPARAM(1); M0 = LORPARAM(2); A0 = LORPARAM(3);
    !PRINT*, 'TESTX0:', R0, M0, A0
    DO TRY = 1,MAXNTRY
       ! Sample R from the lorentzian
       U =  GRND() ! uniform variate
       R = R0 - SQRT(A0)/TAN(PI*U)
       !PRINT*, 'TESTX1:', U, PI*U, R
       
       ! uniform sample to decide whether to reject
       U2 = GRND() 
       RATIO = R*EXP(-A*(R-B)**2)/(M0/(1+(R-R0)**2/A0))
       IF (U2.LT.RATIO) THEN
          EXIT
       ENDIF
    ENDDO

    IF (TRY.GE.MAXNTRY) THEN
       PRINT*, 'ERROR IN SAMPLECYLNORMAL: failed to generate successful trial'
       STOP 1
    ENDIF

    RVAL = R; NTRY = TRY;
    !PRINT*, R, TRY
  END SUBROUTINE SAMPLECYLNORMAL

  SUBROUTINE SAMPLERRHOPHI(DEL,EB,EPERPH,ETA,R,RHO,PHI,NTRIAL)
    ! rejection sampling for R,Rho, PHI coordinates
    ! enveloping distribution is one with uniform phi
    ! sample R and RHO from the enveloping distribution that is independent of PHI
    ! EPERPH corresponds to eperp_hat in the notes, but EPERP in the code
    ! ETA is -EC/LP in the simulation code and EB=LP
    ! see notes from 4/24/2013   
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: DEL, EB, EPERPH, ETA
    DOUBLE PRECISION, INTENT(OUT) :: R, RHO, PHI
    INTEGER, INTENT(OUT) :: NTRIAL
    DOUBLE PRECISION :: A,B, C, U, LORPARAM(3)
    INTEGER :: TRY, NTRYR
    DOUBLE PRECISION :: SP, RATIO
    ! maximum rejection trials before giving up
    INTEGER, PARAMETER :: MAXNTRY = 1000

    A = EPERPH/2/DEL
    B = ETA*EB/EPERPH
    C = EB/2/DEL

    CALL GETLORENTZENVELOPE(A,B,LORPARAM)

    DO TRY = 1,MAXNTRY
       ! sample rho from a truncated exponential
       U = GRND()
       RHO = 1/C*LOG(EXP(-C)+2*U*SINH(C))

       ! sample phi uniformly
       PHI = GRND()*2*PI
       SP = cos(PHI)

       ! sample R from cylindrical normal
       CALL SAMPLECYLNORMAL(A,B,LORPARAM,R,NTRYR)
       
       RATIO = EXP(ETA*EB/DEL*R*(SP*SQRT(1-RHO**2)-1))

       U= GRND()
       IF (U.LT.RATIO) EXIT          
    END DO

    IF (TRY.GE.MAXNTRY) THEN
       PRINT*, 'ERROR IN SAMPLERRHOPHI: failed to generate accepted sample'
       stop 1
    ENDIF

    NTRIAL = TRY
  END SUBROUTINE SAMPLERRHOPHI

  SUBROUTINE SAMPLERELCOORDS(DEL,EB,GAM,EPAR,EPERPH,ETA,NSAMP,COORDS,SHEARABLE)
    ! sample relative coordinates Z, R, RHO, PHI for segment junctions
    ! generates NSAMP samples and stores them in the NSAMPx4 array COORDS
    ! see notes from 4/24/2013
    ! uses rejection sampling for the coupled R,RHO,PHI coordinates
    ! EPERPH corresponds to eperp_hat in the notes, but EPERP in the simulation code
    ! ETA is -EC/LP in the simulation code and EB=LP
    ! see notes from 4/24/2013   
    ! if SHEARABLE is set to false do the infinite shear modulus limit

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: DEL, EB, GAM,EPAR,EPERPH, ETA
    INTEGER, INTENT(IN) :: NSAMP
    DOUBLE PRECISION, INTENT(OUT) :: COORDS(4,NSAMP)
    LOGICAL, INTENT(IN) :: SHEARABLE
    DOUBLE PRECISION :: A,B, C, D, U, LORPARAM(3), UN, R, RHO, PHI, Z
    INTEGER :: TRY, NTRYR, SC
    DOUBLE PRECISION :: CP, RATIO
    
    ! maximum rejection trials before giving up
    INTEGER, PARAMETER :: MAXNTRY = 100000

    ! parameters for the cylindrical normal sampling of R
    A = EPERPH/2/DEL
    B = ETA*EB/EPERPH
    ! parameter for the truncated exponential rho sampling
    C = EB/DEL
    ! parameter for the normal z sampling
    D = SQRT(DEL/EPAR)

    ! get enveloping lorentz distribution for R sampling
    CALL GETLORENTZENVELOPE(A,B,LORPARAM)

    DO SC = 1,NSAMP
       ! sample Z from a normal distribution
       UN = RNORM()
       Z = UN*D + GAM*DEL

       ! rejection sampling from the coupled r,rho,phi distribution
       DO TRY = 1,MAXNTRY
          ! sample rho from a truncated exponential
          U = GRND()
          RHO = 1/C*LOG(EXP(-C)+2*U*SINH(C))
                    
          ! sample phi uniformly
          PHI = GRND()*2*PI
          CP = COS(PHI)
          
          IF (SHEARABLE) THEN
             ! sample R from cylindrical normal          
             CALL SAMPLECYLNORMAL(A,B,LORPARAM,R,NTRYR)
          ELSE
             R = 0
          ENDIF
          
          ! decide whether to accept
          RATIO = EXP(ETA*EB/DEL*R*(CP*SQRT(1-RHO**2)-1))         
          U= GRND()
          IF (U.LT.RATIO) EXIT          
       END DO

       IF (TRY.GE.MAXNTRY) THEN
          PRINT*, 'ERROR IN SAMPLERELCOORDS: failed to generate accepted sample'
          stop 1
       ENDIF

       COORDS(:,SC) = (/Z,R,RHO,PHI/)
    ENDDO

  END SUBROUTINE SAMPLERELCOORDS

  SUBROUTINE SAMPLERELCOORDSNOSHEAR(DEL,EB,GAM,EPAR,NSAMP,COORDS,DOSTRETCH)
    ! sample relative coordinates Z, RHO for segment junctions assuming no shear
    ! generates NSAMP samples and stores them in the NSAMPx4 array COORDS
    ! uses direct normal and exponential sampling
    ! PHI is sampled uniformly

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: DEL, EB, GAM,EPAR
    INTEGER, INTENT(IN) :: NSAMP
    DOUBLE PRECISION, INTENT(OUT) :: COORDS(4,NSAMP)
    LOGICAL, INTENT(IN) :: DOSTRETCH
    DOUBLE PRECISION :: C, D, U, LORPARAM(3), UN, R, RHO, PHI, Z
    INTEGER ::  SC
    DOUBLE PRECISION :: CP, RATIO
    ! maximum rejection trials before giving up
    INTEGER, PARAMETER :: MAXNTRY = 100000

    ! parameter for the truncated exponential rho sampling
    C = EB/DEL
    ! parameter for the normal z sampling
    D = SQRT(DEL/EPAR)

    DO SC = 1,NSAMP
       IF (DOSTRETCH) THEN
          ! sample Z from a normal distribution
          UN = RNORM()
          Z = UN*D + GAM*DEL
       ELSE
          Z = GAM*DEL
       ENDIF
       
       ! sample rho from a truncated exponential
       U = GRND()
       RHO = 1/C*LOG(EXP(-C)+2*U*SINH(C))
                    
       PHI = 2*PI*GRND()
       
       COORDS(:,SC) = (/Z,0D0,RHO,PHI/)
    ENDDO

  END SUBROUTINE SAMPLERELCOORDSNOSHEAR

  SUBROUTINE SAMPLERELCOORDSGAUSSIAN(DEL,EPAR,EPERP,NSAMP,COORDS)
    ! sample relative coordinates segment junctions assuming 
    ! chain is gaussian

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: DEL, EPAR,EPERP
    INTEGER, INTENT(IN) :: NSAMP
    DOUBLE PRECISION, INTENT(OUT) :: COORDS(4,NSAMP)
    DOUBLE PRECISION :: C, D, U, LORPARAM(3), UN, R, RHO, PHI, Z, X1, X2
    INTEGER ::  SC
    DOUBLE PRECISION :: CP, RATIO
    ! maximum rejection trials before giving up
    INTEGER, PARAMETER :: MAXNTRY = 100000

    ! parameter for the truncated exponential rho sampling
    C = SQRT(DEL/EPERP)
    ! parameter for the normal z sampling
    D = SQRT(DEL/EPAR)

    DO SC = 1,NSAMP
       UN = RNORM()
       Z = UN*D
       UN = RNORM()
       X1 = UN*C
       UN = RNORM()
       X2 = UN*C
       R = SQRT(X1**2+X2**2)
       RHO = 0D0
       PHI = 0D0

       COORDS(:,SC) = (/Z,R,RHO,PHI/)
    ENDDO

  END SUBROUTINE SAMPLERELCOORDSGAUSSIAN

  SUBROUTINE SAMPLERELCOORDSMC(DEL,EB,GAM,EPAR,EPERPH,ETA,NSAMP,SAMPEVERY,&
& INITSAMP,COORDS,FINALRANGES,STARTCOORDS)
    ! sample a sequence of relative coordinates using monte carlo sampling
    ! for the coupled coordinates r,rho,phi
    ! and ordinary gaussian for the stretch coordinate z
    ! FINALRANGES is the final step sizes
    ! optionally, STARTCOORDS(1:3) gives the starting r,rho,phi coordinates
    ! and STARTCOORDS(4:6) gives the starting step ranges for r,rho,phi
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: DEL, EB, GAM,EPAR,EPERPH, ETA
    INTEGER, INTENT(IN) :: NSAMP,SAMPEVERY,INITSAMP
    DOUBLE PRECISION, INTENT(OUT) :: COORDS(4,NSAMP), FINALRANGES(3)    
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: STARTCOORDS(6)
    DOUBLE PRECISION :: D, Z, UN, FACC
    INTEGER :: SC

    ! parameter for the normal z sampling
    D = SQRT(DEL/EPAR)

    ! sample the independent stretch coordinate
    DO SC = 1,NSAMP
       ! sample Z from a normal distribution
       UN = RNORM()
       Z = UN*D + GAM*DEL
       COORDS(1,SC) = Z
    ENDDO

    ! sample the coupled coordinates
    IF (PRESENT(STARTCOORDS)) THEN
       CALL MCSAMPLERRHOPHI(DEL,EB,EPERPH,ETA,NSAMP,SAMPEVERY,INITSAMP,COORDS(2:4,:),FACC,FINALRANGES,STARTCOORDS)
    ELSE
       CALL MCSAMPLERRHOPHI(DEL,EB,EPERPH,ETA,NSAMP,SAMPEVERY,INITSAMP,COORDS(2:4,:),FACC,FINALRANGES)
    ENDIF

  END SUBROUTINE SAMPLERELCOORDSMC

  SUBROUTINE MCSAMPLERRHOPHI(DEL,EB,EPERPH,ETA,NSAMP,SAMPEVERY,INITSAMP,COORDS,FACC,FINALRANGES,STARTCOORDS)
    ! sample the coupled relative coordinates R, RHO, PHI
    ! using metropolis monte carlo
    ! should be better than the rejection sampling above for very stiff parameters
    ! get a total of NSAMP samples
    ! skip the first INITSAMP steps, then sample every SAMPEVERY steps
    ! COORDS(:,I) contains values of R,RHO,PHI for the Ith samp
    ! FACC gives the overall fraction accepted
    ! optionally, STARTCOORDS(1:3) gives the starting r,rho,phi coordinates
    ! and STARTCOORDS(4:6) gives the starting step ranges for r,rho,phi

    USE KEYS, ONLY : MCPRINTFREQ,ADJUSTEVERY,FACCTARGET,FACCTOL,ADJUSTSCL
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: DEL,EB,EPERPH,ETA
    INTEGER, INTENT(IN) :: NSAMP,SAMPEVERY,INITSAMP
    DOUBLE PRECISION, INTENT(OUT) :: COORDS(3,NSAMP), FACC, FINALRANGES(3)
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: STARTCOORDS(6)
    DOUBLE PRECISION :: R, RHO, PHI,ENERGY, EPREV, RPREV,RHOPREV,PHIPREV
    DOUBLE PRECISION :: DRMIN,DRHOMIN,DRHOMAX,RSTEP,RHOSTEP,PHISTEP
    DOUBLE PRECISION :: EPD,EBD,ETAD, DELE, TMP, U, U2, DRHO
    INTEGER :: STEP, TOTSTEPS, SAMPCT, TOTACCEPT
    LOGICAL :: ACCEPT, ADJUSTED

    TOTACCEPT= 0
    SAMPCT = 0      

    EPD = EPERPH/2/DEL
    EBD = EB/DEL
    ETAD = ETA*EB/DEL

   ! print*, 'testx1:',EB, ETA,  ebd, etad, EPD

    IF (PRESENT(STARTCOORDS)) THEN
       R = STARTCOORDS(1); RHO = STARTCOORDS(2); PHI = STARTCOORDS(3);
       RSTEP = STARTCOORDS(4); RHOSTEP = STARTCOORDS(5); PHISTEP = STARTCOORDS(6)
    ELSE
       PHI = 0D0; RHO = 1D0-min(sqrt(1/ebd)*0.1,1d0); R = SQRT(1D0/EPD)*0.1
    ENDIF
    !EPREV = EPD*R**2
    EPREV = -LOG(R)+EPD*R**2 - ETAD*R*COS(PHI)*SQRT(1-RHO**2)  - EBD*RHO
    !print*, 'testx0:', eprev, r,rho,phi

    RSTEP = SQRT(1D0/EPD);    
    RHOSTEP = MIN(SQRT(1D0/EBD),1D0);
    PHISTEP= MIN(SQRT(1D0/ETAD),2*PI);        

    TOTSTEPS = INITSAMP+NSAMP*SAMPEVERY
    DO STEP = 1,TOTSTEPS
       RPREV = R; RHOPREV = RHO; PHIPREV = PHI;

       ! take a step, keeping R, RHO, and PHI within bounds
       !DRMIN = MIN(R,RSTEP)
       U = GRND()
       !R = R + U*(RSTEP+DRMIN)-DRMIN
       R = R + U*2*RSTEP - RSTEP

       PHI = PHI + GRND()*2*PHISTEP - PHISTEP

       IF (PHI.GT.2*PI) THEN
          PHI = PHI - 2*PI
       ELSEIF (PHI.LT.0) THEN
          PHI = PHI + 2*PI
       ENDIF

       ! DRHO = MIN(MIN(1-RHO,RHO+1),RHOSTEP)         
       ! U2 = GRND()
       ! PRINT*, 'TESTX1:', RHO, DRHO, RHOSTEP, U2*2*DRHO-DRHO
       ! RHO = RHO + U2*2*DRHO-DRHO
       ! IF (ABS(RHO).GT.1) THEN
       !    PRINT*, 'ERROR IN MCRELSAMPLE: BAD RHO', RHO
       !    STOP 1
       ! ENDIF
       !DRHOMIN = MIN(RHO+1,RHOSTEP)
       !DRHOMAX = MIN(1-RHO,RHOSTEP)
       U2 = GRND()       
       !PRINT*, 'TESTX2:', RHO, drhomin,drhomax,U2*(DRHOMIN+DRHOMAX)-DRHOMIN
       !RHO = RHO + U2*(DRHOMIN+DRHOMAX)-DRHOMIN
       RHO = RHO + U2*2*RHOSTEP - RHOSTEP
       !IF (RHO.LT.-1D0) THEN
       !   RHO = -1D0+EPSILON(1D0)
       !ELSEIF(RHO.GT.1D0) THEN
       !   RHO = 1D0+EPSILON(1D0)
       !ENDIF       

       IF (ABS(RHO).GT.1.OR.R.LT.0) THEN
         ACCEPT = .FALSE.
       ELSE
         ACCEPT = .TRUE.
       ENDIF

       IF (ACCEPT) THEN
          ! calculate change in energy
          !ENERGY = EPD*R**2
          ENERGY = -LOG(R)+EPD*R**2 - ETAD*R*COS(PHI)*SQRT(1-RHO**2)  - EBD*RHO

          DELE = ENERGY - EPREV

          ! decide whether to accept
          IF (DELE.LT.0) THEN
             ACCEPT = .TRUE.
          ELSE
             TMP = GRND()
             ACCEPT = (TMP.LT.EXP(-DELE))
          ENDIF
       ENDIF

       IF (ACCEPT) THEN
          EPREV = ENERGY

          TOTACCEPT = TOTACCEPT+1
       ELSE
          ! restore old coords
          R = RPREV; RHO = RHOPREV; PHI = PHIPREV
          ENERGY = EPREV
       ENDIF

       FACC = TOTACCEPT/DBLE(STEP)

       IF (MOD(STEP,MCPRINTFREQ).EQ.0) THEN
          PRINT '(A,I20,7G20.10)', 'MCSTEP:', STEP, ENERGY, FACC, R, RHO, PHI, U, RHOSTEP
       ENDIF

       IF (STEP.GT.INITSAMP.AND.MOD(STEP-INITSAMP,SAMPEVERY).EQ.0) THEN
          SAMPCT = SAMPCT+1
          COORDS(:,SAMPCT) = (/R,RHO,PHI/)
       ENDIF

       IF (ADJUSTEVERY.GT.0.AND.MOD(STEP,ADJUSTEVERY).EQ.0) THEN
          ! check whether ranges need adjusting
          IF (FACC.LT.FACCTARGET-FACCTOL) THEN
             ! PRINT*, 'TESTX1:', RSTEP, RSTEP/ADJUSTSCL
             RSTEP = RSTEP/ADJUSTSCL
             RHOSTEP = RHOSTEP/ADJUSTSCL
             PHISTEP = PHISTEP/ADJUSTSCL
             ADJUSTED = .TRUE.
          ELSEIF (FACC.GT.FACCTARGET+FACCTOL) THEN
             RSTEP = RSTEP*ADJUSTSCL
             RHOSTEP = MIN(RHOSTEP*ADJUSTSCL,2D0)
             PHISTEP = PHISTEP*ADJUSTSCL
             ADJUSTED = .TRUE.
          ELSE
             ADJUSTED = .FALSE.
          ENDIF

          !IF (ADJUSTED) PRINT*, 'ADJUSTED RANGES:', RSTEP, RHOSTEP, PHISTEP
       ENDIF
    ENDDO

    FINALRANGES = (/RSTEP,RHOSTEP,PHISTEP/)
  END SUBROUTINE MCSAMPLERRHOPHI

  SUBROUTINE SAMPLERRHOMVN(DEL,EB,EPERPH,ETA,NSAMP,RVALS,RHOVALS,AVGTRY)
    ! rejection sampling of R and RHO coordinates using multivariate normal
    ! much more efficient for high values of eta than the old version of SAMPLERRHOMVN
    ! see notes from 4/30/2013
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: DEL,EB,EPERPH,ETA
    INTEGER, INTENT(IN) :: NSAMP
    DOUBLE PRECISION, INTENT(OUT) :: RVALS(NSAMP),RHOVALS(NSAMP), AVGTRY
    INTEGER, PARAMETER :: MAXNTRY = 1000
    DOUBLE PRECISION :: A,B,C,D, C2, E
    DOUBLE PRECISION :: R,V,U, COV(2,2),MU(2), DISC, PAIR(2), RATIO, I0VAL
    INTEGER :: TRY, SC, INFO, TOTTRY
    

    ! Set up the multivariate normal envelope distribution for sampling R and V
    ! V = sqrt(1-rho)

    A = EPERPH/2/DEL
    B = EB/DEL
    C = ETA*EB/DEL
    C2 = SQRT(2D0)*C
    D = SQRT(2*A)
    E = SQRT(2*B)

    !PRINT*, 'TESTX0:', A,B,C,D

    ! get mean vector
    DISC = 4*A*B-C2*C2
    IF (DISC.LE.0) THEN
       PRINT*, 'ERROR IN SAMPLERELCOORDSMVN: negative discriminant. Not a well defined distribution.', A, B, C, DISC
       STOP 1
    END IF
    MU = (/2*B*D+C2*E,C2*D+2*A*E/)/DISC
    ! get covariance matrix (only lower triangular half is nonzero)
    COV(1,:) = (/2*B,0D0/)
    COV(2,:) = (/C2,2*A/)
    COV = COV/DISC
    
    ! get cholesky decomposition of covariance matrix
    CALL DPOTRF('L',2,COV,2,INFO)
    IF (INFO.NE.0) THEN
       PRINT*, 'ERROR IN SAMPLERELCOORDSMVN: Cholesky decomposition of covariance matrix failed:', INFO
       STOP 1
    ENDIF

    TOTTRY = 0
    DO SC = 1,NSAMP
       DO TRY = 1,MAXNTRY
          ! sample a multivariate normal pair
          PAIR = MVNORM(2,MU,COV)
          R=PAIR(1); V = PAIR(2)
          
          !print*, 'testx2:', sc, try, r, v
          IF (V.LE.0.OR.V.GE.2.OR.R.LE.0) THEN
             ! reject if R or V are out of bounds
             CYCLE
          ENDIF
          
          ! uniform deviate
          U = GRND()
          ! accept if below ratio
          CALL BESSF_I0(C*R*V*SQRT(2-V*V),I0VAL)
          RATIO = R*V*I0VAL*D*E*EXP(-C2*R*V - D*R - E*V + 2)
          !print*, 'TESTX1:', sc, TRY, RATIO, u
          IF (U.LT.RATIO) THEN             
             EXIT
          ENDIF
       END DO     
       TOTTRY = TOTTRY + TRY
       IF (TRY.GE.MAXNTRY) THEN
          PRINT*, 'ERROR IN SAMPLERRHOMVN: failed to generate acceptable sample'
          STOP 1
       ENDIF

       RVALS(SC) = R
       RHOVALS(SC) = 1 - V*V
    ENDDO
    
    AVGTRY = TOTTRY/NSAMP
  END SUBROUTINE SAMPLERRHOMVN

  SUBROUTINE SAMPLEPHICOND(NSAMP,COEFF,RVALS,RHOVALS,PHIVALS,AVGTRY)
    ! sample PHI from a conditional distribution with given R and RHO values
    ! use rejection sampling with normal distrib envelope

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: NSAMP
    DOUBLE PRECISION, INTENT(IN) :: RVALS(NSAMP), RHOVALS(NSAMP), COEFF
    DOUBLE PRECISION, INTENT(OUT) :: PHIVALS(NSAMP), AVGTRY
    INTEGER, PARAMETER :: MAXNTRY = 1000
    INTEGER :: SC, TRY, tottry
    DOUBLE PRECISION :: D, PHI, U, RATIO

    TOTTRY = 0
    DO SC = 1,NSAMP
       D = COEFF*RVALS(SC)*SQRT(1 - RHOVALS(SC)**2)

       DO TRY = 1,MAXNTRY
          ! sample from normal distrib
          PHI = RNORM()/SQRT(0.4*D)

          IF (ABS(PHI).GT.PI) CYCLE

          ! sample from uniform
          U = GRND()
          ! decide whether to reject
          RATIO = EXP(D*(COS(PHI)-1+0.2*PHI*PHI))

!          print*, 'testx1:', phi, u, ratio, try
          IF (U.LT.RATIO) THEN
             EXIT
          ENDIF
       ENDDO
       IF (TRY.GE.MAXNTRY) THEN
          PRINT*, 'ERROR IN SAMPLEPHICOND: failed to generate acceptable sample.'
          STOP 1
       ENDIF
       PHIVALS(SC) = PHI
       TOTTRY = TOTTRY + TRY
    ENDDO
    
    AVGTRY = dble(TOTTRY)/NSAMP
  END SUBROUTINE SAMPLEPHICOND

  SUBROUTINE SAMPLERELCOORDSMVN(DEL,EB,GAM,EPAR,EPERPH,ETA,NSAMP,COORDS)
    ! sample relative coordinates Z, R, RHO, PHI for segment junctions
    ! generates NSAMP samples and stores them in the NSAMPx4 array COORDS
    ! uses rejection sampling with multivariate normal (see notes from 4/30/2013)
    ! Should be significantly more efficient for high eta values
     ! EPERPH corresponds to eperp_hat in the notes, but EPERP in the simulation code
    ! ETA is -EC/LP in the simulation code and EB=LP
    

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: DEL, EB, GAM,EPAR,EPERPH, ETA
    INTEGER, INTENT(IN) :: NSAMP
    DOUBLE PRECISION, INTENT(OUT) :: COORDS(4,NSAMP)
    DOUBLE PRECISION :: ZSIG, U, AVGTRY
    INTEGER :: SC

    ! sample the Z values from a normal distribution
    ZSIG = SQRT(DEL/EPAR)
    DO SC = 1,NSAMP
       U = RNORM()
       COORDS(1,SC) = U*ZSIG + GAM*DEL
    ENDDO

    ! sample the R and RHO values using rejection sampling with multivariate normal
    CALL SAMPLERRHOMVN(DEL,EB,EPERPH,ETA,NSAMP,COORDS(2,:),COORDS(3,:),AVGTRY)

    ! sample the PHI values conditional on the r, rho
    ! using rejection sample with normal envelope
    CALL SAMPLEPHICOND(NSAMP,ETA*EB/DEL,COORDS(2,:),COORDS(3,:),COORDS(4,:),AVGTRY)
  END SUBROUTINE SAMPLERELCOORDSMVN
END MODULE SAMPLEUTIL
