C! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calculation  of the  total cross-section integrad qq-> hh
c      In the UEHiggsY   model or SM if the Wilson coefficients are set to zero
C   written by Lina Alasfar, HU Berlin  04.04.2019                          c
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C-- qq > hh
      SUBROUTINE UPEVNT
C...All real arithmetic in double precision.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C...User process event common block.
      INTEGER MAXNUP, FLAV,RF, FLAV2, mod
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      DOUBLE PRECISION  MH,MW,MZ,lambdaEW,lambdaNP, gammah, m
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,IDUP(MAXNUP),
     &ISTUP(MAXNUP),MOTHUP(2,MAXNUP),ICOLUP(2,MAXNUP),PUP(5,MAXNUP),
     &VTIMUP(MAXNUP),SPINUP(MAXNUP)
      SAVE /HEPEUP/

C...The user's own transfer of information.
      COMMON/MYCOMM/ECM,PTMIN,INIT,mod
      SAVE/MYCOMM/
        COMMON/KIN/sh,uh,th,MH,M,MQ
        COMMON/PARAM/WILSON, WC, ghff,ghhff,gammah,ghhh, lambdanp,lambdaEW,xkl,xkf
C      common/ flavour/rf,fp1,fp2
C...Local arrays and parameters.
      DIMENSION  wc(6,6),fp1(-6:6),fp2(-6:6),xkf(6), xs(5,5)
       DOUBLE PRECISION, DIMENSION(6) :: MQ
      DATA PI/3.141592653589793D0/
      DATA CONV/0.3894D9/
C----------------------------first: set up the parameters------------------------------------
C...Default flag for event rejection.
      IREJ=0
C--  initilise the kappa_f values w.r.t the beauty coupling
C          do 556, i = 1,4
C            xkf(i) = pmas(5,1)/pmas(i,1)
C556        continue
C           xkf(2)= 2000d00





C--    Other parameters
C          xkl = 1d0                            !   NP contibution to the trilinear cpupling
          gammah = PMAS(25,2)                  !     Higgs width
          gf = 1.166378D-05                    ! Fermi's constant [GeV^-2]
         lambdaEW = (dSQRT(2D00)* GF)**(-0.5) ! the vev or electroweak scale [GeV]
         lambdaNP = 1300D00                   ! new physics scale ~ TeV scale  [GeV]
         MW =PMAS(24,1)                       ! W mass  [GeV]
         mh=PMAS(25,1)                        ! Higgs mass
         g=sqrt(8*mw**2*gf/sqrt(2d0))         ! weak coupling cnstant
         ghhh=(xkl)*(6d0*g*mh**2)/(4d0*mw)         ! Higgs trilinear coupling
         hwidth = cmplx(0d0,-gammah*mh)
C...Default return value.
      XWGTUP=0D0         ! the event weight
      dsiggma =0d0       ! storing the total matrix element
       do 657, j = 1,5,1
         do 656 ,i =1,5,1
           xs(i,j) = 0d0    ! flavour matrix element
656      continue
657      continue
C--  Setting up the Wilson coeffecients
       DO 15, ncol = 1,6
         DO 16, nrow = 1,6
            wc(ncol,nrow) = 0.0D00    !  off diagonal wilson , FCNC
           if( (ncol.eq.nrow).and.(ncol.lt.5) )then
                 ! first one k is in terms of beauty couplings ratio
              wc(ncol,nrow) =
     .        (lambdaNP/lambdaEW)**2 *
     .        (pmas(4,1)/lambdaEW) * (xkf(nrow)-pmas(nrow,1)/pmas(4,1))
                  ! this one , k is in terms of ratio of yukawa to the SM
C               wc(ncol,nrow) =
C     .              (lambdaNP/lambdaEW)**2 *
C     .              (pmas(nrow,1)/lambdaEW) * (XKF(nrow)-1)
                endif
16           CONTINUE
15         CONTINUE
C---------------------------Second: Phase space and mandalstam variables--------------------------------
      SHATMIN = (2D0*MH)**2! +PTMIN**2
      SHATMAX = ECM**2
      TAUM = SHATMIN/ECM**2
      TAU=  TAUM**PYR(0) !(2D0*HM)**2/ECM**2 !4D0*PTMIN**2/ECM**2
       X2 =   TAU**PYR(0)
       X1= TAU/X2
       s=tau*ECM**2
       sh = s

      T2=1.d0/2.D0*(-s+dSqrt(s)*2d0*dSQRT(s/4d0-mh**2))+mh**2
      T1=1.d0/2.D0*(-s-dSqrt(s)*2d0*dSQRT(s/4d0-mh**2))+mh**2

      T=T1+(T2-T1)*PYR(0)


C...Derive other kinematical quantities.
      U=2D00*MH**2-S-T
      PT2=(T*U-S*mh**2)/S
C...Pick Q2 scale (which involves some arbitrariness)
      Q2= 0.25d0*S!PT2
      SCALUP=SQRT(Q2)
C--   EVOLVE PDF'S
      call evolvePDF(x1,SCALUP,fp1)
      call evolvePDF(x2,SCALUP,fp2)

C--   the jacobian
      PHSPV=(-dLOG(TAUM))*(-dLOG(TAU))*(T2-T1) !dlog(T2/T1)*that

C------------------------Third: compute the matrix elements ----------------------------------
C SELECT FLAVOUR AT RANDOM
      RF = FLOOR( (5-1)*PYR(0)+1 ) ! just for generation asthetics
C--   The flavour diagonal matrix elements
      do 666, i = 1,5,1
80     ghff =  (pmas(i,1)/lambdaEW)!*xkf(i)
C       ghhff=m/lambdaEW**2*(3d0/2d0*(1d0-mq(5)/(m)))
C       ghhff = 3*(pmas(i,1)/lambdaEW**2)*(xkf(i)-1)/2d0
       ghhff = (pmas(i,1)/lambdaEW**2)*xkf(i)


C-- t and u channels with their interference
        xs(i,i) =(ghff)**4*((2*mh**2*s)/t**2 - (2*s*(mh**2 - t))/t**2
     .    + (2*mh**2*s)/(t*u**2)
     .    - (2*s*(mh**2 - u))/(t*u**2)
     .    + (s*(-2*mh**2 + s))/(t**2*u)
     .    - (mh**2 - t)**2/(t**2*u) + (s*(-2*mh**2 + s))/(t*u)
     .    - (2*s*(mh**2 - t))/(t*u)
     .    + (mh**2 - t)**2/(t*u) - (2*s*(mh**2 - u))/(t**2*u)
     .     + (mh**2 - u)**2/(t**2*u)
     .    - (mh**2 - u)**2/(t*u))
C-- s , with NP channels with their interference
      xs(i,i) = xs(i,i)+(2d0)*s*((ghff*ghhh / (-mh**2 + s)) +2d0* ghhff)**2
      xs(i,i) =  xs(i,i)/s**2
C-- weight by the parton distibution functions
      xs(i,i) = xs(i,i)*fp1(i)*fp2(-i)
     .                +xs(i,i)*fp2(i)*fp1(-i)
666     CONTINUE
C  FCNC elements :
      do 777 i = 1,5,1
         do 778 j=1,5,1
            if(i.ne.j) then  ! only off diagonal
            ghff=(lambdaEW / lambdaNP)**2 * wc(i,j) ! only NP part of Yukawa
            ghhff = -3d0* (lambdaEW / lambdaNP**2) * wc(i,j)
          ! t and u channels are not allowed with FCNC
             xs(i,j) = 2*s*((ghff*ghhh / (-mh**2 + s))+ ghhff)**2
C             xs(i,j) =  xs(i,j)/s**2
             xs(i,j) = 0d0!xs(i,j)*fp1(i)*fp2(-j)
C     .                +xs(i,j)*fp2(j)*fp1(-i)
            end if
778       continue
777       continue


C...The CONV factor converts from GeV**(-2) to pb.
      factor1 =1d0/(16d0*PI) *(1d0/24d0)  ! prefactor
      dsiggma = abs(sum(xs))   ! the total matrix elemets
      XWGTUP= CONV*PHSPV*dsiggma*factor1 !the differential cross section
C...Done when initializing for maximum cross section.
      IF(INIT.EQ.1) RETURN
C...Store scale choice etc.


C---------------------Fourth: setting up the event------------------------------

C...The third step is to set up the partonic process that is selected.

C...Define number of partons - two incoming and two outgoing.
      NUP=4

C...Flavour codes for entries. Note that definition of t-hat
C...means quark either is 1 and 3 or 2 and 4.
      IDUP(1)=1
      IDUP(2)=-1
      IDUP(3)=25
      IDUP(4)=35

C...Status codes.
      ISTUP(1)=-1
      ISTUP(2)=-1
      ISTUP(3)=1
      ISTUP(4)=1

C...Mother codes.
      MOTHUP(1,1)=0
      MOTHUP(2,1)=0
      MOTHUP(1,2)=0
      MOTHUP(2,2)=0
      MOTHUP(1,3)=1
      MOTHUP(2,3)=2
      MOTHUP(1,4)=1
      MOTHUP(2,4)=2

C...Colour flow.
C...Colour stretched from initial quark to  antiquark
        ICOLUP(1,1)  =501
        ICOLUP(2,1)=0
        ICOLUP(1,2)=0
        ICOLUP(2,2)=501
        ICOLUP(1,3)  =0
        ICOLUP(2,3)=0
        ICOLUP(1,4)=0
        ICOLUP(2,4)=0


C...Reset momenta to zero.
      DO 130 I=1,4
        DO 120 J=1,5
          PUP(J,I)=0D0
120   CONTINUE
130   CONTINUE
CC
C...Masses of final state entries; initial assumed massless.
      PUP(5,3)= PYMASS(IDUP(3))
      PUP(5,4)=PYMASS(IDUP(4))
CC        sqrt2 = dsqrt(2d00)
C...Four-momenta of the incoming partons simple.
      PUP(4,1)=X1*ECM*0.50d00
      PUP(3,1)=PUP(4,1)
      PUP(4,2)=ECM*X2*0.50d00
      PUP(3,2)=-PUP(4,2)

C...Energies and absolute momentum of the outgoing partons in
C...the subsystem frame.
      RTSHAT=SQRT(X1*X2)*ECM
      PABS=0.5D0*SQRT(MAX(0D0,(RTSHAT**2-PUP(5,3)**2-
     &PUP(5,4)**2)**2-4D0*PUP(5,3)**2*PUP(5,4)**2))/RTSHAT
      PE3=0.5D0*(RTSHAT**2+PUP(5,3)**2-PUP(5,4)**2)/RTSHAT
      PE4=RTSHAT-PE3
C...Subsystem scattering angle defined
      COSTHE=(t-u)/s
      SINTHE=SQRT(MAX(0D0,1D0-COSTHE**2))
C...Azimuthal angle at random.
      PHI=2D0*PI*PYR(0)
C...Momenta of outgoing partons in the subsystem frame.
      PUP(1,3)=PABS*SINTHE*COS(PHI)
      PUP(2,3)=PABS*SINTHE*SIN(PHI)
      PZ3=PABS*COSTHE
      PUP(1,4)=-PUP(1,3)
      PUP(2,4)=-PUP(2,3)
      PZ4=-PZ3
C...Longitudinal boost of outgoing partons to cm frame.
      BETA=(X1-X2)/(X1+X2)
      GAMMA=0.5D0*(X1+X2)/SQRT(X1*X2)
      PUP(3,3)=GAMMA*(PZ3+BETA*PE3)
      PUP(4,3)=GAMMA*(PE3+BETA*PZ3)
      PUP(3,4)=GAMMA*(PZ4+BETA*PE4)
      PUP(4,4)=GAMMA*(PE4+BETA*PZ4)
      RETURN
      END
CC


********************************************************************
