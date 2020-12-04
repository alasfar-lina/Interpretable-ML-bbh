C! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calculation  of the  total cross-section integrad qq-> hh
c      In the UEHiggsY   model or SM if the Wilson coefficients are set to zero
C   written by Lina Alasfar, HU Berlin  04.04.2019                          c
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C --- gg > hh
      SUBROUTINE UPEVNT
       use COLLIER
C...All real arithmetic in double precision.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C-- comple stuff
      double complex :: C0,C01,C02,C03,C04,C05
      double complex :: D0,D01,D02,D03,D04
      double complex :: cm2,cmh2,cs,ct,cu,cz
      complex*16 Amp1,Amp2,Amp3,Ampli3z,Ampli2z,Ampli1z,
     -A1,a2,a3,cghhff
C...User process event common block.
      INTEGER MAXNUP, FLAV,RF, fi
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      DOUBLE PRECISION :: MH,MW,MZ,lambdaEW,lambdaNP, gammah, m
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,IDUP(MAXNUP),
     &ISTUP(MAXNUP),MOTHUP(2,MAXNUP),ICOLUP(2,MAXNUP),PUP(5,MAXNUP),
     &VTIMUP(MAXNUP),SPINUP(MAXNUP)
      SAVE /HEPEUP/
          DOUBLE PRECISION, DIMENSION(6) :: MQ
C...The user's own transfer of information.
      COMMON/MYCOMM/ECM,PTMIN,INIT
      SAVE/MYCOMM/
      COMMON/KIN/SHAT,UHAT,THAT,MH,M,MQ,pt2
      COMMON/PARAM/WILSON, WC, ghff,ghhff,gammah,ghhh,lambdanp,lambdaEW,xkl
       COMMON/ALS/XLAMBDA,AMC,AMB,AMT,N0
C...Local arrays and parameters.
       DIMENSION  wc(6,6),fp1(-6:6),fp2(-6:6),xkf(6), xs(5,5),xf2(6)
      DATA PI/3.141592653589793D0/
      DATA CONV/0.3894D9/

      !   in loop quark masses , use these for numerical stability
             MQ(6) = 173.03D0              ! top mass (on shell)  [GeV]
             MQ(5) = 4.18D0               ! beauty mass   [GeV]
             MQ(4) = 0d0!1.275D0              ! charm mass    [GeV]
             MQ(3) =0d0! 97.00D-03                      ! strange mass [GeV]
             MQ(2) = 0d0!4.725D-03                     ! down mass [GeV]
             MQ(1) =  0d0!2.225D-03                      ! up mass [GeV]

C-------------------Zeroth : initlise the alphas_running--------------------------------------
        XLAMBDA = 0.08782725D00       ! select such that alpha_s at MZ is 0.118
        N0 = 5
        AMT= 173.03D0
        AMB = 4.18D0
        AMC =1.275D0
        ACC = 1.0D-08
       CALL ALSINI(ACC)
C        COLLIER stuff
        call Init_cll(2,2,"output_cll", .TRUE.) ! initilise the Collier subroutines
         call InitCacheSystem_cll(1,2)
         call InitEvent_cll
C----------------------------first: set up the parameters------------------------------------
C...Default flag for event rejection.
      IREJ=0
C--  initilise the kappa_f values w.r.t the beauty coupling
            xkf(5) = 1d0
            xkf(6) = 1d0
          do 555, i = 1,4
            xkf(i) = 1d0! pmas(5,1) /pmas(i,1)
555        continue
        xmH = 600.0d0
        gammaH = 3.901509d0

C  non-linear ghhff couplings
         xf2(1) = 0.046402d0
         xf2(2) = 0.021720d0
         xf2(3) = 0.334588d0
         xf2(4) = 5.524813d0
         xf2(5) = 10.284945d0
         xf2(6) = 18.397843d0
C--    Other parameters
C         cpupling
          gammah = PMAS(25,2)                  !     Higgs width
          gf = 1.166378D-05                    ! Fermi's constant [GeV^-2]
         lambdaEW = (dSQRT(2D00)* GF)**(-0.5) ! the vev or electroweak scale [GeV]
         lambdaNP = 1300D00                   ! new physics scale ~ TeV scale  [GeV]
         MW =PMAS(24,1)                       ! W mass  [GeV]
         mh=PMAS(25,1)                        ! Higgs mass
         g=sqrt(8*mw**2*gf/sqrt(2d0))         ! weak coupling cnstant
         ghhh=(xkl)*(6d0*g*mh**2)/(4d0*mw)         ! Higgs trilinear coupling
         hwidth = cmplx(0d0,-gammah*mh)

!
C...Default return value.
      XWGTUP=0D0         ! the event weight
      dsiggma =0d0       ! storing the total matrix element

C---------------------------Second: Phase space and mandalstam variables--------------------------------
      SHATMIN = (2D0*MH)**2! +PTMIN**2
      SHATMAX = ECM**2
      TAUM = SHATMIN/ECM**2
      TAU=  TAUM**PYR(0) !(2D0*HM)**2/ECM**2 !4D0*PTMIN**2/ECM**2
       X2 =   TAU**PYR(0)
       X1= TAU/X2
       s=tau*ECM**2

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

CC
C---------------------------Third: Matrix elements -------------------------------------------------
      alps=alphas(SCALUP,1)  ! alpha s running ( not by Pythia )
C
       Ampli3z=(0d0,0d0)
      Ampli2z=(0d0,0d0)
      Ampli1z=(0d0,0d0)
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        m = mq(6)
        fi = 6
80     ghff =  (m/lambdaEW)*xkf(fi)

C      cghhff = xf2(fi)/cmplx(s-xmH**2,-xmH*gammaH)
      ghhff = 3*(m/lambdaEW**2)*(xkf(fi)-1)
      ! complixified quantities for the loop integrals
        cz = 0d0
         if(m.eq.mq(6)) then
         cm2 = m**2
          else
             cm2 = mq(fi)**2
           end if
        cmh2 = mh**2
        cs = s
        ct = t
        cu = u

      ! three point functions
         call C0_cll(C0,cz,cs,cz,cm2,cm2,cm2) ! trinagle and  3pt function s channel
         call C0_cll(C01,cmh2,cu,cz,cm2,cm2,cm2) ! box 3pt function  u channel
         call C0_cll(C02,cmh2,ct,cz,cm2,cm2,cm2) ! box 3pt function t channel
         call C0_cll(C03,cz,cmh2,cu,cm2,cm2,cm2) ! box 3pt function u channel
         call C0_cll(C04,cmh2,cs,cmh2,cm2,cm2,cm2) ! box 3pt function s channel ii
         call C0_cll(C05,cz,cmh2,ct,cm2,cm2,cm2) ! box 3pt function t channel ii

      ! four point functions
        call D0_cll(D0,cz,cz,cmh2,cmh2,cs,ct,cm2,cm2,cm2,cm2) ! st channel
        call D0_cll(D01,cz,cz,cmh2,cmh2,cs,cu,cm2,cm2,cm2,cm2) ! su channel
        call D0_cll(D02,cz,cz,cmh2,cmh2,ct,cu,cm2,cm2,cm2,cm2) ! tu channel

!   box diagram, spin 0
         Amp1=
     -        (1d0/s*ghff**2*
     - (4*s+8*s*m**2*C0-2*s*
     - (s+2*mh**2-8*m**2)*m**2*(D0+D01+D02)+
     - (2*mh**2-8*m**2)*(2*(t-mh**2)*C02+2*(u-mh**2)
     - *C01-(t*u-mh**4)*D02)) )

! box diagram, spin 2
       Amp2=
     -      (ghff**2*
     - 1d0/(t*u-mh**4)*((t**2+mh**4-8*t*m**2)*
     -(s*C0+2*(t-mh**2)*C05
     - -s*t*D0)+(u**2+mh**4-8*u*m**2)
     - *(s*C0+2*(u-mh**2)*C03-s*u*D01)-(t**2+u**2-2*mh**4)*
     - (t+u-8*m**2)*C04-2*(t+u-8*m**2)*(t*u-
     - mh**4)*m**2*(D0+D01+D02)) )
! Triangle diagram
      Amp3=
     -      (2d0*ghff*ghhh*
     -   m*(2 + (4*m**2 - s)*C0)) / (-mh**2 + s+hwidth) +
     -   2d0*ghhff*
     -   m*(2 + (4*m**2 - s)*C0) ! effective hh ff coupling dim 6 operator
c--- Amplitude together.... triangel has got same projector as factor1 of box
         Ampli3z=Ampli3z+Amp3
         Ampli2z=Ampli2z+Amp2
         Ampli1z=Ampli1z+Amp1
C  Include all the quark flavours
      DO 666, FLAV = 0,4,1
          IF(M.EQ.pmas(6-FLAV,1)) THEN
            M=pmas(5-FLAV,1)
          fi = fi -1
            GOTO 80
       ENDIF
666     CONTINUE
      A1=Ampli1z
      A2=Ampli2z
      A3=Ampli3z
      Ampli=(A1)*Conjg(A1)+ (A2)*Conjg(A2)+
     . (A3)*Conjg(A3)+A3*Conjg(A1)+A1*Conjg(A3)

        AQEDUP=AEM
        AQCDUP=ALPS
        factor=-1.d0/(s**2)*1d0/1024d0*1d0/pi**3*ALPS**2/8d0
        XWGTUP=-CONV*PHSPV*Ampli*factor*fp1(0)*fp2(0)
C---------------------------Fourth: Setting up the event  -------------------------------------------------
CC
C...The third step is to set up the partonic process that is selected.
CC
C...Define number of partons - two incoming and two outgoing.
      NUP=4
CC
C...Flavour codes for entries. Note that definition of t-hat
C...means quark either is 1 and 3 or 2 and 4.
      IDUP(1)=21
      IDUP(2)=21
      IDUP(3)=25
      IDUP(4)=35
CC
C...Status codes.
      ISTUP(1)=-1
      ISTUP(2)=-1
      ISTUP(3)=1
      ISTUP(4)=1
CC
C...Mother codes.
      MOTHUP(1,1)=0
      MOTHUP(2,1)=0
      MOTHUP(1,2)=0
      MOTHUP(2,2)=0
      MOTHUP(1,3)=1
      MOTHUP(2,3)=2
      MOTHUP(1,4)=1
      MOTHUP(2,4)=2
CC
C...Colour flow.
C...Colour stretched from initial quark to  antiquark
        ICOLUP(1,1)  =501
        ICOLUP(2,1)=502
        ICOLUP(1,2)=502
        ICOLUP(2,2)=501
        ICOLUP(1,3)  =0
        ICOLUP(2,3)=0
        ICOLUP(1,4)=0
        ICOLUP(2,4)=0
CC
CC
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
