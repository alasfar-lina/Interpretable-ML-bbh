ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                               HH production                          c
c               With PYTHIA 6 Implementation , SM and EFT BSM          c
c                       WRITTEN BY : LINA ALASFAR                      c
c                          HU-BERLIN   25.03.2019                      C
C                        CURRENT REVISION :: 09.04.2019                c
c  REQUIRMENTS :  LHAPDF, COLLIER FOR SCALAR LOOP INTEGRALS            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                       VERSION 0.1                                   c
c                      Version 0.2 :: added progress bar,             c
c                       Verson  0.3  ::  added qqbar with FCNC        c
c                       Version 1.0  :: Interface with ROOT
c                       Version 1.1 add acceptence      pre cuts   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C-----------------------------------------------------------------

C...Preamble: declarations.

C...All real arithmetic in double precision.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C...Three Pythia functions return integers, so need declaring.
      INTEGER PYK,PYCHGE,PYCOMP,L
        DOUBLE PRECISION :: MH,MW,MZ,lambdaEW,lambdaNP, gammah,m
        DIMENSION  wc(6,6),fp1(-6:6),fp2(-6:6),xkf(6), xs(5),XGG(4),XBB(4), W(4), bins(7)

       DOUBLE PRECISION, DIMENSION(6) :: MQ
       CHARACTER(LEN = 50) :: file
       CHARACTER(LEN=4) :: BORIG
       CHARACTER(LEN=1) ::BX,BY,BZ,BT
       CHARACTER(LEN=1) :: BMASS
       CHARACTER(LEN=2) :: BVX, BVY, BVZ

C...EXTERNAL statement links PYDATA on most machines.
      EXTERNAL PYDATA

C...Commonblocks.
      INTEGER MAXNUP,ii,jj
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,IDUP(MAXNUP),
     &ISTUP(MAXNUP),MOTHUP(2,MAXNUP),ICOLUP(2,MAXNUP),PUP(5,MAXNUP),
     &VTIMUP(MAXNUP),SPINUP(MAXNUP)
      SAVE /HEPEUP/
C...The event record.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
C...Parameters.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
C...Parameters.C...Particle properties + some flavour parameters.
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
C...Decay information.
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000)
     &  ,KFDP(8000,5)
       COMMON/PYINT1/MINT(400),VINT(400)
C...Selection of hard scattering subprocesses.

C      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200),ISUB
C...Parameters.
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
C...settings transfer of information.
       COMMON/MYCOMM/ECM,PTMIN,INIT,mod
       SAVE/MYCOMM/
        COMMON/KIN/sh,uh,th,MH,M,MQ
        COMMON/PARAM/WILSON, WC, ghff,ghhff,gammah,ghhh, lambdanp,lambdaEW,xkl,xkf

C...Local arrays and parameters.
C-----------------------------------------------------------------

C...First section: initialization.

C...Number of events and cm energy.
       nev = 10000
C      write(*,*) "Input the number of events"
C                      read(*,*)  NEV
      W=  (/13d3, 14d3,27d3,100d3 /)
C      bins = (/-15d0, -10d0,-5d0,-2d0,-1d0,0d0,1d0,2d0,3d0,5d0
C     .,10d0,15d0 /)

      do 555, i = 1,6
       xkf(i) = 1d0
555        continue
        xkl =1d0
      bins= (/2.5d1, 5.0d1,1d2,1.5d2,1.75d2,2.5d2,2d2 /)
      do ii = 1, 4, 1
        do jj = 1, 7, 1

      ECM= W(ii)!14d03       ! COM energy
      xkf(3)= bins(jj)!15d0                 !   NP contibution to the trilinear
C...Initialize ROOT .
        write(file,'(a,i0,a,i0,a)') './root_files/strange/css',INT(bins(jj)),'-',INT(W(ii)/1000),'TeV.root'
C
C        print * , file
       call  ROOTINIT(file)
C       call  ROOTINIT('dum.root')
       call ROOTWRITE
      CALL PYINIT('USER',' ',' ',0D0)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C..  ::  MASTER PYTHIA  AND PHYSICS SENTTINGS ::                       C
C             ALL MASSES ARE IN GEV                                    C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
CC...Switch off unnecessary aspects: initial and  final state
C...showers, multiple interactions, hadronization.
C...(Optional for faster simulation of the parton-level
C...processes only.)
C...Higgs mass, and width
        PMAS(25,1)= 125d00 ! mass
        PMAS(25,2)= 0.3d0  ! width
        PMAS(35,1)= 125d00
        PMAS(35,2)= 0.3d0
C.. Higgs Decay channels (switched off)
          do 464 , i = 211,226,1
            MDME(i,1)=0
464       continue
C other Higgs channels off
          do 474 , i = 335,357,1
            MDME(i,1)=0
474       continue
C Higgs branching ratios, made equal for H and h ..
C
CBRAT(IDC)
C on channels
C        MSTU(41) = 1
        MDME(223,1) =1 !h > gamma gamma
        MDME(338,1) =1! b bbar
C.. other physics SETTINGS
        MSEL=1           ! select low pT events
        MSTU(112) = 5   ! NUMBER OF FLAVOURS 5-flavour scheme
        MSTP(3) = 1
        MSTP(32) = 4
        PARP(1) = 0.146D0   ! LAMBDA
        MDCY(6,1) =  1 !ALLOW TOP DECAY
c-- turn off some final state decays
       MDCY(PYCOMP(310),1) = 0 ! K0 not decaying
       MDCY(PYCOMP(3112),1) = 0 ! sigma-
       MDCY(PYCOMP(3122),1) = 0 ! lambda
       MDCY(PYCOMP(3222),1) = 0 ! sigma+
        MDCY( PYCOMP(3312),1) = 0 ! Xi+
        MDCY( PYCOMP(3322),1) = 0 ! Xi0
        MDCY(PYCOMP(3334),1) = 0 ! omega-


        MSTP(48) = 1  ! TOP DECAY BEFORE FRAGMINTATION
        MSTP(61) = 1! INITIAL STATE QCD AND QED RADIATION ON
        MSTP(71) = 1 !FINAL STATE QCD AND QED RADIATION ON
        MSTJ(41)=1 ! trun off QED bransstraliung, only QCD showering
        ! IF YOU WANT TO HAVE HIGGS DECAY AND HADRONISATION
        MDCY(25,1) = 1 ! ON / OFF HIGGS DECAY
        MSTP(81) = 1! MULTIPLE INTERACTION
        MSTP(111) = 2 !FREGMENTATION AND DECAY
        MSTJ(1)  = 1
        MSTJ(24)  = 2    ! particle masses from BW
C        MSTJ(11) = 3 ! PETERSON FREGMENTATION  FUNCTION FOR  b and c
C        PARJ(54) = -0.07D00  ! FREGMENTATION OF FLAVOUR c
C        PARJ(55) = -0.006    ! FREGMENTATION OF FLAVOUR b
        MSTP(42) = 0        ! OFF SHELL BOSONS PRODUCTION
        MSTP(123) = 1 ! reaction to violation of maximum
        ! QUARK NASSES
        PMAS(3,1) = 97.00D-03   ! STRANGE
        PMAS(4,1)=   1.275D0    ! CHARM
        PMAS(5,1)= 4.18D0 ! BEAUTY
        PMAS(6,1) =173.03D0    ! TOP
        !LIGHT FLAVOURS
        PMAS(2,1) = 2.225D-03 ! UP
        PMAS(1,1) =  4.725D-03   ! DOWN
        ! LEPTON MASSES
        PMAS(13,1) = 0.1056583D00 ! MUON
        PMAS(15,1) = 1.7771D00    ! TAU
C                  PARU(101) = 137.03598D00    ! 1/ALPHA QED
        ! W, Z MASSES
        PMAS(23,1) = 91.187D00   ! Z
        PMAS(24,1) = 80.330D0     ! W+-
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


C-----------------------------------------------------------------
C...Second section: event loop.
C      call pylist(12)  !call me for help
C...Generate events and look at first few.
C Some brach names
C!!!!!!!!!!!!!!!!!!!!!!!
       BMASS = 'evt_weight'  ! mass branch
       BX = 'x'     ! Production vertex x-ccordinated
       BY= 'y'      ! Production vertex y-ccordinated
       BZ ='z'      ! Production vertex z-ccordinated
       BT = 't'     ! proper lifetime
       BVX ='vx'    ! Decay vertex, x -coordinates
       BVY = 'vy'   ! Decay vertex, y -coordinates
       BVZ = 'vz'   ! Decay vertex, z -coordinates
C!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO 200 IEV=1,NEV
        CALL PYEVNT
C         CALL PYSPHE(SPH,APL)
C        CALL PYCLUS(NJET)   jet clustering
        IF(IEV.LT.2)  CALL PYLIST(7)
        IF(IEV.LT.2) CALL PYLIST(1)
C       progress bar with a time
C       call progress_bar_time(IEV,NEV)

C--   LOOP OVER PARTICLES
       DO 454 , I = 1, N, 1
               ! cut on accepetnce of the ATLAS detector
                 ETA = ATANH(P(I,3)/ DSQRT(P(I,3)**2+P(I,2)**2+P(I,1)**2))
C                 IF(ABS(ETA).GT.4.5) GOTO 454
C-- saftey check
            IF(K(I,2).NE.K(I,2)) PRINT *, "ID HIT A NAN"
            IF(P(I,1).NE.P(I,1)) PRINT *, "PX HIT A NAN", K(I,2)
            IF(P(I,2).NE.P(I,2)) PRINT *, "PY HIT A NAN", K(I,2)
             IF(P(I,3).NE.P(I,3)) PRINT *, "PZ HIT A NAN", K(I,2)
             IF(P(I,4).NE.P(I,4)) PRINT *, "E HIT A NAN", K(I,2)

C--                              PDG-ID  Mother  px     py     pz    E
            CALL ROOTADDPARTICLE(K(I,2),K(I,3),P(I,1),P(I,2),P(I,3),P(I,4))

C             CALL ROOTADDDOUBLE(V(I,1),BX)
C             CALL ROOTADDDOUBLE(V(I,2),BY)
C             CALL ROOTADDDOUBLE(V(I,3),BZ)
C             CALL ROOTADDDOUBLE(V(I,5),BT)
C             DVX = V(I,1)+V(I,5)*P(I,1)/P(I,5)
C             DVY = V(I,2)+V(I,5)*P(I,2)/P(I,5)
C             DVZ = V(I,3)+V(I,5)*P(I,3)/P(I,5)
C             CALL ROOTADDDOUBLE(DVX,BVX)
C             CALL ROOTADDDOUBLE(DVY,BVY)
C             CALL ROOTADDDOUBLE(DVZ,BVZ)
454     CONTINUE
C-         EVENT WEIGHT
C        print *,  PARI(7)*1d12 , XWGTUP*1d3
        IF(PARI(7).NE.PARI(7)) PRINT *, "DSIGMA HIT A NAN"
         CALL ROOTADDEVENT(XWGTUP*1d3)  ! convert from pb to fb ; PARI(7) is in mb
200     CONTINUE
      print*,''
      print *, " total cross section = ", pari(1)*1d12, '+-', pari(1)*1d12/sqrt(dble(nev)), '(stat)', '+-', pari(1)*1d12 *0.063, '(sys)' ,"fb"
      CALL   ROOTCLOSE
      end do
      end do
      END



CC*********************************************************************

      SUBROUTINE UPINIT

C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
C...Selection of hard scattering subprocesses.
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
C...User process initialization commonblock.

      DOUBLE PRECISION :: MH,MW,MZ,lambdaEW,lambdaNP, gammah, m
      DIMENSION  wc(6,6),fp1(-6:6),fp2(-6:6),xkf(6), xs(5)
      DOUBLE PRECISION, DIMENSION(6) :: MQ
      INTEGER IDBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,LPRUP
      DOUBLE PRECISION EBMUP,XSECUP,XERRUP,XMAXUP
         INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,IDUP(MAXNUP),
     &ISTUP(MAXNUP),MOTHUP(2,MAXNUP),ICOLUP(2,MAXNUP),PUP(5,MAXNUP),
     &VTIMUP(MAXNUP),SPINUP(MAXNUP)
      SAVE /HEPEUP/
      COMMON/HEPRUP/IDBMUP(2),EBMUP(2),PDFGUP(2),PDFSUP(2),
     &IDWTUP,NPRUP,XSECUP(1),XERRUP(1),XMAXUP(1),
     &LPRUP(1)
C
      SAVE /HEPRUP/

C...The user's own transfer of information.
      COMMON/MYCOMM/ECM,PTMIN,INIT,mod
      SAVE/MYCOMM/
      COMMON/KIN/sh,uh,th,MH,M,MQ
      COMMON/PARAM/WILSON, WC, ghff,ghhff,gammah,ghhh, lambdanp,lambdaEW,xkl,xkf
C...Local arrays and parameters.
C...Set up incoming beams. (proton proton beam )
      IDBMUP(1)= 2212
      IDBMUP(2)= 2212
      EBMUP(1) = 0.5D0*ECM  ! use CMS system
      EBMUP(2) = 0.5D0*ECM

C...Set up the external process.
      IDWTUP   = 4 ! weight info (4 IS THE BEST OPTION FOR THE CURRENT USE)
      NPRUP    = 1  ! number of external processes keep one in general
      LPRUP(1) = 303      ! a unique integer identifier of each external process
                          ! this falls into the double Higgs production
C...(This is done by a cheating call to UPEVNT with special switch.)
C calling PDF's
        PDFGUP= -1  ! do not use Pythia's PDF's
        PDFSUP = -1
        INIT=1
      CALL InitPDFsetByName("NNPDF30_lo_as_0118")  ! LHA PDF library
      CALL InitPDF(0)
      CALL UPEVNT  ! call the event generation routine with matrix element
      INIT=0
      ! no need for these variables
C      XMAXUP(1) = 5d0*XWGTUP
C        XSECUP(1) = XWGTUP
      RETURN
      END
