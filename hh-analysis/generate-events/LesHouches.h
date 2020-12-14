c -*-Fortran-*-

      integer maxpup
      parameter(maxpup=100)
      integer idbmup,pdfgup,pdfsup,idwtup,nprup,lprup
      double precision ebmup,xsecup,xerrup,xmaxup
      COMMON/HEPRUP/IDBMUP(2),EBMUP(2),PDFGUP(2),PDFSUP(2),
     &IDWTUP,NPRUP,XSECUP(1),XERRUP(1),XMAXUP(1),
     &LPRUP(1)
      integer maxnup
      parameter (maxnup=500)
      integer nup,idprup,idup,istup,mothup,icolup
      double precision xwgtup,scalup,aqedup,aqcdup,pup,vtimup,spinup
      common/hepeup/nup,idprup,xwgtup,scalup,aqedup,aqcdup,
     &              idup(maxnup),istup(maxnup),mothup(2,maxnup),
     &              icolup(2,maxnup),pup(5,maxnup),vtimup(maxnup),
     &              spinup(maxnup)
      save /hepeup/
