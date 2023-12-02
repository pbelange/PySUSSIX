! === NO CAPTAL LETTERS :)
      subroutine sussixnoo(x, xp, y, yp, s, sp,
     & tunex, tuney, tunez, amplitude, phase, ox, ax, oy, ay, os, as,
     & n_points)

!       program sussix
!======================================================================
!
! SUSSIX CONTAINS:
!
!
! ROUTINES FOR HIGH PRECISION TUNE CALCULATION:
!
!         TUNELASR (NO WINDOW)
!         TUNENEWT (HANNING WINDOW)
!         BOTH USING ZFUNR E CALCR
!
! ROUTINES FOR FREQUENCY ANALYSIS OF SIGNALS:
!
!         SPECTRUM      !SPECTRUM COMPUTATION
!
!         ORDRES        !ORDERING OF FREQUENCIES
!         READRES       !READING THE ORDERED FREQUENCIES
!
! ROUTINES FOR POST PROCESSING OF SIGNALS:
!
!         SUSRES        !SUBTRACTION OF NEXT TO LEADING FREQUENCIES
!         READSME       !SMEAR AS QUALITY FACTOR
!         READINV       !INVARIANT CALCULATION
!         READDT3       !THIRD ORDER CONJUGATING FUNCTION COEFFICIENTS
!         READDT4       !FOURTH ORDER CONJUGATING FUNCTION COEFFICIENTS
!
! ROUTINES INTERFACE WITH SIXTRACK OUTPUT
!
!         READRIC       !READS SIXTRACK TRACKING OUTPUT
!         WRITERIC      !WRITES POSTPROCESSED DATA IN SIXTRACK OUTPUT
!                        FORMAT
!
! AUTHOR: R.BARTOLINI CERN AND BOLOGNA UNIVERSITY
!
!======================================================================
!
!  MAIN PROGRAM FOR THE SPECTRAL ANALYSIS OF TRACKING OR BEAM DATA
!
!  TWO TYPES OF DATA MAY BE ANALYSED ACCORDING TO THE ISIX OPTION:
!
!  ISIX=1 ---> SIXTRACK BINARY OUTPUT,
!              input data have to be provided according to sixtrack
!              format:
!				fort.90 -> fort.59
!  ISIX=0 ---> TRACKING DATA FROM A USER PROVIDED ASCII FILE
!              input data have to be provided in the files:
!				bpm.001 -> bpm.999
!              the files must have two columns for each plane
!
!  IF ISIX=1 IT TRANSFORMS THE SIXTRACK BINARY OUTPUT INTO AN ASCII
!  FILE (SUBROUTINE READRIC),
!
!  MORE THAN ONE FILE CAN BE TREATED:
!  NTOT IS THE TOTAL NUMBER OF FILES TO BE PROCESSED
!
!  ONCE THE DATA ARE READ THE SPECTRUM IS CALCULATED, ORDERED AND
!  THE OUTPUT OF ALL THE CASES IS WRITTEN IN THE FILE lines
!  (SUBROUTINE DATSPE ---> SUBROUTINE RSPECTRUM, ORDRES, ETC.)
!
!  DIFFERENT KIND OF POSTPROCESSING ARE AVALIABLE AN ARE SWITCHED
!  ON WITH THE CORRESPONDING FLAG, BY READING FILE lines WITH
!  THE SUBROUTINE  READRES.
!  NSUS GE 1 ---> SUBTRACTS THE NEXT TO LEADING FREQUENCIES
!  (SUBROUTINE SUSRES). IN THIS CASE THE PROGRAM ALSO WRITES THE
!  MODIFIED DATA BACK INTO THE STARTING SIXTRACK BYNARY OUTPUT
!  (SUBROUTINE WRITERIC). WARNING: IT OVERWRITES THE FILES.
!  ISME=1 ---> SMEAR CALCULATION
!  (SUBROUTINE READSME)
!  INV=1 ---> INVARIANT CALCULATION
!  (SUBROUTINE READINV)
!
!  N.B.: IF THE OUTPUT IN FILE lines IS ALREADY AVALIABLE THE OPTION
!  IANA=0 ALLOWS TO SKIP DATA ANALYSIS, STARTING DIRECTLY WITH THE
!  lines ANALYSIS.
!
!  DESCRIPTION OF ALL INPUT ITEMS:
!
!  ISIX FLAG FOR SIXTRACK (1) OR ASCII (0) DATA
!  NTOT TOTAL NUMBER OF DATA FILES TO BE ANALYZED
!  IANA FLAG FOR FULL ANALYSIS OF DATA
!  ICONV IS THE FLAG FOR LINEAR TRASFORMATION (YES=1)
!  NT1,NT2 initial & final turn number
!  NARM THE MUNBER OF HARMONIC TO BE CALCULATED
!  ISTUNE IS THE FLAG FOR FUNDAMENTALFREQUENCIES
!  ISTUNE=0 => off
!  ISTUNE=1 uses Qx,y,z values as guess values
!  ISTUNE=2 takes the Qx,y,z values as fundamental frequencies
!  etune(3) allowed distance to Qx,y,z
!  TUNEX,TUNEY,TUNEZ guess or fundamental tunes
!  NSUS THE NUMBER OF HARMONIC TO BE SUBTRACTED
!  IDAM IS THE DIMENSION OF PHASE SPACE
!  NTWIX IS A FLAG FOR THE TWIN PARTICLES (1 or 2)
!  IR IS A FLAG FOR REAL SIGNAL (YES=1)
!  IMETH IS A FLAG FOR WINDOWING (HANNING=1, NO FILTER=0)
!  NRC IS THE MAXIMUM ORDER OF LINEAR COMBINATION OF FREQUENCIES
!  EPS IS THE TOLERANCE ON THE IDENTIFICATION OF FREQUENCIES
!  NLINE IS THE NUMBER OF LINES TO BE LOOKED FOR
!  LR MR KR SPECIFY THE LINE
!  IDAMX SELECT THE PLANE TO ANALYZE
!  IFIN IS THE UNIT FOR THE FINAL OUTPUT
!  ISME FLAG FOR SMEAR CALCULATION
!  IUSME UNIT FOR SMEAR OUTPUT
!  INV FLAG FOR INVARIANTS CALCULATION
!  IINV UNIT FOR INVARIANTS OUTPUT
!  ICF FLAG FOR INVARIANTS CALCULATION
!  IICF UNIT FOR INVARIANTS OUTPUT
!
!  THE SPECTRUM IS WRITTEN TO FILE lines
!
!  reserved file units:
!       10 input file (filename: sussix.input)
!       30 main output (filename: lines)
!       90 -> 59 SIXTRACK input data
!       91/92 SIXTRACK aux files (filename: aux1/aux2)
!       90 ascii input (filename: bpmXXX)
! the other output units can be chosen by the user in the input
! files. defaults are:
!       50 -> 50+nline spectral line ordered (filename: resonXX)
!       20 for smear calculation (filename: smearXX)
!       25 for invariant calculation (filename: invarXX)
!       35 for conjugating function calculation (filename: conjuXX)
!
! AUTHORS: R.BARTOLINI/F.SCHMIDT -- DIAMOND/CERN
! Copyright (C) Riccardo Bartolini and Frank Schmidt
!
!  MODIFIED 16/09/1996: ADDED THE ICONV OPTION
!  MODIFIED 17/09/1996: ADDED THE INPUT FROM FILE
!  MODIFIED 20/08/1997: ADDED THE INVARIANT CALCULATION
!  MODIFIED 17/09/1997: ADDED THE TREATEMENT OF ASCII BEAM DATA
!  MODIFIED 30/05/1998: ADDED DEFAULT VALUES, AND FFT
!  MODIFIED 12/10/1999: TOTAL CLEAN-UP
!  LAST MODIFIED:16/01/2005: Portability F77/F90, IO <=99
!
! ==============================================================================
! 2012 KEVIN: CHANGES MADE:
! FORWARD FILE OUTPUT TO OUTPUT VARIABLES
! ADDED INPUT VARIABLES: N_POINTS, X, XP, Y, YP, S, SP
! ADDED PUTPUT VARIABLES: AMPLITUDE, PHASE, OX, AX, OY, AY, OS, AS
! ==============================================================================


      implicit none
      integer i,iana,icf,iconv,idam,idamx,ifi,ifin,iicf,iinv,imeth,ini, &
     &inv,iouk,ir,isix,isme,istune,iunit,iusme,k,kr,lr,mr,n,narm,nf,    &
     &nline,nlst,nrc,nsus,nt1,nt2,ntot,nturn,ntwin,ntwix
      double precision eps,etune,tunex,tuney,tunez
      integer mterm
      parameter(mterm=300)
      integer mbpm
      double precision pieni,zero,one,two,pi,duepi,piph
      parameter(mbpm=1000,pieni=1d-17,zero=0d0,one=1d0,two=2d0)
      double precision tsa,txa,tya
      common/tune/txa(mterm),tya(mterm),tsa(mterm)
      dimension lr(100),mr(100),kr(100),etune(3)
      character*200 ch,ch1
      character*8 filename


      integer n_points
      double precision x(n_points), xp(n_points),
     & y(n_points), yp(n_points), s(n_points), sp(n_points),
     & ox(600), ax(600), oy(600), ay(600), os(600), as(600),
     & xyz(6*n_points), tunexyz(3), amplitude(14), phase(14)

Cf2py intent(out) tunex
Cf2py intent(out) tuney
Cf2py intent(out) tunez
Cf2py intent(out) amplitude
Cf2py intent(out) phase
Cf2py intent(out) ox
Cf2py intent(out) ax
Cf2py intent(out) oy
Cf2py intent(out) ay
Cf2py intent(out) os
Cf2py intent(out) as
Cf2py depend(n_points) x, xp, y, yp, s, sp

!===================
!.....INITIALIZATION
!===================
!=========================
!.....DEFAULT VALUES FIRST
!=========================
      isix=1
      ntot=1
      iana=1
      iconv=1
      nt1=1
      nt2=1024
      narm=1
      istune=0
      etune(1)=1d-2
      etune(2)=1d-2
      etune(3)=1d-2
      tunex=0.28
      tuney=0.31
      tunez=0.006
      nsus=0
      idam=2
      ntwix=1
      ir=0
      imeth=2
      nrc=10
      eps=1d-6
      nline=1
      do i=1,100
        lr(i)=0
        mr(i)=0
        kr(i)=0
      enddo
      lr(1)=1
      idamx=1
      ifin=50
      isme=0
      iusme=20
      inv=0
      iinv=25
      icf=0
      iicf=35
      do i=1,mterm
        txa(i)=zero
        tya(i)=zero
        tsa(i)=zero
      enddo
!=========================
!.....READS THE INPUT FILE
!=========================
 
      open(10,file='sussix.inp',form='formatted',status='unknown')
      read(10,'(4(/))')
      read(10,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)isix
      read(10,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)ntot
      read(10,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)iana
      read(10,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)iconv
      read(10,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)nt1,nt2
      read(10,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)narm
      read(10,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)istune,etune(1),etune(2),etune(3)
      read(10,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)tunex,tuney,tunez
      read(10,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)nsus
      read(10,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)idam
      read(10,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)ntwix
      read(10,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)ir
      read(10,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)imeth
      read(10,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)nrc
      read(10,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)eps
      read(10,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)nline
      do k=1,nline
        read(10,'(A)') ch
        ch1=ch(9:80)//' / '
        read(ch1,*) lr(k),mr(k),kr(k)
      enddo
      read(10,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)idamx
      read(10,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)ifin
      read(10,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)isme
      read(10,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)iusme
      read(10,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)inv
      read(10,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)iinv
      read(10,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)icf
      read(10,'(A)') ch
      ch1=ch(9:80)//' / '
      read(ch1,*)iicf
 
!.....CHECKS
      if(narm.le.0) then
        write(6,*) 'NARM too small'
        close(10)
        stop
      endif
      if(narm.gt.mterm) then
        write(6,*) 'NARM too big => reduced to maximum: ',mterm
        narm=mterm
      endif
      if(nt1.le.0) then
        write(6,*) 'NT1 too small'
        close(10)
        stop
      endif
      if(nt2.le.nt1) then
        write(6,*) 'NT2 smaller than NT1'
        close(10)
        stop
      endif
      if(idam.ne.1.and.idam.ne.2.and.idam.ne.3) then
        write(6,*) 'The order of phase space IDAM must be 1, 2 or 3'
        close(10)
        stop
      endif
      if(isix.eq.1) then
        if(ntot.gt.32) then
          write(6,*) 'NTOT too big => reduced to maximum: ',32
          ntot=32
        endif
      endif
      if(ntwix.ne.1.and.ntwix.ne.2) then
        write(6,*) 'NTWIX ill defined, set to: ',1
        ntwix=1
      endif
 
!==============================
!.....END OF THE INITIALIZATION
!==============================
!===========================
!.....STARTING DATA ANALYSIS
!===========================
 
!.....1) CHECK IF FILE lines IS ALREADY PRODUCED
!.....IF NOT (IANA=1 or 2)
!........2) CHECK THE ISIX OPTION
!........3) PROCEED WITH SPECTRUM CALCULATION
!........4) FFT IF IANA=2
!.....IF YES (IANA=0)
!........2) PROCEED WITH POSTPROCESSING


! ==============================================================================
! FOR NoO: commented all output of files '6' (screen) and '30'
! ==============================================================================

 
!       open(30,file='lines',form='formatted',status='unknown')
      if(iana.eq.1.or.iana.eq.2) then
        if(isix.eq.1) then
!@@@@@@@@@@@@@@@@@@@@@@
!.....SIXTRACK DATA   @
!@@@@@@@@@@@@@@@@@@@@@@
          nlst=90-ntot+1
          open(91,file='aux1',form='formatted',status='unknown')
          open(92,file='aux2',form='formatted',status='unknown')
          do n=90,nlst,-1
            if(n.lt.59) goto 100
            filename='fort.'
            write(filename(6:7),'(i2.2)') n
            open(n,file=filename,form='unformatted',status='unknown')
            call readric(n,idam,ntwin,iconv)
!.....NOW READRIC HAS CREATED NTWIN ASCII FILES aux1 AND
!.....EVENTUALLY aux2 NTWIN IS AN OUTPUT PARAMETER OF READRIC!!
            if(ntwix.lt.ntwin) then
!               write(6,*)'WARNING: THE TWIN PARTICLE IS IGNORED'
              if(nsus.ge.2) then
!                 write(6,*)'ERROR: TWIN PARTICLE NEEDED'
                close(10)
!                 close(30)
                stop
              endif
            endif
            do nf=1,ntwix
              iunit=90+nf
              call datspe(iunit,idam,ir,nt1,nt2,nturn,imeth,narm,iana,
     &         x,xp,y,yp,s,sp,n_points)
!.....N.B. NTURN IS AN OUTPUT PARAMETER OF DATSPE
              call ordres(eps,narm,nrc,idam,iunit,nturn,
     & -tunex,-tuney,-tunez,istune,etune,
     & amplitude,phase,ox,ax,oy,ay,os,as)
              if(nsus.ge.2) then
!.....SUBTRACTS AND OVERWRITE UNIT IUNIT
                call susres(iunit,nsus,nturn,3)
              endif
            enddo
 
            if(nsus.ge.2) then
!.....N.B. WRITERIC NEEDS BOTH THE FILES TREATED WITH SUSRES
!.....     IT OVERWRITES THE INITIAL SIXTRACK OUTPUT!!!!!!!!
              call writeric(n,ntwin,nturn,iconv)
            endif
          enddo
        elseif(isix.eq.0) then
!@@@@@@@@@@@@@@@@@@@@
!......ASCII DATA   @
!@@@@@@@@@@@@@@@@@@@@
          do n=1,ntot
            iunit=90
            filename='bpm.'
            write(filename(5:8),'(i3.3)') n
!             open(iunit,file=filename,form='formatted',status='unknown')
            call datspe(iunit,idam,ir,nt1,nt2,nturn,imeth,narm,iana,
     &       x,xp,y,yp,s,sp,n_points)
            call ordres(eps,narm,nrc,idam,n,nturn,
     & -tunex,-tuney,-tunez,istune,etune,
     & amplitude,phase,ox,ax,oy,ay,os,as)
!             close(iunit)
          enddo
        endif
      endif
 100  continue
      close(91)
      close(92)
!=============================================
!.....STARTING POSTPROCESSING OF ANALYZED DATA
!=============================================
 
!......1) SELECTION OF LINES
!......2) SMEAR CALCULATION
!......3) INVARIANT CALCULATION
 
      if(isix.eq.1) then
!.......NTOT*NTWIX CASES ANALYZED FROM FILE lines
        ini=90
        ifi=ini-ntot*ntwix+1
      elseif(isix.eq.0) then
!.......NTOT CASES ANALYZED FROM FILE lines
        ini=1
        ifi=ntot
      endif
 
!       write(6,*)' '
!       write(6,*)'************************************* '
!       write(6,*)'STARTING FILE lines ANALYSIS OF CASES:',ini,ifi
!       write(6,*)'************************************* '
 
!.....SELECTION OF NLINE LINES SPECIFIED IN THE ARRAY LR,MR,KR
!.....THE RESULTS ARE WRITTEN IN THE UNIT IOUK, IDAMX IS THE
!.....PLANE TO BE ANALYZED
      do k=1,nline
        iouk=ifin+k
        if(iouk.gt.99) then
!           write(6,*) 'Unit for resonance output  must be below 100'
          close(10)
!           close(30)
          stop
        endif
        filename='reson'
        write(filename(6:8),'(i2.2)') iouk
        open(iouk,file=filename,form='formatted',status='unknown')
        call readres(ini,ifi,lr(k),mr(k),kr(k),narm,                    &
     &idam,idamx,iouk)
        close(iouk)
      enddo
!.....SMEAR CALCULATION: OUTPUT IN THE UNIT IUSME
      if(isme.eq.1) then
!         write(6,*)'SMEAR CALCULATION: OUTPUT IN THE UNIT         ',iusme
        if(iusme.gt.99) then
!           write(6,*) 'Unit for smear calculation must be below 100'
          close(10)
!           close(30)
          stop
        endif
        filename='smear'
        write(filename(6:8),'(i2.2)') iusme
        open(iusme,file=filename,form='formatted',status='unknown')
        call readsme(ini,ifi,narm,idam,idamx,iusme)
        close(iusme)
      endif
!.....INVARIANT CALCULATION: OUTPUT IN THE UNIT IINV
      if(inv.eq.1) then
        if(iinv.gt.99) then
!           write(6,*) 'Unit for invariant calculation must be below 100'
          close(10)
!           close(30)
          stop
        endif
        filename='invar'
        write(filename(6:8),'(i2.2)') iinv
        open(iinv,file=filename,form='formatted',status='unknown')
!         write(6,*)'INVARIANT CALCULATION: OUTPUT IN THE UNIT    ',iinv
        call readinv(ini,ifi,narm,idam,iinv)
        close(iinv)
      endif
!.....3-RD ORDER CONJUGATING FUNCTION CALCULATION: OUTPUT IN THE UNIT ICF
      if(icf.eq.1) then
!         write(6,*)'3-RD ORDER CONJUG. FUNC.: OUTPUT IN THE UNIT',iicf
        if(iicf.gt.99) then
!           write(6,*) 'Unit for conjugating function must be below 100'
          close(10)
!           close(30)
          stop
        endif
        filename='conju'
        write(filename(6:8),'(i2.2)') iicf
        open(iicf,file=filename,form='formatted',status='unknown')
        call readdt3(ini,ifi,narm,idam,iicf)
        call readdt4(ini,ifi,narm,idam,iicf)
        close(iicf)
      endif
      close(10)
!       close(30)
      return
      end
      subroutine spectrum(x,xp,maxn,tune,zpesi,narm,meth)
!=======================================================================
!
! SUBROUTINE SPECTRUM
!
! COMPUTE THE MAIN FREQUENCY
! X, XP ARE THE COORDINATES OF THE ORBIT AND MAXN IS THE
! LENGTH IF THE ORBIT.
! WITHOUT ORTHOGONALIZATION OF GRAM-SCHMIDT
!
! METH SELECTS THE WINDOW:
!     1 --> HANNING WINDOW
!     2 --> RECTANGULAR WINDOW
!
! AUTHOR:  R. BARTOLINI 9/1/1996
!          A. BAZZANI
!
!=======================================================================
      implicit none
      integer maxn,meth,n,na,narm
      double precision freq,tune,tunelasr,tunenewt,x,xp
      double complex z,zef,zgs,zpesi,zw,zx,zz
      integer mterm
      parameter(mterm=300)
      integer maxiter
      parameter(maxiter=100000)
      integer mbpm
      double precision pieni,zero,one,two,pi,duepi,piph
      parameter(mbpm=1000,pieni=1d-17,zero=0d0,one=1d0,two=2d0)
      dimension x(maxiter),xp(maxiter)
      dimension z(maxiter),zz(maxiter)
      dimension tune(mterm),zpesi(mterm),zgs(maxiter)
 
!===============
! INITIALIZATION
!===============
      duepi=8d0*atan(1d0)
 
      do n=1,maxiter
        z(n)=dcmplx(0d0,0d0)
        zz(n)=z(n)
      enddo
      zw=dcmplx(0d0,0d0)
      if(maxn.gt.maxiter) then
        write(6,*) 'ERROR IN SPECTRUM: MAXN TOO LARGE'
        close(10)
        close(30)
        stop
      endif
      if(narm.lt.2) then
        write(6,*) 'ERROR IN SPECTRUM: NA SMALLER THAN 2'
        close(10)
        close(30)
        stop
      endif
 
      do n=1,maxn
        z(n)=dcmplx(x(n),xp(n))
        zz(n)=z(n)
      enddo
 
      do na=1,narm
        if(meth.eq.1) then
          tune(na)=tunenewt(x,xp,maxn,zw)
        elseif(meth.eq.2) then
          tune(na)=tunelasr(x,xp,maxn,zw)
        endif
!.....beginning of subtraction procedure
        freq=tune(na)
        zpesi(na)=zw/dble(maxn)
        zef=exp(dcmplx(0d0,freq*duepi))
        zx=1
        zgs(1)=zpesi(na)*zx
 
        do n=2,maxn
          zx=zx*zef
          zgs(n)=zpesi(na)*zx
        enddo
 
        do n=1,maxn
          z(n)=z(n)-zgs(n)
        enddo
        do n=1,maxn
          x(n)=dble(z(n))
          xp(n)=dimag(z(n))
        enddo
 
      enddo
 
!.....restore the original signal in x, xp.
      do n=1,maxn
        x(n)=dble(zz(n))
        xp(n)=dimag(zz(n))
      enddo
 
      return
      end
      double precision function tunenewt(x,xp,maxn,zw)
!=======================================================================
!
! SUBROUTINE TUNENEWT
!
! COMPUTES THE TUNE USING A DISCRETE VERSION OF LASKAR METHOD.
! IT INCLUDES A NEWTON METHOD FOR THE SEARCH OF THE FREQUENCY.
! X, XP ARE THE COORDINATES OF THE ORBIT AND MAXN IS THE
! LENGTH IF THE ORBIT.
!
! AUTHOR:     A. BAZZANI - BOLOGNA UNIVERSITY
!             R. BARTOLINI - CERN HAS INTRODUCED SOME MODIFICATIONS
!
!=======================================================================
      implicit none
      integer i,maxn,maxn2,mf,mft,nft,nftmax,npoint
      double precision deltat,ftmax,step,tune,tune1,tunefou,x,xp
      double complex z,zw
      complex zsing
      integer maxiter
      parameter(maxiter=100000)
      integer mbpm
      double precision pieni,zero,one,two,pi,duepi,piph
      parameter(mbpm=1000,pieni=1d-17,zero=0d0,one=1d0,two=2d0)
      dimension x(maxiter),xp(maxiter),zsing(maxiter)
      dimension z(maxiter)
 
! INITIALIZATION
      duepi=atan(1d0)*8d0
      tune=0d0
 
!.............................................................
!    ESTIMATION OF TUNE WITH FFT
!.............................................................
      mft=int(log(dble(maxn))/log(2d0))
      npoint=2**mft
      maxn2=maxn/2
      step=duepi/maxn
      do mf=1,maxn
        z(mf)=dcmplx(x(mf),xp(mf))*(1d0+cos(step*(mf-maxn2)))
        zsing(mf)=z(mf)
      enddo
      call cfft(zsing,-mft)
!.............................................................
!   SEARCH FOR MAXIMUM OF FOURIER SPECTRUM
!.............................................................
      ftmax=0d0
      nftmax=0
      do nft=1,npoint
        if (abs(zsing(nft)).gt.ftmax) then
          ftmax=abs(zsing(nft))
          nftmax=nft
        end if
      enddo
      tunefou=dble(nftmax-1)/dble(npoint)
      if(tunefou.ge.0.5d0) tunefou=-(1d0-tunefou)
      deltat=1d0/npoint
      tune1=tunefou-deltat
      call zfunr(tune,zw,z,maxn,tune1,deltat)
      tunenewt=tune
 
!............................................................
      return
!............................................................
      end
      double precision function tunelasr(x,xp,maxn,zw)
!=======================================================================
!
! SAME AS TUNENEWT BUT NO HANNING FILTER
!
! AUTHOR:     A. BAZZANI - BOLOGNA UNIVERSITY
!             R. BARTOLINI - CERN HAS INTRODUCED SOME MODIFICATIONS
!
!=======================================================================
      implicit none
      integer i,maxn,maxn2,mf,mft,nft,nftmax,npoint
      double precision deltat,ftmax,step,tune,tune1,tunefou,x,xp
      double complex z,zw
      complex zsing
      integer maxiter
      parameter(maxiter=100000)
      integer mbpm
      double precision pieni,zero,one,two,pi,duepi,piph
      parameter(mbpm=1000,pieni=1d-17,zero=0d0,one=1d0,two=2d0)
      dimension x(maxiter),xp(maxiter),zsing(maxiter)
      dimension z(maxiter)
 
! INITIALIZATION
      duepi=atan(1d0)*8d0
      tune=0d0
      zw=dcmplx(0d0,0d0)
 
!.............................................................
!    ESTIMATION OF TUNE WITH FFT
!.............................................................
      mft=int(log(dble(maxn))/log(2d0))
      npoint=2**mft
      maxn2=maxn/2
      step=duepi/maxn
      do mf=1,maxn
        z(mf)=dcmplx(x(mf),xp(mf))  ! no filter
        zsing(mf)=z(mf)
      enddo
      call cfft(zsing,-mft)
!.............................................................
!   SEARCH FOR MAXIMUM OF FOURIER SPECTRUM
!.............................................................
      ftmax=0d0
      nftmax=0
      do nft=1,npoint
        if (abs(zsing(nft)).gt.ftmax) then
          ftmax=abs(zsing(nft))
          nftmax=nft
        end if
      enddo
      tunefou=dble(nftmax-1)/dble(npoint)
      if(tunefou.ge.0.5d0) tunefou=-(1d0-tunefou)
      deltat=1d0/npoint
      tune1=tunefou-deltat
      call zfunr(tune,zw,z,maxn,tune1,deltat)
      tunelasr=tune
 
      return
      end
      subroutine zfunr(tune,zw,z,maxn,tunea1,deltat)
!=======================================================================
! AUXILIARY ROUTINE USED BY TUNENEWT.
!
! AUTHOR:     A. BAZZANI - BOLOGNA UNIVERSITY
!
!=======================================================================
      implicit none
      integer maxn,nc,ncont,nd,ntest,num
      double precision deltat,dtune1,dtune2,dtune3,dtunea1,dtunea2,     &
     &err,ratio,tune,tune1,tune2,tune3,tunea1,tunea2,tunetest,tuneval,  &
     &tunevmax
      double complex z,zd,zf,zfd,ztune,ztune1,ztune2,ztune3,zu,zw
      integer maxiter
      parameter(maxiter=100000)
      integer mbpm
      double precision pieni,zero,one,two,pi,duepi,piph
      parameter(mbpm=1000,pieni=1d-17,zero=0d0,one=1d0,two=2d0)
      dimension z(maxiter),zd(maxiter),tunetest(10),tuneval(10)
 
! INITIALIZATION
      duepi=atan(1d0)*8d0
      err=1d-10
      zu=dcmplx(0d0,1d0)
 
!............................................................
!.... WE DIVIDE DELTAT IN 5 PARTS
!............................................................
      deltat=deltat/5.d0
!............................................................
      do nd=1,maxn
        zd(nd)=zu*nd*z(nd)
      enddo
      do nd=1,10
        tunetest(nd)=0d0
        tuneval(nd)=0d0
      enddo
!............................................................
      ztune1=exp(-zu*duepi*tunea1)
      call calcr(ztune1,zf,z,maxn)
      call calcr(ztune1,zfd,zd,maxn)
      dtunea1=dble(zf)*dble(zfd)+dimag(zf)*dimag(zfd)
      num=1
      do ntest=1, 10
        tunea2=tunea1+deltat
        ztune2=exp(-zu*duepi*tunea2)
        call calcr(ztune2,zf,z,maxn)
        call calcr(ztune2,zfd,zd,maxn)
        dtunea2=dble(zf)*dble(zfd)+dimag(zf)*dimag(zfd)
        if ((dtunea1.le.0d0).and.(dtunea2.ge.0d0)) then
          tune1=tunea1
          tune2=tunea2
          dtune1=dtunea1
          dtune2=dtunea2
          do ncont=1,100
!            ratio=-dtune1/dtune2
            if(abs(dtune2).gt.0) then
              ratio=-dtune1/dtune2
            else
              ratio=0d0
            endif
            tune3=(tune1+ratio*tune2)/(1.d0+ratio)
            ztune3=exp(-zu*duepi*tune3)
            call calcr(ztune3,zf,z,maxn)
            call calcr(ztune3,zfd,zd,maxn)
            dtune3=dble(zf)*dble(zfd)+dimag(zf)*dimag(zfd)
            if (dtune3.le.0d0) then
              if(tune1.eq.tune3) goto 100
              tune1=tune3
              dtune1=dtune3
            else
              if(tune2.eq.tune3) goto 100
              tune2=tune3
              dtune2=dtune3
            endif
            if (abs(tune2-tune1).le.err) goto 100
          enddo
 100      tunetest(num)=tune3
          tuneval(num)=cdabs(zf)
          num=num+1
        endif
        tunea1=tunea2
        dtunea1=dtunea2
      enddo
      tune=tunetest(1)
      tunevmax=tuneval(1)
      do nc=2, num-1
        if(tunevmax.le.tuneval(nc)) then
          tunevmax=tuneval(nc)
          tune=tunetest(nc)
        endif
      enddo
      ztune=exp(-zu*duepi*tune)
      call calcr(ztune,zw,z,maxn)
!............................................................
      return
!............................................................
      end
      subroutine calcr(zv,zpp,zp,maxd)
!=======================================================================
! AUXILIARY ROUTINE USED BY TUNENEWT.
!
! AUTHOR:     A. BAZZANI - BOLOGNA UNIVERSITY
!
!=======================================================================
      implicit none
      integer maxd,np
      double complex zp,zpp,zv
      dimension zp(*)
      zpp=zp(maxd)
!............................................................
      do np=maxd-1,1, -1
        zpp=zpp*zv+zp(np)
      enddo
!............................................................
      return
!............................................................
      end
      subroutine susres(nfile,nsus,max,idams)
!=======================================================================
!
! SUBROUTINE SUSRES
!
! SUBTRACTS TO THE SIGNALS THE NEXT TO LEADING FREQUENCIES
! LEAVING THE TUNE IN.
!
! NSUS = HARMONICS TO BE SUBTRACTED
! MAX  = LENGTH OF THE SIGNAL
! NFILE = THE OUTPUT FILE WITH THE SUBTRACTED SIGNAL
!
! THIS ROUTINE MUST BE CALLED ONLY AFTER THE EXECUTION OF
! THE ROUTINE DATSPE WHICH FILLS THE COMMONS
!
! AUTHOR: R.BARTOLINI
! LAST MODIFIED: 01/04/1996
!
!=======================================================================
      implicit none
      integer idams,j,k,max,nfile,nsus
      double precision ss,ssp,xx,xxp,yy,yyp
      double complex zpots,zpotx,zpoty,zts,ztx,zty,zxs,zys,zss
      integer mterm
      parameter(mterm=300)
      integer maxiter
      parameter(maxiter=100000)
      integer mbpm
      double precision pieni,zero,one,two,pi,duepi,piph
      parameter(mbpm=1000,pieni=1d-17,zero=0d0,one=1d0,two=2d0)
      dimension zxs(maxiter),zys(maxiter),zss(maxiter)
      double precision s,sp,x,xp,y,yp
      common/data/x(maxiter),y(maxiter),xp(maxiter),yp(maxiter),        &
     &s(maxiter),sp(maxiter)
      double precision tsa,txa,tya
      common/tune/txa(mterm),tya(mterm),tsa(mterm)
      double complex zspes,zxpes,zypes
      common/fcoe/zxpes(mterm),zypes(mterm),zspes(mterm)
 
! INITIALIZATION
      duepi=8d0*atan(1d0)
 
      write(6,*)'STARTING THE SUBTRACTION PROCEDURE IN FILE',nfile
 
!....BUILD THE COMPLEX SIGNAL Z = X + i PX
 
      do j=1,max
        zxs(j)=dcmplx(x(j),xp(j))
        zys(j)=dcmplx(y(j),yp(j))
        zss(j)=dcmplx(s(j),sp(j))
      enddo
 
      do j=2,nsus
        zpotx=1d0
        ztx=exp(dcmplx(0d0,txa(j)*duepi))
        zpoty=1d0
        zty=exp(dcmplx(0d0,tya(j)*duepi))
        zpots=1d0
        zts=exp(dcmplx(0d0,tsa(j)*duepi))
        do k=1,max
          zxs(k)=zxs(k)-zxpes(j)*zpotx
          zys(k)=zys(k)-zypes(j)*zpoty
          zss(k)=zss(k)-zspes(j)*zpots
          zpotx=zpotx*ztx
          zpoty=zpoty*zty
          zpots=zpots*zts
        enddo
      enddo
 
!....WRITES THE REMAINING SIGNAL
 
      do k=1,max
        xx=dble(zxs(k))
        xxp=dimag(zxs(k))
        yy=dble(zys(k))
        yyp=dimag(zys(k))
        ss=dble(zss(k))
        ssp=dimag(zss(k))
        if(idams.eq.1) then
          write(nfile,'(2(G21.15,1X))')xx,xxp
        elseif(idams.eq.2) then
          write(nfile,'(4(G21.15,1X))')xx,xxp,yy,yyp
        elseif(idams.eq.3) then
          write(nfile,'(6(G21.15,1X))')xx,xxp,yy,yyp,ss,ssp
        endif
      enddo
 
      return
 
      end
      subroutine ordres(eps,narm,nr,idam,iunit,nturn,tunex,
     & tuney,tunez,istune,etune,
     & amplitude,phase,ox,ax,oy,ay,os,as)
!=======================================================================
!
!  ORDERS THE HARMONICS FOUND BY SPECTRUM
!
!  INPUT PARAMETERS:
!
!  NARM : NUMBER OF HARMONICS TO BE OREDERED
!  NR   : MAXIMUM ORDER OF HARMONIC TO BE LOOKED FOR IN THE LINEAR
!         COMBINATIONS. (NR<10 ALTRIMENTI L'OUTPUT PUO INCASINARSI)
!  EPS  : MAXIMUM ERROR ACCEPTED TO FIND THE LINEAR COMBINATIONS.
!  TUNEX,TUNEZ,TUNES ARE THE EXPECTED TUNES WITHIN 0.01
!
!  MIND :  THE INPUT DATA COME FROM THE COMMON, SO THIS ROUTINE
!          MUST BE USED ONLY AFTER THE CALL TO DATSPE
!          THE ARRAYS TX,TY,TZ ARE USED IN ORDER NOT TO CHANGE
!          THE ARRAYS TXA,TYA,TSA.
!
!  THE OUTPUT IS PLACED IN THE FILE lines WHICH IS NOT CLOSED AT THE
!  END OF THE SUBROUTINE IN ORDER TO COLLECT THE RESULTS OF DIFFERENT
!  SIGNALS IN ONE FILE ONLY. THE DIFFERENT SIGNALS ARE NUMBERED BY IOU.
!
!  A CHECK IS PERFORMED TO SEE IF THE HARMONICS ARE CLOSER THAN 1/NTURN
!
!  IOU  : INDEX WHICH IDENTIFIES THE CASE ANALIZED
!
!  AUTHOR: R.BARTOLINI 12/02/1996
!  LAST MODIFIED: 21/09/1997
!
!=======================================================================
! ==============================================================================
! 2012 KEVIN: CHANGES MADE:
! ADDED OUTPUT VARIABLES: AMPLITUDE, PHASE, OX, AX, OY, AY, OS, AS
! COMMENTED MOST WRITE(30)
! REPLACED ALL write(30,100) -> CORRESPONDING VARIABLE ASSIGNMENTS I.E. OX, AX
! INSERTED SOME LINES FROM ROGELIO'S SUSSIX4DRIVEXX
! ==============================================================================


      implicit none
      integer idam,imiss,imissx,imissy,imissz,isca,iscax,iscay,iscaz,   &
     &isearx,iseary,isearz,istune,iunit,j,j1,k,k1,l,l1,m,m1,n,narm,     &
     &narm2,nr,nt,nturn,ntx,nty,ntz
      double precision az,checkn,dt,dtunex,dtuney,dtunez,dtxz,dty,dtyz, &
     &eps,epsx,epsy,epsz,etune,ex,ey,ez,fx,fxt,fy,fyt,fz,fzt,ordc,ordcx,&
     &ordcy,ordcz,px,pxt,pxti,pxtr,py,pyt,pyti,pytr,pz,pzt,pzti,pztr,   &
     &tunex,tuney,tunez,tx,txt,ty,tyt,tz,tzt
      double complex zpx,zpy,zpz
      integer mterm
      parameter(mterm=300)
      integer mbpm
      double precision pieni,zero,one,two,pi,duepi,piph
      parameter(mbpm=1000,pieni=1d-17,zero=0d0,one=1d0,two=2d0)
      double precision tsa,txa,tya
      common/tune/txa(mterm),tya(mterm),tsa(mterm)
      double complex zspes,zxpes,zypes
      common/fcoe/zxpes(mterm),zypes(mterm),zspes(mterm)
      dimension tx(mterm),ty(mterm),tz(mterm),etune(3)

      ! insertion P.Belanger
      !--------------------
      double precision j_vec,k_vec,l_vec,m_vec,order_vec
      common/jklm/j_vec(mterm),k_vec(mterm),l_vec(mterm),m_vec(mterm)
      common/jklm/order_vec(mterm)
      !--------------------

      double precision amplitude(14), phase(14), ox(600), ax(600),
     & oy(600), ay(600), os(600), as(600)

 
! INITIALIZATION
      pi=4d0*atan(1d0)
      checkn=1d0/dble(nturn)
      imissx=0
      imissy=0
      imissz=0
      iscax=0
      iscay=0
      iscaz=0
 
!       if(nr.gt.10) then
! !         write(6,*)'ERROR IN ORDRES: NR LARGER THAN 10'
! !         close(10)
!         close(30)
! !         close(iunit)
!         stop
!       endif
 
      do j=1,narm
        tx(j)=txa(j)
        ty(j)=tya(j)
        tz(j)=tsa(j)
      enddo
 
!
!.. TUNES PARAMETERS AND EVENTUAL CHECK FOR THE EXPECTED TUNES
!
 
      if(istune.ge.1) then
!.....CHECK X TUNE
        dtunex=abs(abs(tx(1))-abs(tunex))
        if(dtunex.gt.etune(1).or.(tx(1)*tunex).lt.0) then
          write(6,*)'X TUNE DIFFERENT FROM EXPECTED'
          write(6,*)-tx(1),-tunex
          ntx=1
          do nt=2,narm
            dtunex=abs(abs(tx(nt))-abs(tunex))
            if(dtunex.le.etune(1).and.(tx(nt)*tunex).gt.0) then
              ntx=nt
              goto 7
            endif
          enddo
 7        if(ntx.gt.1) then
            write(6,*)'EXPECTED TUNE X FOUND AT LINE',ntx
          elseif(ntx.eq.1) then
            write(6,*)'EXPECTED TUNE X NOT FOUND'
            if(istune.eq.1) write(6,*)'LINE 1 ASSUMED AS TUNE!!!'
          endif
          tx(1)=tx(ntx)
          tx(ntx)=txa(1)
          zpx=zxpes(1)
          zxpes(1)=zxpes(ntx)
          zxpes(ntx)=zpx
        endif
      endif
 
      if(istune.eq.2) then
        txt=tunex
      else
        txt=tx(1)
      endif
      pxt=abs(zxpes(1))
      if (abs(pxt).gt.0) then
        pxtr=dble(zxpes(1))/pxt
        pxti=dimag(zxpes(1))/pxt
        if(pxti.eq.zero.and.pxtr.eq.zero) then
          fxt=zero
        else
          fxt=atan2(pxti,pxtr)
        endif
      else
        pxtr=0d0
        pxti=0d0
        fxt=0d0
      endif
      fxt=fxt/pi*1.8d2
      tyt=9999.
      tzt=8888.
      if(idam.ge.2) then
        if(istune.ge.1) then
!.....CHECK Y TUNE
          dtuney=abs(abs(ty(1))-abs(tuney))
          if(dtuney.gt.etune(2).or.(ty(1)*tuney).lt.0) then
            write(6,*)'Y TUNE DIFFERENT FROM EXPECTED'
            write(6,*)-ty(1),-tuney
            nty=1
            do nt=2,narm
              dtuney=abs(abs(ty(nt))-abs(tuney))
              if(dtuney.le.etune(2).and.(ty(nt)*tuney).gt.0) then
                nty=nt
                goto 8
              endif
            enddo
 8          if(nty.gt.1) then
              write(6,*)'EXPECTED TUNE Y FOUND AT LINE',nty
            elseif(nty.eq.1) then
              write(6,*)'EXPECTED TUNE Y NOT FOUND'
              if(istune.eq.1) write(6,*)'LINE 1 ASSUMED AS TUNE!!!'
            endif
            ty(1)=ty(nty)
            ty(nty)=tya(1)
            zpy=zypes(1)
            zypes(1)=zypes(nty)
            zypes(nty)=zpy
          endif
        endif
 
        if(istune.eq.2) then
          tyt=tuney
        else
          tyt=ty(1)
        endif
        pyt=abs(zypes(1))
        if(abs(pyt).gt.0) then
          pytr=dble(zypes(1))/pyt
          pyti=dimag(zypes(1))/pyt
          if(pyti.eq.zero.and.pytr.eq.zero) then
            fyt=zero
          else
            fyt=atan2(pyti,pytr)
          endif
        else
          pytr=0d0
          pyti=0d0
          fyt=0d0
        endif
        fyt=fyt/pi*1.8d2
        dty=abs(tyt-txt)
        if(dty.le.eps) then
          write(6,*)'TUNEX AND TUNEY ARE TOO CLOSE'
          tyt=9999.
        endif
      else if(idam.lt.2.and.istune.eq.2) then
        tyt=tuney
      endif
      if(idam.eq.3) then
        if(istune.ge.1) then
!.....CHECK Z TUNE
          dtunez=abs(abs(tz(1))-abs(tunez))
          if(dtunez.gt.etune(3).or.(tz(1)*tunez).lt.0) then
            write(6,*)'Y TUNE DIFFERENT FROM EXPECTED'
            write(6,*)-tz(1),-tunez
            ntz=1
            do nt=2,narm
              dtunez=abs(abs(tz(nt))-abs(tunez))
              if(dtunez.le.etune(3).and.(tz(nt)*tunez).gt.0) then
                ntz=nt
                goto 9
              endif
            enddo
 9          if(ntz.gt.1) then
              write(6,*)'EXPECTED TUNE S FOUND AT LINE',ntz
            elseif(ntz.eq.1) then
              write(6,*)'EXPECTED TUNE S NOT FOUND'
              if(istune.eq.1) write(6,*)'LINE 1 ASSUMED AS TUNE!!!'
            endif
            tz(1)=tz(ntz)
            tz(ntz)=tsa(1)
            zpz=zspes(1)
            zspes(1)=zspes(ntz)
            zspes(ntz)=zpz
          endif
        endif
 
        if(istune.eq.2) then
          tzt=tunez
        else
          tzt=tz(1)
        endif
        pzt=abs(zspes(1))
        if(abs(pzt).gt.0) then
          pztr=dble(zspes(1))/pzt
          pzti=dimag(zspes(1))/pzt
          if(pzti.eq.zero.and.pztr.eq.zero) then
            fzt=zero
          else
            fzt=atan2(pzti,pztr)
          endif
        else
          pztr=0d0
          pzti=0d0
          fzt=0d0
        endif
        fzt=fzt/pi*1.8d2
        dtxz=abs(tzt-txt)
        if(dtxz.le.eps) then
          write(6,*)'TUNEX AND TUNES ARE TOO CLOSE'
          tzt=8888.
        endif
        dtyz=abs(tzt-tyt)
        if(dtyz.le.eps) then
          write(6,*)'TUNEY AND TUNES ARE TOO CLOSE'
          tzt=8888.
        endif
      else if(idam.lt.3.and.istune.eq.2) then
        tzt=tunez
      endif
 
!
!.. X SIGNAL PROCESSING
!
 
!       write(30,*)
!       write(30,*)'ANALYSIS OF X SIGNAL, CASE:',iunit
!       write(30,*)'TUNE X = ',-txt
!       write(30,*)
!       write(30,'(3a)')'        Line Frequency      Amplitude',          &
!      &'             Phase       ',                                      &
!      &'      Error         mx  my  ms  p'
!       write(30,*)
 
      iscax=0
      imissx=0
      do n=1,narm
!.....CHECK WITH PREVIOUS HARMONICS
        do k=1,n-1
          dt=abs(tx(k)-tx(n))
          if(dt.le.checkn) then
            iscax=1
          endif
        enddo
!.....SEARCH LINEAR COMBINATIONS
        isearx=0
        ordcx=3000 ! large number needed here
        do l=-nr,nr
          do m=-nr,nr
            do k=-nr,nr
              do j=-nr,nr
                az=l*txt+m*tyt+k*tzt+j
                ex=abs(az-tx(n))
                if(ex.lt.eps) then
                  isearx=isearx+1
!.....check for the lowest possible order of combination
                  ordc=abs(l)+abs(m)+abs(k)
                  if(ordc.lt.ordcx) then
                    l1=l
                    m1=m
                    k1=k
                    j1=j
                    epsx=ex
                    ordcx=ordc

                    ! insertion P.Belanger
                    !--------------------
                    j_vec(n)=l
                    k_vec(n)=m
                    l_vec(n)=k
                    m_vec(n)=j
                    order_vec(n)=ordc
                    !--------------------
                  endif
                endif
              enddo
            enddo
          enddo
        enddo
        if(isearx.ge.1) then
          px=abs(zxpes(n))
          if(abs(px).gt.0) then
            pxtr=dble(zxpes(n))/px
            pxti=dimag(zxpes(n))/px
            if(pxti.eq.zero.and.pxtr.eq.zero) then
              fx=zero
            else
              fx=atan2(pxti,pxtr)
            endif
          else
            pxtr=0d0
            pxti=0d0
            fx=0d0
          endif
          fx=fx/pi*1.8d2

c...  BEGIN INSERTION ROGELIO, to return parameters
!             if(l1.eq.1.and.m1.eq.0.and.k1.eq.0.and.
!      &           j1.eq.0.and.flagad(1).eq.0) then
!               amplitude(1)=px
!               phase(1)=-fx
!               tunexy(1)=-tx(n)
!               flagad(1)=1
!             endif
!             if(l1.eq.-2.and.m1.eq.0.and.k1.eq.0.and.
!      &           flagad(2).eq.0) then
!               amplitude(2)=px
!               phase(2)=-fx
!               flagad(2)=1
!             endif
! 
!             if(flagad(3).eq.0.and.l1.eq.0.and.m1.eq.1.and.
!      &           k1.eq.0)then
!               amplitude(3)=px
!               phase(3)=-fx
!               flagad(3)=1
!             endif
! 
!             if(flagad(7).eq.0.and.l1.eq.-3.and.m1.eq.0.and.
!      &           k1.eq.0)then
!               amplitude(7)=px
!               phase(7)=-fx
!               flagad(7)=1
!             endif
! 
!             if(flagad(8).eq.0.and.l1.eq.-4.and.m1.eq.0.and.
!      &           k1.eq.0)then
!               amplitude(8)=px
!               phase(8)=-fx
!               flagad(8)=1
!             endif
! 
!             if(flagad(9).eq.0.and.l1.eq.-5.and.m1.eq.0.and.
!      &           k1.eq.0.and.j1.eq.0)then
!               amplitude(9)=px
!               phase(9)=-fx
!               flagad(9)=1
!             endif
! 
!             if(flagad(10).eq.0.and.l1.eq.2.and.m1.eq.0.and.
!      &           k1.eq.0)then
!               amplitude(10)=px
!               phase(10)=-fx
!               flagad(10)=1
!             endif
! 
!             if(flagad(11).eq.0.and.l1.eq.0.and.m1.eq.0.and.
!      &           k1.eq.0)then
!               amplitude(11)=px
!               phase(11)=-fx
!               flagad(11)=1
!             endif
! 
!             if(flagad(13).eq.0.and.l1.eq.0.and.m1.eq.2.and.
!      &           k1.eq.0.and.j1.eq.0)then
!               amplitude(13)=px
!               phase(13)=-fx
!               flagad(13)=1
!             endif
c...  END INSERTION ROGELIO
          
!           write(30,100)n,-tx(n),px,-fx,epsx,l1,m1,k1,j1
          ox(n) = -tx(n)
          ax(n) = px
        elseif(isearx.eq.0) then
          imissx=imissx+1
          px=abs(zxpes(n))
          if(abs(px).gt.0) then
            pxtr=dble(zxpes(n))/px
            pxti=dimag(zxpes(n))/px
            if(pxti.eq.zero.and.pxtr.eq.zero) then
              fx=zero
            else
              fx=atan2(pxti,pxtr)
            endif
          else
            pxtr=0d0
            pxti=0d0
            fx=0d0
          endif
          fx=fx/pi*1.8d2
!           write(30,100)n,-tx(n),px,-fx,eps
          ox(n) = -tx(n)
          ax(n) = px
        endif
 10     continue
      enddo
 
!
!.. Y SIGNAL PROCESSING
!
 
      if(idam.lt.2) then
!         write(30,*)
!         write(30,*) 'NO Y SIGNAL'
        if(istune.eq.2) then
!           write(30,*) 'TUNE Y = ',-tyt
!           write(30,*) 'TUNE S = ',-tzt
          goto 51
        else
!           write(30,*)
!           write(30,*)
          goto 50
        endif
      endif
 
!       write(30,*)
!       write(30,*) 'ANALYSIS OF Y SIGNAL, CASE:',iunit
      if(istune.eq.2) then
!         write(30,*) 'TUNE Y = ',-tyt
      else
!         write(30,*) 'TUNE Y = ',-ty(1)
      endif
!       write(30,*)
!       write(30,'(3a)')'        Line Frequency      Amplitude',          &
!      &'             Phase       ',                                      &
!      &'      Error         mx  my  ms  p'
!       write(30,*)
 
      iscay=0
      imissy=0
      do n=1,narm
!.....CHECK WITH PREVIOUS HARMONICS
        do k=1,n-1
          dt=abs(ty(k)-ty(n))
          if(dt.le.checkn) then
            iscay=1
          endif
        enddo
!.....SEARCH LINEAR COMBINATIONS
        iseary=0
        ordcy=3000 ! large number needed here
        do l=-nr,nr
          do m=-nr,nr
            do k=-nr,nr
              do j=-nr,nr
                az=l*txt+m*tyt+k*tzt+j
                ey=abs(az-ty(n))
                if(ey.lt.eps) then
                  iseary=iseary+1
!.....check for the lowest possible order of combination
                  ordc=abs(l)+abs(m)+abs(k)
                  if(ordc.lt.ordcy) then
                    l1=l
                    m1=m
                    k1=k
                    j1=j
                    epsy=ey
                    ordcy=ordc
                  endif
                endif
              enddo
            enddo
          enddo
        enddo
        if(iseary.ge.1) then
          py=abs(zypes(n))
          if(abs(py).gt.0) then
            pytr=dble(zypes(n))/py
            pyti=dimag(zypes(n))/py
            if(pyti.eq.zero.and.pytr.eq.zero) then
              fy=zero
            else
              fy=atan2(pyti,pytr)
            endif
          else
            pytr=0d0
            pyti=0d0
            fy=0d0
          endif
          fy=fy/pi*1.8d2

c...  BEGIN INSERTION ROGELIO, to return parameters
!             if(n.eq.1) then
!               amplitude(4)=py
!               write(*,*)"p3 ", fy, py, ty(n), tunexy(1)
!               phase(4)=-fy
!               tunexy(2)=-ty(n)
!             endif
!             if(l1.eq.0.and.m1.eq.1.and.k1.eq.0.and.
!      &           j1.eq.0.and.flagad(4).eq.0) then
!               amplitude(4)=py
!               phase(4)=-fy
!               tunexy(2)=-ty(n)
! c              write(*,*)"p32 ", fy, py, ty(n)
!               flagad(4)=1
!             endif
!             if(l1.eq.0.and.m1.eq.-2.and.k1.eq.0.and.
!      &           flagad(5).eq.0) then
!               amplitude(5)=py
!               phase(5)=-fy
!               flagad(5)=1
!             endif
! 
!             if(flagad(6).eq.0.and.l1.eq.1.and.m1.eq.0.and.
!      &         k1.eq.0.and.j1.eq.0)then
!               amplitude(6)=py
!               phase(6)=-fy
!               flagad(6)=1
!             endif
! 
!             if(flagad(12).eq.0.and.l1.eq.0.and.m1.eq.-3.and.
!      &         k1.eq.0.and.j1.eq.0)then
!               amplitude(12)=py
!               phase(12)=-fy
!               flagad(12)=1
!             endif
! 
!             if(flagad(14).eq.0.and.l1.eq.-1.and.m1.eq.-1.and.
!      &           k1.eq.0.and.j1.eq.0)then
!               amplitude(14)=py
!               phase(14)=-fy
!               flagad(14)=1
!             endif
c...  END INSERTION ROGELIO

!           write(30,100)n,-ty(n),py,-fy,epsy,l1,m1,k1,j1
          oy(n) = -ty(n)
          ay(n) = py
        elseif(iseary.eq.0) then
          imissy=imissy+1
          py=abs(zypes(n))
          if(abs(py).gt.0) then
            pytr=dble(zypes(n))/py
            pyti=dimag(zypes(n))/py
            if(pyti.eq.zero.and.pytr.eq.zero) then
              fy=zero
            else
              fy=atan2(pyti,pytr)
            endif
          else
            pytr=0d0
            pyti=0d0
            fy=0d0
          endif
          fy=fy/pi*1.8d2
!           write(30,100)n,-ty(n),py,-fy,eps
          oy(n) = -ty(n)
          ay(n) = py
        endif
      enddo
 
!
!.. S SIGNAL PROCESSING
!
 
 51   if(idam.lt.3) then
!         write(30,*)
!         write(30,*) 'NO S SIGNAL'
        if(istune.eq.2) then
!           write(30,*) 'TUNE S = ',-tzt
        else
!           write(30,*)
        endif
!         write(30,*)
        goto 50
      endif
 
!       write(30,*)
!       write(30,*) 'ANALYSIS OF S SIGNAL, CASE:',iunit
      if(istune.eq.2) then
!         write(30,*) 'TUNE S = ',-tzt
      else
!         write(30,*) 'TUNE S = ',-tz(1)
      endif
!       write(30,*)
!       write(30,'(3a)')'        Line Frequency      Amplitude',          &
!      &'             Phase       ',                                      &
!      &'      Error         mx  my  ms  p'
!       write(30,*)
 
      iscaz=0
      imissz=0
      do n=1,narm
!.....CHECK WITH PREVIOUS HARMONICS
        do k=1,n-1
          dt=abs(tz(k)-tz(n))
          if(dt.le.checkn) then
            iscaz=1
          endif
        enddo
!.....SEARCH LINEAR COMBINATIONS
        isearz=0
        ordcz=3000 ! large number needed here
        do l=-nr,nr
          do m=-nr,nr
            do k=-nr,nr
              do j=-nr,nr
                az=l*txt+m*tyt+k*tzt+j
                ez=abs(az-tz(n))
                if(ez.lt.eps) then
                  isearz=isearz+1
!.....check for the lowest possible order of combination
                  ordc=abs(l)+abs(m)+abs(k)
                  if(ordc.lt.ordcz) then
                    l1=l
                    m1=m
                    k1=k
                    j1=j
                    epsz=ez
                    ordcz=ordc
                  endif
                endif
              enddo
            enddo
          enddo
        enddo
        if(isearz.ge.1) then
          pz=abs(zspes(n))
          if(abs(pz).gt.0) then
            pztr=dble(zspes(n))/pz
            pzti=dimag(zspes(n))/pz
            if(pzti.eq.zero.and.pztr.eq.zero) then
              fz=zero
            else
              fz=atan2(pzti,pztr)
            endif
          else
            pztr=0d0
            pzti=0d0
            fz=0d0
          endif
          fz=fz/pi*1.8d2
!           write(30,100)n,-tz(n),pz,-fz,epsz,l1,m1,k1,j1
          os(n) = -tz(n)
          as(n) = pz
        elseif(isearz.eq.0) then
          imissz=imissz+1
          pz=abs(zspes(n))
          if(abs(pz).gt.0) then
            pztr=dble(zspes(n))/pz
            pzti=dimag(zspes(n))/pz
            if(pzti.eq.zero.and.pztr.eq.zero) then
              fz=zero
            else
              fz=atan2(pzti,pztr)
            endif
          else
            pztr=0d0
            pzti=0d0
            fz=0d0
          endif
          fz=fz/pi*1.8d2
!           write(30,100)n,-tz(n),pz,-fz,eps
          os(n) = -tz(n)
          as(n) = pz
        endif
 30     continue
      enddo

 50     continue
! 50    write(6,*)'================================================'
!       write(6,*)'HARMONICS NON-IDENTIFIED IN X = ',imissx
!       write(6,*)'HARMONICS NON-IDENTIFIED IN Y = ',imissy
!       write(6,*)'HARMONICS NON-IDENTIFIED IN S = ',imissz
      narm2=narm/2d0
      imiss=imissx+imissy+imissz
      if(imiss.ge.narm) then
!         write(6,*)'WARNING: CHECK EPS'
      endif
      isca=iscax+iscay+iscaz
      if(isca.ge.1) then
!         write(6,*)'WARNING: TOO CLOSE BY HARMONICS DETECTED'
!         write(6,*)'WARNING: TRY A LARGER NUMBER OF TURNS'
      endif
!       write(6,*)'WROTE TO lines THE SPECTRAL ANALYSIS OF CASE:',iunit
!       write(6,*)'============================================='
 
100   format(i3,1x,e19.13,2(1x,e18.12),1x,e17.11,:,4(1x,i3))
 
      return
      end
      subroutine datspe(iunit,idam,ir,nt1,nt2,nturn,imeth,narm,iana,
     & x, xp, y, yp, s, sp, n_points)
!=======================================================================
!
!  SUBROUTINE DATSPE
!
!  THIS PROGRAM CALCULATES THE SPECTRUM OF A SINGLE FILE
!  OF TRACKING DATA WRITTEN IN A STANDARD ASCII FROM IN 2*IDAM COLUMNS.
!
!  WITH THE DATA SOME INFORMATIONS ARE REQUESTED:
!
!    IUNIT = FORTRAN UNIT OF INPUT DATA
!    IDAM  = DIMENSION OF PHASE SPACE
!    IR    = FLAG FOR REAL SIGNAL (1=REAL)
!    NT1   = INITIAL TURN TO BE ANALYZED
!    NT2   = FINAL TURN TO BE ANALYZED
!
!  THE TUNE AND THE LINES ARE CALCULATED WITH THE ROUTINE SPECTRUM
!
!    NARM : THE NUMBER OF HARMONIC TO BE CALCULATED
!    IMETH: THE CHIOCE ON THE WINDOWING
!           1 HANNING WINDOW (TUNENEWT)
!           2 NO WINDOW (TUNELASR)
!
!  N.B.: ONLY ONE FILE IS ANALIZED AND THE OUTPUT IS PLACED IN
!        THE COMMONS TUNE AND FCOE TOGETHER WITH THE TRACKING DATA.
!  N.B.: NTURN IS AN OUTPUT PARAMETER!
!
!  AUTHOR: R.BARTOLINI 21/08/1996
!  MODIFIED 17/09/1996: ADDED THE IR OPTION
!  MODIFIED 30/05/1998: ADDED FFT
!
!=======================================================================
! ==============================================================================
! 2012 KEVIN: CHANGES MADE:
! ADDED INPUT VARIABLES: N_POINTS, X, XP, Y, YP, S, SZ
! COMMENTED MOST WRITE(6,*) (TO SCREEN)
! ==============================================================================


      implicit none
      integer iana,idam,imeth,ir,iunit,j,k,narm,nt1,nt2,nturn,nturn2
      complex zsing
      integer mterm
      parameter(mterm=300)
      integer maxiter
      parameter(maxiter=100000)
!       double precision s,sp,x,xp,y,yp
!       common/data/x(maxiter),y(maxiter),xp(maxiter),yp(maxiter),
!      &s(maxiter),sp(maxiter)
      double precision tsa,txa,tya
      common/tune/txa(mterm),tya(mterm),tsa(mterm)
      double complex zspes,zxpes,zypes
      common/fcoe/zxpes(mterm),zypes(mterm),zspes(mterm)
      dimension zsing(maxiter)


      integer n_points
      double precision x(n_points), xp(n_points),
     & y(n_points), yp(n_points), s(n_points), sp(n_points)


!.....READ INPUT FILE
!  
!       write(6,*)'****************************'
!       write(6,*)'ANALYZING UNIT',iunit
!       write(6,*)'****************************'
!       if(idam.eq.1) then
!         do j=1,maxiter
!           read(iunit,*,end=990)x(j),xp(j)
!         enddo
!       elseif(idam.eq.2) then
!         do j=1,maxiter
!           read(iunit,*,end=990)x(j),xp(j),y(j),yp(j)
!         enddo
!       elseif(idam.eq.3) then
!         do j=1,maxiter
!           read(iunit,*,end=990)x(j),xp(j),y(j),yp(j),s(j),sp(j)
!         enddo
!       endif
! 990   nturn=j-1       ! check this it is always larger by 1
!       write(6,*)'NUMBER OF TURNS DETECTED IN THE INPUT',nturn
!       if(nturn.eq.0) then
!         write(6,*)'No data detected form bpm data files'
!         write(6,*)'Check filenames program stops!'
!         stop
!       endif
!       rewind(iunit)
! ==============================================================================
! FOR NoO
! ==============================================================================
      nturn = n_points
c      write(6,*)'NUMBER OF TURNS DETECTED IN THE INPUT',nturn

!.....CHECK FOR REQUIRED SECTIONING OF THE SIGNAL
      nturn2=nt2-nt1+1
      if(nturn.gt.nturn2) then
        if(idam.eq.1) then
          do k=1,nturn2
            x(k)=x(nt1+k-1)
            xp(k)=xp(nt1+k-1)
          enddo
        elseif(idam.eq.2) then
          do k=1,nturn2
            x(k)=x(nt1+k-1)
            xp(k)=xp(nt1+k-1)
            y(k)=y(nt1+k-1)
            yp(k)=yp(nt1+k-1)
          enddo
        elseif(idam.eq.3) then
          do k=1,nturn2
            x(k)=x(nt1+k-1)
            xp(k)=xp(nt1+k-1)
            y(k)=y(nt1+k-1)
            yp(k)=yp(nt1+k-1)
            s(k)=s(nt1+k-1)
            sp(k)=sp(nt1+k-1)
          enddo
        endif
        nturn=nturn2
        write(6,*)'REDUCTION OF THE SIGNAL PERFORMED.'
        write(6,*)'NEW NUMBER OF TURNS ANALYZED =',nturn
        write(6,*)'INTERVAL',nt1,nt2
      endif
 
!.....FLAG FOR THE ANALYSIS OF THE REAL PART ONLY
 
      if(ir.eq.1) then
        do j=1,nturn
          xp(j)=0d0
          yp(j)=0d0
          sp(j)=0d0
        enddo
!         write(6,*)'WARNING: ONLY THE REAL PART OF THE SIGNAL IS'
!         write(6,*)'         CONSIDERED FOR SPECTRUM CALCULATION'
      endif
 
!.....SPECTRUM CALCULATION
 
      if(idam.eq.1) then
        call spectrum(x,xp,nturn,txa,zxpes,narm,imeth)
c        write(6,*)'TURNS ANALYZED AND TUNE'
c        write(6,*)nturn,-txa(1)
        if(iana.eq.2) then
          call fftr(x,xp,nturn,zsing,imeth)
        endif
      else if(idam.eq.2) then
        call spectrum(x,xp,nturn,txa,zxpes,narm,imeth)
        call spectrum(y,yp,nturn,tya,zypes,narm,imeth)
c        write(6,*)'TURNS ANALYZED AND TUNE'
c        write(6,*)nturn,-txa(1),-tya(1)
        if(iana.eq.2) then
          call fftr(x,xp,nturn,zsing,imeth)
          call fftr(y,yp,nturn,zsing,imeth)
        endif
      elseif(idam.eq.3) then
        call spectrum(x,xp,nturn,txa,zxpes,narm,imeth)
        call spectrum(y,yp,nturn,tya,zypes,narm,imeth)
        call spectrum(s,sp,nturn,tsa,zspes,narm,imeth)
c        write(6,*)'TURNS ANALYZED AND TUNE'
c        write(6,*)nturn,-txa(1),-tya(1),-tsa(1)
        if(iana.eq.2) then
          call fftr(x,xp,nturn,zsing,imeth)
          call fftr(y,yp,nturn,zsing,imeth)
          call fftr(s,sp,nturn,zsing,imeth)
        endif
      endif
 
      end
      subroutine fftr(x,xp,maxn,zsing,meth)
!=======================================================================
!
! SUBROUTINE FFTR
!
! COMPUTES THE FFT.
! X, XP ARE THE COORDINATES OF THE ORBIT AND MAXN IS THE
! LENGTH IF THE ORBIT.
!
! METH SELECTS THE HANNING WINDOW (1) OR NOT (2)
!
! FOURIER COEFFICIENTS ARE GIVEN IN ZSING
! AND THE SPECTRUM IS WRITTEN IN FORT.31
!
! AUTHOR:     R. BARTOLINI
!
!=======================================================================
      implicit none
      integer maxn,maxn2,meth,mf,mft,nf,npoint
      double precision amp,omnf,pha,step,x,xp
      double complex z
      complex zsing
      integer maxiter
      parameter(maxiter=100000)
      integer mbpm
      double precision pieni,zero,one,two,pi,duepi,piph
      parameter(mbpm=1000,pieni=1d-17,zero=0d0,one=1d0,two=2d0)
      dimension x(maxiter),xp(maxiter),zsing(maxiter)
      dimension z(maxiter)
 
! INITIALIZATION
      duepi=8d0*atan(1d0)
 
!.............................................................
!    ESTIMATION OF TUNE WITH FFT
!.............................................................
      mft=int(log(dble(maxn))/log(2d0))
      npoint=2**mft
      maxn2=maxn/2
      step=duepi/maxn
      if (meth.eq.1) then
        do mf=1,maxn
          z(mf)=dcmplx(x(mf),xp(mf))*(1d0+cos(step*(mf-maxn2)))
          zsing(mf)=z(mf)
        enddo
      elseif (meth.eq.2) then
        do mf=1,maxn
          z(mf)=dcmplx(x(mf),xp(mf))
          zsing(mf)=z(mf)
        enddo
      endif
 
      call cfft(zsing,-mft)
 
      do nf=1,npoint
        omnf=(nf-1d0)/npoint
        amp=abs(zsing(nf))
        if(aimag(zsing(nf)).eq.zero.and.real(zsing(nf)).eq.zero) then
          pha=zero
        else
          pha=atan2(aimag(zsing(nf)),real(zsing(nf)))/duepi*360
        endif
        write(31,*)omnf,amp/npoint,pha
      enddo
 
!............................................................
      return
!............................................................
      end
      subroutine readres(imin,imax,lr,mr,kr,narm,idam,idamx,iout)
!=======================================================================
!
!  SUBROUTINE READRES
!
!  READ THE HARMONICS FROM THE FILE lines
!  AND WRITE SEPARATELY THE SELECTED ONE
!
!  INPUT PARAMETERS:
!
!  IMIN,IMAX ENUMERATE THE DIFFERENT CASES STORED IN FILE lines
!  NARM : NUMERO DI ARMONICHE DA ORDINARE
!  IR   : FLAG FOR COMPLEX OR REAL SIGNALS:
!         1 = INPUT SIGNAL IS REAL.
!
!  N.B. :  THE STRENGTH OF NEXT TO LEADING LINES ARE NORMALIZED AT THE
!          TUNE STRENGTH.
!  MIND :  THE ARRAYS TX,TY,TZ ARE USED IN ORDER NOT TO CHANGE
!          THE ARRAYS TUNEX,TUNEY,TUNEZ.
!
!  IOU  : FORTRAN UNIT OF THE OUTPUT, IF 0 NO OUTPUT
!
!  AUTHOR: R.BARTOLINI 21/08/1996
!  LAST MODIFIED 08/10/1997
!
!=======================================================================
      implicit none
      integer i,idam,idamx,ifoun,imax,imin,iout,istep,j,k,ki,kr,l,lr,m, &
     &mr,n,narm,ni
      double precision ex,fx,px,tx
      integer mterm
      parameter(mterm=300)
      integer mbpm
      double precision pieni,zero,one,two,pi,duepi,piph
      parameter(mbpm=1000,pieni=1d-17,zero=0d0,one=1d0,two=2d0)
      double precision dtfx,dtpx,dttx
      common/dt/dtpx(mbpm),dtfx(mbpm),dttx(mbpm)
      dimension tx(mterm)
      character*200 ch,ch1
 
! INITIALIZATION
      ch1=' '
 
      if(imin.le.imax) then
        istep=1
      else
        istep=-1
      endif
 
      rewind(30)
      do i=imin,imax,istep
        write(6,*)'ANALYZING FILE lines FOR LINE IDENTIFICATION CASE:',i
        do ki=1,idam
 1        read(30,'(A)',end=2) ch
! Find start for each plane
          if(index(ch(1:10),'ANALYSIS').gt.0d0) then
            read(30,'()',end=2)
            read(30,'()',end=2)
            read(30,'()',end=2)
            read(30,'()',end=2)
            ifoun=0
            do ni=1,narm
              read(30,'(A)',end=2) ch
              ch(198:200)=' / '
! Detect end of data and exit to 3 to study next plane
              if(ch.eq.ch1) goto 3
              l=-9999
              m=-9999
              k=-9999
              j=-9999
              read(ch,*)n,tx(n),px,fx,ex,l,m,k,j
              if(ki.eq.idamx) then
                if((l.eq.lr.and.m.eq.mr.and.k.eq.kr).and.               &
     &ifoun.eq.0) then
                  if(iout.ne.0) then
                    write(iout,101)tx(1),tx(ni),px,fx,ex,l,m,k,j,       &
     &' CASE:',i
                  endif
! Fills the vector for eventual DT calculation
                  dtpx(i)=px
                  dtfx(i)=fx
                  dttx(i)=tx(ni)
                  ifoun=1
                endif
              endif
            enddo
            if(ki.eq.idamx.and.ifoun.eq.0.and.iout.ne.0) then
              write(iout,101)0d0,0d0,0d0,0d0,0d0,lr,mr,kr,0,' CASE:',i
            endif
          else
            goto 1
          endif
        enddo
 3      continue
      enddo
 2    continue
      rewind(30)
 
101   format(4(1x,e18.12),1x,e17.12,:,4(1x,i3),a6,1x,i3)
 
      return
      end
      subroutine readdt3(imin,imax,narm,idam,iout)
!=======================================================================
!
!  SUBROUTINE READDT3
!
!  READ THE HARMONICS FROM THE FILE lines
!  AND CONVERT THEM INTO DRIVING TERMS (SEE LHC Project Report 132)
!
!  INPUT PARAMETERS:
!
!  IMIN,IMAX ENUMERATE THE DIFFERENT CASES STORED IN FILE lines
!  NARM : NUMERO DI ARMONICHE DA ORDINARE
!  IR   : FLAG FOR COMPLEX OR REAL SIGNALS:
!         1 = INPUT SIGNAL IS REAL.
!
!
!  IOUT  : FORTRAN UNIT OF THE OUTPUT
!
!  AUTHOR: R.BARTOLINI 21/08/1996
!
!=======================================================================
      implicit none
      integer i,idam,imax,imin,iout,istep,narm
      double precision a0012,a0030,a1002,a1011,a1020,a1110,a1200,a2001, &
     &a2010,a3000,aaux,c0012,c0030,c1002,c1011,c1020,c1110,c1200,c2001, &
     &c2010,c3000,f0012,f0030,f1002,f1011,f1020,f1110,f1200,f2001,f2010,&
     &f3000,faux,fx,fy,px,py,s0012,s0030,s1002,s1011,s1020,s1110,s1200, &
     &s2001,s2010,s3000,tx,ty,y1,y2
      double complex z1011,z2100,zfli,zline
      integer mbpm
      double precision pieni,zero,one,two,pi,duepi,piph
      parameter(mbpm=1000,pieni=1d-17,zero=0d0,one=1d0,two=2d0)
      double precision dtfx,dtpx,dttx
      common/dt/dtpx(mbpm),dtfx(mbpm),dttx(mbpm)
      dimension tx(mbpm),px(mbpm),fx(mbpm),ty(mbpm),py(mbpm),fy(mbpm),  &
     &aaux(mbpm),faux(mbpm)
 
! INITIALIZATION
      pi=4d0*atan(1d0)
      piph=pi/180d0
      if(imin.le.imax) then
        istep=1
      else
        istep=-1
      endif
 
      do i=imin,imax,istep
        write(iout,*)
        write(iout,*)'THIRD ORDER CONJUGATING FUNCTION COEFFICIENTS',   &
     &', CASE: ',i
        write(iout,*)
        write(iout,'(2a)')'jklm          Amplitude              Phase', &
     &'           Cos             Sin'
        write(iout,*)
! tune lines
        dtpx(i)=0d0
        dtfx(i)=0d0
        dttx(i)=0d0
        call readres(imin,i,1,0,0,narm,idam,1,0)
        tx(i)=dttx(i)
        px(i)=dtpx(i)
        fx(i)=dtfx(i)
        if(idam.ge.2) then
 
!
! attenzione se analizzo il moto verticale ci
! va la y. cosi come e' faccio un overwrite su dtpx,dtfx,sttx
!
! inoltre la linea eccitata in orizzontale
! e' 1-j+k,m-l e non l-m
! idem per il vertical avro' k-j e non j-k
! correggere.
!
 
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,0,1,0,narm,idam,2,0)
          ty(i)=dttx(i)
          py(i)=dtpx(i)
          fy(i)=dtfx(i)
        endif
 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! One dimensional resonances horizontal C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
! -2,0 line and h3000
        dtpx(i)=0d0
        dtfx(i)=0d0
        dttx(i)=0d0
        call readres(imin,i,-2,0,0,narm,idam,1,0)
        a3000=dtpx(i)/2d0/3d0/px(i)**2
        f3000=dtfx(i)-fx(i)+90d0
        c3000=a3000*cos(f3000*piph)
        s3000=a3000*sin(f3000*piph)
! 2,0 line and h1200
        dtpx(i)=0d0
        dtfx(i)=0d0
        dttx(i)=0d0
        call readres(imin,i,2,0,0,narm,idam,1,0)
        a1200=dtpx(i)/2d0/px(i)**2
        f1200=dtfx(i)-fx(i)+90d0
        c1200=a1200*cos(f1200*piph)
        s1200=a1200*sin(f1200*piph)
! stores this data for eventual subresonance calculation
        aaux(i)=a1200
        faux(i)=f1200
 
! other lines in case of 4D motion
        if(idam.gt.1) then
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! One dimensional vertical resonances (vertical motion needed) C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! 0,-2 line and h0030
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,0,-2,0,narm,idam,2,0)
          a0030=dtpx(i)/2d0/3d0/py(i)**2
          f0030=dtfx(i)-fy(i)+90d0
          c0030=a0030*cos(f0030*piph)
          s0030=a0030*sin(f0030*piph)
! 0,2 line and h0012
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,0,2,0,narm,idam,2,0)
          a0012=dtpx(i)/2d0/py(i)**2
          f0012=dtfx(i)-fy(i)+90d0
          c0012=a0012*cos(f0012*piph)
          s0012=a0012*sin(f0012*piph)
 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Normal resonances from horizontal motion C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
! 0,-2 line and h1020
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,0,-2,0,narm,idam,1,0)
          a1020=dtpx(i)/2d0/py(i)**2
          f1020=dtfx(i)-fx(i)+90d0
          c1020=a1020*cos(f1020*piph)
          s1020=a1020*sin(f1020*piph)
! 0,2 line and h1002
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,0,2,0,narm,idam,1,0)
          a1002=dtpx(i)/2d0/py(i)**2
          f1002=dtfx(i)-fx(i)+90d0
          c1002=a1002*cos(f1002*piph)
          s1002=a1002*sin(f1002*piph)
! 0,0 line and h1011 --subresonance: h2100 needed---
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,0,0,0,narm,idam,1,0)
          zline=dtpx(i)*exp(dcmplx(0d0,dtfx(i)*piph))
!.......rebulding the contribution to the (0,0) from h2100
          zfli=4*px(i)**2*dcmplx(0d0,-1d0)
          z2100=aaux(i)*exp(dcmplx(0d0,-faux(i)*piph))*zfli
!.......subtracting the contribution from h2100 to the (0,0)
          z1011=zline-z2100
!.......building h1011 from z1011
          a1011=abs(z1011)/2d0/py(i)**2
          y1=dble(z1011)
          y2=dimag(z1011)
          if(y1.eq.zero.and.y2.eq.zero) then
            f1011=zero
          else
            f1011=atan2(y2,y1)*180/pi-fx(i)+90d0
          endif
          c1011=a1011*cos(f1011*piph)
          s1011=a1011*sin(f1011*piph)
 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Skew resonances from horizontal motion C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! -1,-1 line and h2010
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,-1,-1,0,narm,idam,1,0)
          a2010=dtpx(i)/2d0/2d0/px(i)/py(i)
          f2010=dtfx(i)-fx(i)+90d0
          c2010=a2010*cos(f2010*piph)
          s2010=a2010*sin(f2010*piph)
!  1,-1 line and h1110
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,1,-1,0,narm,idam,1,0)
          a1110=dtpx(i)/2d0/px(i)/py(i)
          f1110=dtfx(i)-fx(i)+90d0
          c1110=a1110*cos(f1110*piph)
          s1110=a1110*sin(f1110*piph)
! -1,1 line and h2001
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,-1,1,0,narm,idam,1,0)
          a2001=dtpx(i)/2d0/2d0/px(i)/py(i)
          f2001=dtfx(i)-fx(i)+90d0
          c2001=a2001*cos(f2001*piph)
          s2001=a2001*sin(f2001*piph)
 
        endif
        write(iout,100)'3000',a3000,f3000,c3000,s3000
        write(iout,200)'2100',a1200,-f1200
        write(iout,100)'1200',a1200,f1200,c1200,s1200
        write(iout,200)'0300',a3000,-f3000
        write(iout,*)
        write(iout,100)'1020',a1020,f1020,c1020,s1020
        write(iout,100)'1011',a1011,f1011,c1011,s1011
        write(iout,100)'1002',a1002,f1002,c1002,s1002
        write(iout,200)'0120',a1002,-f1002
        write(iout,200)'0111',a1011,-f1011
        write(iout,200)'0102',a1020,-f1020
        write(iout,*)
        write(iout,100)'0030',a0030,f0030,c0030,s0030
        write(iout,200)'0021',a0012,-f0012
        write(iout,100)'0012',a0012,f0012,c0012,s0012
        write(iout,200)'0003',a0030,-f0030
        write(iout,*)
        write(iout,100)'2010',a2010,f2010,c2010,s2010
        write(iout,100)'1110',a1110,f1110,c1110,s1110
        write(iout,200)'0210',a2001,-f2001
        write(iout,100)'2001',a2001,f2001,c2001,s2001
        write(iout,200)'1101',a1110,-f1110
        write(iout,200)'0201',a2010,-f2010
 
      enddo
 
 100  format(1x,a4,4(1x,e19.13))
 200  format(1x,a4,2(1x,e19.13))
      return
      end
      subroutine readdt4(imin,imax,narm,idam,iout)
!=======================================================================
!
!  SUBROUTINE READDT4
!
!  READ THE HARMONICS FROM THE FILE lines
!  AND CONVERT THEM INTO DRIVING TERMS (SEE LHC Project Report 132)
!
!  IT IS ASSUMED THAT THE DATA ARE ALREADY CLEANED FROM THIRD ORDER
!
!  INPUT PARAMETERS:
!
!  IMIN,IMAX ENUMERATE THE DIFFERENT CASES STORED IN FILE lines
!  NARM : NUMERO DI ARMONICHE DA ORDINARE
!  IR   : FLAG FOR COMPLEX OR REAL SIGNALS:
!         1 = INPUT SIGNAL IS REAL.
!
!
!  IOUT  : FORTRAN UNIT OF THE OUTPUT
!
!  AUTHOR: R.BARTOLINI 21/08/1996
!
!=======================================================================
      implicit none
      integer i,idam,imax,imin,iout,istep,narm
      double precision a0013,a0040,a1003,a1012,a1021,a1030,a1120,a1201, &
     &a1210,a1300,a2002,a2011,a2020,a3001,a3010,a4000,aaux,aauy,c0013,  &
     &c0040,c1003,c1012,c1021,c1030,c1120,c1201,c1210,c1300,c2002,c2011,&
     &c2020,c3001,c3010,c4000,f0013,f0040,f1003,f1012,f1021,f1030,f1120,&
     &f1201,f1210,f1300,f2002,f2011,f2020,f3001,f3010,f4000,faux,fauy,  &
     &fx,fy,px,py,s0013,s0040,s1003,s1012,s1021,s1030,s1120,s1201,s1210,&
     &s1300,s2002,s2011,s2020,s3001,s3010,s4000,tx,ty,y1,y2
      double complex z1012,z1021,z2011,z2101,z2110,z3100,zfli,zline
      integer mbpm
      double precision pieni,zero,one,two,pi,duepi,piph
      parameter(mbpm=1000,pieni=1d-17,zero=0d0,one=1d0,two=2d0)
      double precision dtfx,dtpx,dttx
      common/dt/dtpx(mbpm),dtfx(mbpm),dttx(mbpm)
      dimension tx(mbpm),px(mbpm),fx(mbpm),ty(mbpm),py(mbpm),fy(mbpm),  &
     &aaux(mbpm),faux(mbpm),aauy(mbpm),fauy(mbpm)
 
! INITIALIZATION
      pi=4d0*atan(1d0)
      piph=pi/180d0
      if(imin.le.imax) then
        istep=1
      else
        istep=-1
      endif
 
      do i=imin,imax,istep
        write(iout,*)
        write(iout,*)'FOURTH ORDER CONJUGATING FUNCTION COEFFICIENTS',  &
     &', CASE: ',i
        write(iout,*)
        write(iout,'(2a)')'jklm          Amplitude              Phase', &
     &'           Cos             Sin'
        write(iout,*)
! tune lines
        dtpx(i)=0d0
        dtfx(i)=0d0
        dttx(i)=0d0
        call readres(imin,i,1,0,0,narm,idam,1,0)
        tx(i)=dttx(i)
        px(i)=dtpx(i)
        fx(i)=dtfx(i)
        if(idam.ge.2) then
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,0,1,0,narm,idam,2,0)
          ty(i)=dttx(i)
          py(i)=dtpx(i)
          fy(i)=dtfx(i)
        endif
 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! One dimensional resonances horizontal C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
! -3,0 line and h4000
        dtpx(i)=0d0
        dtfx(i)=0d0
        dttx(i)=0d0
        call readres(imin,i,-3,0,0,narm,idam,1,0)
        a4000=dtpx(i)/2d0/4d0/px(i)**3
        f4000=dtfx(i)-fx(i)+90d0
        c4000=a4000*cos(f4000*piph)
        s4000=a4000*sin(f4000*piph)
! 3,0 line and h1300
        dtpx(i)=0d0
        dtfx(i)=0d0
        dttx(i)=0d0
        call readres(imin,i,3,0,0,narm,idam,1,0)
        a1300=dtpx(i)/2d0/px(i)**3
        f1300=dtfx(i)-fx(i)+90d0
        c1300=a1300*cos(f1300*piph)
        s1300=a1300*sin(f1300*piph)
! stores this data for eventual subresonance calculation
        aaux(i)=a1300
        faux(i)=f1300
 
! other lines in case of 4D motion
        if(idam.gt.1) then
 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! One dimensional vertical resonances (vertical motion needed) C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! 0,-3 line and h0040
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,0,-3,0,narm,idam,2,0)
          a0040=dtpx(i)/2d0/4d0/py(i)**3
          f0040=dtfx(i)-fy(i)+90d0
          c0040=a0040*cos(f0040*piph)
          s0040=a0040*sin(f0040*piph)
! 0,3 line and h0013
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,0,3,0,narm,idam,2,0)
          a0013=dtpx(i)/2d0/py(i)**3
          f0013=dtfx(i)-fy(i)+90d0
          c0013=a0013*cos(f0013*piph)
          s0013=a0013*sin(f0013*piph)
! stores this data for eventual subresonance calculation
          aauy(i)=a0013
          fauy(i)=f0013
 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Normal resonances from horizontal motion C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
! -1,-2 line and h2020
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,-1,-2,0,narm,idam,1,0)
          a2020=dtpx(i)/2d0/2d0/px(i)/py(i)**2
          f2020=dtfx(i)-fx(i)+90d0
          c2020=a2020*cos(f2020*piph)
          s2020=a2020*sin(f2020*piph)
! -1,2 line and h2002
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,-1,2,0,narm,idam,1,0)
          a2002=dtpx(i)/2d0/2d0/px(i)/py(i)**2
          f2002=dtfx(i)-fx(i)+90d0
          c2002=a2002*cos(f2002*piph)
          s2002=a2002*sin(f2002*piph)
! 1,-2 line and h1120
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,1,-2,0,narm,idam,1,0)
          a1120=dtpx(i)/2d0/px(i)/py(i)**2
          f1120=dtfx(i)-fx(i)+90d0
          c1120=a1120*cos(f1120*piph)
          s1120=a1120*sin(f1120*piph)
! -1,0 line and h2011 --subresonance: h3100 needed---
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,-1,0,0,narm,idam,1,0)
          zline=dtpx(i)*exp(dcmplx(0d0,dtfx(i)*piph))
!.......rebulding the contribution to the (-1,0) from h3100
          zfli=2d0*3d0*px(i)**3*dcmplx(0d0,-1d0)
          z3100=aaux(i)*exp(dcmplx(0d0,-faux(i)*piph))*zfli
!.......subtracting the contribution from h3100 to the (-1,0)
          z2011=zline-z3100
!.......building h2011 from z2011
          a2011=abs(z2011)/2d0/2d0/px(i)/py(i)**2
          y1=dble(z2011)
          y2=dimag(z2011)
          if(y1.eq.zero.and.y2.eq.zero) then
            f2011=zero
          else
            f2011=atan2(y2,y1)*180/pi-fx(i)+90d0
          endif
          c2011=a2011*cos(f2011*piph)
          s2011=a2011*sin(f2011*piph)
 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Skew resonances from horizontal motion C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
! -2,-1 line and h3010
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,-2,-1,0,narm,idam,1,0)
          a3010=dtpx(i)/2d0/3d0/px(i)**2/py(i)
          f3010=dtfx(i)-fx(i)+90d0
          c3010=a3010*cos(f3010*piph)
          s3010=a3010*sin(f3010*piph)
!  2,1 line and h1201
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,2,1,0,narm,idam,1,0)
          a1201=dtpx(i)/2d0/px(i)**2/py(i)
          f1201=dtfx(i)-fx(i)+90d0
          c1201=a1201*cos(f1201*piph)
          s1201=a1201*sin(f1201*piph)
!  2,-1 line and h1210
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,2,-1,0,narm,idam,1,0)
          a1210=dtpx(i)/2d0/px(i)**2/py(i)
          f1210=dtfx(i)-fx(i)+90d0
          c1210=a1210*cos(f1210*piph)
          s1210=a1210*sin(f1210*piph)
! -2,1 line and h3001
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,-2,1,0,narm,idam,1,0)
          a3001=dtpx(i)/2d0/3d0/px(i)**2/py(i)
          f3001=dtfx(i)-fx(i)+90d0
          c3001=a3001*cos(f3001*piph)
          s3001=a3001*sin(f3001*piph)
! 0,-3 line and h1030
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,0,-3,0,narm,idam,1,0)
          a1030=dtpx(i)/2d0/py(i)**3
          f1030=dtfx(i)-fx(i)+90d0
          c1030=a1030*cos(f1030*piph)
          s1030=a1030*sin(f1030*piph)
! 0,-1 line and h1021 --subresonance h2110 needed--
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,0,-1,0,narm,idam,1,0)
          zline=dtpx(i)*exp(dcmplx(0d0,dtfx(i)*piph))
!.......rebulding the contribution to the (0,-1) from h2110
          zfli=2d0*2d0*px(i)**2*py(i)*dcmplx(0d0,-1d0)
          z2110=a1201*exp(dcmplx(0d0,-f1201*piph))*zfli
!.......subtracting the contribution from h2110 to the (0,-1)
          z1021=zline-z2110
!.......building h1021 from z1021
          a1021=abs(z1021)/2d0/py(i)**3
          y1=dble(z1021)
          y2=dimag(z1021)
          if(y1.eq.zero.and.y2.eq.zero) then
            f1021=zero
          else
            f1021=atan2(y2,y1)*180/pi-fx(i)+90d0
          endif
          c1021=a1021*cos(f1021*piph)
          s1021=a1021*sin(f1021*piph)
! 0,1 line and h1012 --subresonance h2101 needed--
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,0,1,0,narm,idam,1,0)
          zline=dtpx(i)*exp(dcmplx(0d0,dtfx(i)*piph))
!.......rebulding the contribution to the (0,1) from h2101
          zfli=2d0*2d0*px(i)**2*py(i)*dcmplx(0d0,-1d0)
          z2101=a1210*exp(dcmplx(0d0,-f1210*piph))*zfli
!.......subtracting the contribution from h2101 to the (0,1)
          z1012=zline-z2101
!.......building h1012 from z1012
          a1012=abs(z1012)/2d0/py(i)**3
          y1=dble(z1012)
          y2=dimag(z1012)
          if(y1.eq.zero.and.y2.eq.zero) then
            f1012=zero
          else
            f1012=atan2(y2,y1)*180/pi-fx(i)+90d0
          endif
          c1012=a1012*cos(f1012*piph)
          s1012=a1012*sin(f1012*piph)
! 0,3 line and h1003
          dtpx(i)=0d0
          dtfx(i)=0d0
          dttx(i)=0d0
          call readres(imin,i,0,3,0,narm,idam,1,0)
          a1003=dtpx(i)/2d0/py(i)**3
          f1003=dtfx(i)-fx(i)+90d0
          c1003=a1003*cos(f1003*piph)
          s1003=a1003*sin(f1003*piph)
 
        endif
        write(iout,100)'4000',a4000,f4000,c4000,s4000
        write(iout,200)'3100',a1300,-f1300
        write(iout,100)'1300',a1300,f1300,c1300,s1300
        write(iout,200)'0400',a4000,-f4000
        write(iout,*)
        write(iout,100)'2020',a2020,f2020,c2020,s2020
        write(iout,100)'2011',a2011,f2011,c2011,s2011
        write(iout,100)'2002',a2002,f2002,c2002,s2002
        write(iout,100)'1120',a1120,f1120,c1120,s1120
        write(iout,200)'1102',a1120,-f1120
        write(iout,200)'0220',a2002,-f2002
        write(iout,200)'0211',a2011,-f2011
        write(iout,200)'0202',a2020,-f2020
        write(iout,*)
        write(iout,100)'0040',a0040,f0040,c0040,s0040
        write(iout,200)'0031',a0013,-f0013
        write(iout,100)'0013',a0013,f0013,c0013,s0013
        write(iout,200)'0004',a0040,-f0040
        write(iout,*)
        write(iout,100)'3010',a3010,f3010,c3010,s3010
        write(iout,200)'2110',a1201,-f1201
        write(iout,100)'1210',a1210,f1210,c1210,s1210
        write(iout,200)'0310',a3001,-f3001
        write(iout,100)'3001',a3001,f3001,c3001,s3001
        write(iout,200)'2101',a1210,-f1210
        write(iout,100)'1201',a1201,f1201,c1201,s1201
        write(iout,200)'0301',a3010,-f3010
        write(iout,*)
        write(iout,100)'1030',a1030,f1030,c1030,s1030
        write(iout,100)'1021',a1021,f1021,c1021,s1021
        write(iout,100)'1012',a1012,f1012,c1012,s1012
        write(iout,100)'1003',a1003,f1003,c1003,s1003
        write(iout,200)'0103',a1030,-f1030
        write(iout,200)'0112',a1021,-f1021
        write(iout,200)'0121',a1012,-f1012
        write(iout,200)'0130',a1003,-f1003
      enddo
 
 100  format(1x,a4,4(1x,e19.13))
 200  format(1x,a4,2(1x,e19.13))
      return
      end
      subroutine readsme(imin,imax,narm,idam,idamx,iusm)
!=======================================================================
!
!  READ THE HARMONICS FROM FILE lines
!  GIVES THE SMEAR AS A QUALITY FACTOR FOR THE DIFFERENT CASES
!
!  INPUT PARAMETERS:
!
!  NARM : NUMERO DI ARMONICHE DA ORDINARE
!  IR   : FLAG FOR COMPLEX OR REAL SIGNALS:
!         1 = INPUT SIGNAL IS REAL.
!         IF THE SIGNAL IS REAL IT IS TRASPOSED IN [0,0.5]
!
!  N.B. :  THE STRENGTH OF NEXT TO LEADING LINES ARE NORMALIZED AT THE
!          TUNE STRENGTH.
!  MIND :  THE ARRAYS TX,TY,TZ ARE USED IN ORDER NOT TO CHANGE
!          THE ARRAYS TUNEX,TUNEY,TUNEZ.
!
!  IOU  : FORTRAN UNIT OF THE OUTPUT
!
!  AUTHOR: R.BARTOLINI 21/08/1996
!
!=======================================================================
      implicit none
      integer i,idam,idamx,imax,imin,istep,iusm,j,k,ki,l,m,n,narm,ni
      double precision ex,fx,px,qsme,tx
      integer mterm
      parameter(mterm=300)
      integer mbpm
      double precision pieni,zero,one,two,pi,duepi,piph
      parameter(mbpm=1000,pieni=1d-17,zero=0d0,one=1d0,two=2d0)
      dimension tx(mterm)
      character*200 ch,ch1
 
! INITIALIZATION
      ch1=' '
 
      if(imin.le.imax) then
        istep=1
      else
        istep=-1
      endif
 
      rewind(30)
      do i=imin,imax,istep
        write(6,*)'ANALYZING FILE lines FOR SMEAR CALCULATION, CASE:',i
        qsme=0d0
        do ki=1,idam
 1        read(30,'(A)',end=2) ch
! Find start for each plane
          if(index(ch(1:10),'ANALYSIS').gt.0d0) then
            read(30,'()',end=2)
            read(30,'()',end=2)
            read(30,'()',end=2)
            read(30,'()',end=2)
            do ni=1,narm
              read(30,'(A)',end=2) ch
              if(ch.eq.ch1) goto 3
              read(ch,100,end=2)n,tx(n),px,fx,ex,l,m,k,j
              if(ki.eq.idamx) then
                if(ni.ge.2) then
                  qsme=qsme+px
                endif
              endif
            enddo
          else
            goto 1
          endif
        enddo
 3      continue
        write(iusm,*)i,qsme
      enddo
 2    continue
      rewind(30)
 
100   format(i3,1x,e19.13,2(1x,e18.12),1x,e17.11,:,4(1x,i3))
 
      return
      end
      subroutine readinv(imin,imax,narm,idam,iinv)
!=======================================================================
!
!  READ THE HARMONICS FROM FILE lines
!  GIVES THE INVARIANTS FOR THE 3 PLANES
!
!  INPUT PARAMETERS:
!
!  IMIN : MIN
!  IMAX : &MAX NUMBER OF CASES
!
!  NARM : NUMERO DI ARMONICHE DA ORDINARE
!
!  IDAM : NUMBER OF PLANES
!
!  IINV : OUTPUT UNIT
!
!  THEORY: A.BAZZANI, L.BONGINI, G.TURCHETTI, AIP 329, PAGE 120
!  AUTHOR: R.BARTOLINI/F.Schmidt 14/08/1997
!  last change 22/09/1997
!
!=======================================================================
      implicit none
      integer i,idam,ii,iinv,imax,imin,istep,j,k,ki,kii,kkini,l,lkini,m,&
     &mkini,n,narm,ni
      double precision ex,fx,psinv,px,pxinv,pyinv,tx
      integer mterm
      parameter(mterm=300)
      integer mbpm
      double precision pieni,zero,one,two,pi,duepi,piph
      parameter(mbpm=1000,pieni=1d-17,zero=0d0,one=1d0,two=2d0)
      dimension px(3,mterm)
      dimension l(3,mterm),m(3,mterm),k(3,mterm)
      character*200 ch,ch1
 
! INITIALIZATION
      ch1=' '
      if(imin.le.imax) then
        istep=1
      else
        istep=-1
      endif
 
      do i=1,3
        do ii=1,mterm
          l(i,ii)=0
          m(i,ii)=0
          k(i,ii)=0
        enddo
      enddo
 
      do i=1,3
        do ii=1,mterm
          px(i,ii)=zero
        enddo
      enddo
 
      rewind(30)
! study imax-imin+1 sets of data
      do i=imin,imax,istep
        write(6,*)'ANALYZING UNIT 30. INVARIANT CALCULATION, CASE:',i
! loop over planes
        do ki=1,idam
 1        read(30,'(A)',end=2) ch
! Find start for each plane
          if(index(ch(1:10),'ANALYSIS').gt.0d0) then
            read(30,'()',end=2)
            read(30,'()',end=2)
            read(30,'()',end=2)
            read(30,'()',end=2)
            do ni=1,narm
              read(30,'(A)',end=2) ch
! Detect end of data and exit to 3 to study next plane
              if(ch.eq.ch1) goto 3
              lkini=0
              mkini=0
              kkini=0
              read(ch,100)n,tx,px(ki,ni),fx,ex,lkini,mkini,kkini,j
! Reject new line if already present in analysed (ni-1) lines
              if(ni.gt.1) then
                do kii=1,ni-1
                  if(l(ki,kii).eq.lkini.and.m(ki,kii).eq.mkini.and.     &
     &k(ki,kii).eq.kkini) then
! true when line already exists => exit to 4 to look for new line
                    l(ki,ni)=0
                    m(ki,ni)=0
                    k(ki,ni)=0
                    goto 4
                  endif
                enddo
              endif
! set new coefficients
              l(ki,ni)=lkini
              m(ki,ni)=mkini
              k(ki,ni)=kkini
 4            continue
            enddo
          else
            goto 1
          endif
 3        continue
        enddo
 
!.......INVARIANTS CALCULATION
        pxinv=0d0
        pyinv=0d0
        psinv=0d0
        do n=1,narm
          pxinv=pxinv+l(1,n)*px(1,n)**2+l(2,n)*px(2,n)**2+              &
     &l(3,n)*px(3,n)**2
          pyinv=pyinv+m(1,n)*px(1,n)**2+m(2,n)*px(2,n)**2+              &
     &m(3,n)*px(3,n)**2
          psinv=psinv+k(1,n)*px(1,n)**2+k(2,n)*px(2,n)**2+              &
     &k(3,n)*px(3,n)**2
        enddo
        if(pxinv.lt.zero) then
          write(6,*)'Warning: Horizontal Invariant negative'
          pxinv=zero
        else
          pxinv=sqrt(pxinv)
        endif
        if(pyinv.lt.zero) then
          write(6,*)'Warning: Vertical Invariant negative'
          pyinv=zero
        else
          pyinv=sqrt(pyinv)
        endif
        if(psinv.lt.zero) then
          write(6,*)'Warning: Longitudinal Invariant negative'
          psinv=zero
        else
          psinv=sqrt(psinv)
        endif
 
        write(iinv,'(3(1X,E20.12))')pxinv,pyinv,psinv
 
      enddo
      rewind(30)
 2    continue
 
100   format(i3,1x,e19.13,2(1x,e18.12),1x,e17.11,:,4(1x,i3))
 
      return
      end
      subroutine readric(nfile,idam,ntwin,iconv)
!=======================================================================
!
! SUBROUTINE READRIC
!
! GIVEN A SIXTRACK BINARY FILE IT CONVERTS IT INTO ONE OR TWO ASCII
! ACCORDING TO NTWIN.
!
! NFILE IS THE UNIT OF THE FORT.NFILE TO BE PROCESSED
!
! THE OUTPUT ARE ALWAYS IN aux1 AND aux2
!
! N.B.: NTWIN IS AN OUTPUT PARAMETER AND THE NUMBER OF TURN
!       IS AUTOMATICALLY DETERMINED FROM READING THE WHOLE FILE.
!
!=======================================================================
      implicit none
      integer i,ia,icode,iconv,idam,ifipa,ilapa,iq,itopa,its6d,j,jq,    &
     &nerror,nfile,nfile0,nfile1,ntwin,numl,idummy(6)
      double precision b,c,c1,clo,clop,d,d1,dizu0,di0x,di0z,dip0x,dip0z,&
     &dmmac,dnms,dnumlr,dummy,e,e1,f,f1,g,g1,h,h1,p,p1,qwc,t,ta,tasum,  &
     &txyz,xyzv
      integer mbpm
      double precision pieni,zero,one,two,pi,duepi,piph
      parameter(mbpm=1000,pieni=1d-17,zero=0d0,one=1d0,two=2d0)
      character*80 sixtit,comment
      character*8 cdate,ctime,prgram
      dimension qwc(3),clo(3),clop(3)
      dimension ta(6,6),t(6,6),txyz(6),xyzv(6)
 
! INITIALIZATION
      dummy=zero
 
!---------------------------------------------------------------------
 
      nfile0=91
      nfile1=92
      rewind(nfile)
      rewind(nfile0)
      rewind(nfile1)
 
      read(nfile) sixtit,comment,cdate,ctime,                           &
     &prgram,ifipa,ilapa,itopa,icode,numl,qwc(1),qwc(2),qwc(3),         &
     &clo(1),clop(1),clo(2),clop(2),clo(3),clop(3),                     &
     &di0x,dip0x,di0z,dip0z,dummy,dummy,                                &
     &ta(1,1),ta(1,2),ta(1,3),                                          &
     &ta(1,4),ta(1,5),ta(1,6),                                          &
     &ta(2,1),ta(2,2),ta(2,3),                                          &
     &ta(2,4),ta(2,5),ta(2,6),                                          &
     &ta(3,1),ta(3,2),ta(3,3),                                          &
     &ta(3,4),ta(3,5),ta(3,6),                                          &
     &ta(4,1),ta(4,2),ta(4,3),                                          &
     &ta(4,4),ta(4,5),ta(4,6),                                          &
     &ta(5,1),ta(5,2),ta(5,3),                                          &
     &ta(5,4),ta(5,5),ta(5,6),                                          &
     &ta(6,1),ta(6,2),ta(6,3),                                          &
     &ta(6,4),ta(6,5),ta(6,6),                                          &
     &dmmac,dnms,dizu0,dnumlr,                                          &
     &dummy,dummy,dummy,dummy,dummy,dummy,                              &
     &dummy,dummy,dummy,dummy,dummy,dummy,                              &
     &dummy,dummy,dummy,dummy,dummy,dummy,                              &
     &dummy,dummy,dummy,dummy,dummy,dummy,                              &
     &dummy,dummy,dummy,dummy,dummy,dummy,                              &
     &dummy,dummy,dummy,dummy,dummy,dummy,                              &
     &dummy,dummy,dummy,dummy,dummy,dummy,                              &
     &dummy,dummy,dummy,dummy
 
!--REPAIR PARTIAL OR TOTAL INCOMPLETE TRANSFER MATRIX
 
      if(abs(ta(1,1)).le.pieni.and.abs(ta(2,2)).le.pieni) then
        ta(1,1)=one
        ta(2,2)=one
      endif
      if(abs(ta(3,3)).le.pieni.and.abs(ta(4,4)).le.pieni) then
        ta(3,3)=one
        ta(4,4)=one
      endif
      if(abs(ta(5,5)).le.pieni.and.abs(ta(6,6)).le.pieni) then
        ta(5,5)=one
        ta(6,6)=one
      endif
      its6d=0
 
!--Convert if requested
 
      if(iconv.ne.1) then
        do i=1,6
          do j=1,6
            if(i.ne.j) then
              ta(i,j)=zero
            else
              ta(i,j)=one
            endif
          enddo
        enddo
      endif
 
!--TEST IF TRANSFER MATRIX HAS LONGITUDINAL PART
 
      do i=1,6
        tasum=tasum+abs(ta(i,5))+abs(ta(i,6))
      enddo
      do i=1,4
        tasum=tasum+abs(ta(5,i))+abs(ta(6,i))
      enddo
      tasum=tasum-two
      if(abs(tasum).ge.pieni) its6d=1
 
!--INVERT TRANSFER MATRIX
 
      do i=1,6
        do j=1,6
          t(i,j)=ta(j,i)
        enddo
      enddo
      call dinv(6,t,6,idummy,nerror)
 
!--DETERMINE ORDER OF FREEDOM OF TRACKED CASE
 
!      IF(ICODE.EQ.1.OR.ICODE.EQ.2.OR.ICODE.EQ.4) IDAM=1
!      IF(ICODE.EQ.3.OR.ICODE.EQ.5.OR.ICODE.EQ.6) IDAM=2
!      IF(ICODE.EQ.7) IDAM=3
 
!--CHECK IF 1 OR 2 SETS OF COORDINATES ARE WRITTEN
 
      ntwin=1
      if(ilapa.ne.ifipa) ntwin=2
      write(6,*)'ANALYZING THE BINARY SIXTRACK OUTPUT'
      do i=1,numl+1
        if(ntwin.eq.1)                                                  &
     &read(nfile,end=999) ia,ifipa,b,c,d,e,f,g,h,p
        if(ntwin.eq.2)                                                  &
     &read(nfile,end=999) ia,ifipa,b,c,d,e,f,g,h,p,                     &
     &ilapa,b,c1,d1,e1,f1,g1,h1,p1
 
!--SUBTRACT CLOSED ORBIT
 
        c=c-clo(1)
        d=d-clop(1)
        e=e-clo(2)
        f=f-clop(2)
        g=g-clo(3)
        h=h-clop(3)
 
!--!! SUBTRACT DISPERSION ONLY IF THE LONGITUDINAL
!--PART OF THE TRANSFER MATRIX HAS NOT BEEN
!--CALCULATED BUT THERE ARE SYNCHROTRON OSCILLATIONS
 
        if(icode.ge.4.and.its6d.eq.0) then
          c=c-di0x*h
          d=d-dip0x*h
          e=e-di0z*h
          f=f-dip0z*h
        endif
        xyzv(1)=c
        xyzv(2)=d
        xyzv(3)=e
        xyzv(4)=f
        xyzv(5)=g
        xyzv(6)=h
 
!--TRANSFER FIRST SET OF DATA
 
        do iq=1,6
          txyz(iq)=zero
          do jq=1,6
            txyz(iq)=txyz(iq)+t(jq,iq)*xyzv(jq)
          enddo
        enddo
        c=txyz(1)
        d=txyz(2)
        e=txyz(3)
        f=txyz(4)
        g=txyz(5)
 
!--!!IMPORTANT THE MOMENTUM HAS A TRIVIAL SCALING!!
 
        h=txyz(6)*1d3
        if(ntwin.eq.2) then
          c1=c1-clo(1)
          d1=d1-clop(1)
          e1=e1-clo(2)
          f1=f1-clop(2)
          g1=g1-clo(3)
          h1=h1-clop(3)
 
!--   TREAT SECOND SET OF DATA IF THEY EXIST
 
          if(icode.ge.4.and.its6d.eq.0) then
            c1=c1-di0x*h1
            d1=d1-dip0x*h1
            e1=e1-di0z*h1
            f1=f1-dip0z*h1
          endif
          xyzv(1)=c1
          xyzv(2)=d1
          xyzv(3)=e1
          xyzv(4)=f1
          xyzv(5)=g1
          xyzv(6)=h1
 
!--TRANSFER SECOND SET OF DATA
 
          do iq=1,6
            txyz(iq)=zero
            do jq=1,6
              txyz(iq)=txyz(iq)+t(jq,iq)*xyzv(jq)
            enddo
          enddo
          c1=txyz(1)
          d1=txyz(2)
          e1=txyz(3)
          f1=txyz(4)
          g1=txyz(5)
!--!!IMPORTANT THE MOMENTUM HAS A TRIVIAL SCALING!!
          h1=txyz(6)*1d3
        endif
!--OBEY DEMANDED ORDER OF PHASE SPACE
        if(idam.lt.3) then
          g=zero
          h=zero
          g1=zero
          h1=zero
        endif
        if(idam.eq.1) then
          e=zero
          f=zero
          e1=zero
          f1=zero
        endif
        if(ntwin.eq.1) then
          write(nfile0,'(6(G21.15,1X))') c,d,e,f,g,h
        endif
        if(ntwin.eq.2) then
          write(nfile0,'(6(G21.15,1X))') c,d,e,f,g,h
          write(nfile1,'(6(G21.15,1X))') c1,d1,e1,f1,g1,h1
        endif
      enddo
 999  continue
      rewind(nfile)
      rewind(nfile0)
      rewind(nfile1)
 
      return
      end
      subroutine writeric(nfile,ntwin,numl,iconv)
!=======================================================================
!
! SUBROUTINE WRITERIC
!
! READS THE PROCESSED ASCII FILE FROM aux1 AND aux2 IF THERE
! IS A TWIN PARTICLE.
!
! AND PLUGS THEM INTO A SIXTRACK BINARY FILE IN UNIT NFILE.
!
!
!=======================================================================
      implicit none
      integer i,i11,i22,ia,icode,iconv,ifipa,ilapa,iq,itopa,its6d,j,jq, &
     &ndiff,nfile,ntwin,numl,num1
      double precision b,c,c1,clo,clop,d,d1,di0x,di0z,dip0x,dip0z,dizu0,&
     &dmmac,dnms,dnumlr,dummy,e,e1,f,f1,g,g1,h,h1,p,p1,qwc,t,ta,tasum,  &
     &txyz,xyzv
      integer mbpm
      double precision pieni,zero,one,two,pi,duepi,piph
      parameter(mbpm=1000,pieni=1d-17,zero=0d0,one=1d0,two=2d0)
      character*80 sixtit,comment
      character*8 cdate,ctime,prgram
      dimension qwc(3),clo(3),clop(3)
      dimension ta(6,6),t(6,6),txyz(6),xyzv(6)
 
!.....READS AGAIN THE HEADER AND CHECK FOR THE ACTUAL NUMBER OF TURNS
 
      rewind(nfile)
 
      read(nfile) sixtit,comment,cdate,ctime,                           &
     &prgram,ifipa,ilapa,itopa,icode,num1,qwc(1),qwc(2),qwc(3),         &
     &clo(1),clop(1),clo(2),clop(2),clo(3),clop(3),                     &
     &di0x,dip0x,di0z,dip0z,dummy,dummy,                                &
     &ta(1,1),ta(1,2),ta(1,3),                                          &
     &ta(1,4),ta(1,5),ta(1,6),                                          &
     &ta(2,1),ta(2,2),ta(2,3),                                          &
     &ta(2,4),ta(2,5),ta(2,6),                                          &
     &ta(3,1),ta(3,2),ta(3,3),                                          &
     &ta(3,4),ta(3,5),ta(3,6),                                          &
     &ta(4,1),ta(4,2),ta(4,3),                                          &
     &ta(4,4),ta(4,5),ta(4,6),                                          &
     &ta(5,1),ta(5,2),ta(5,3),                                          &
     &ta(5,4),ta(5,5),ta(5,6),                                          &
     &ta(6,1),ta(6,2),ta(6,3),                                          &
     &ta(6,4),ta(6,5),ta(6,6),                                          &
     &dmmac,dnms,dizu0,dnumlr,                                          &
     &dummy,dummy,dummy,dummy,dummy,dummy,                              &
     &dummy,dummy,dummy,dummy,dummy,dummy,                              &
     &dummy,dummy,dummy,dummy,dummy,dummy,                              &
     &dummy,dummy,dummy,dummy,dummy,dummy,                              &
     &dummy,dummy,dummy,dummy,dummy,dummy,                              &
     &dummy,dummy,dummy,dummy,dummy,dummy,                              &
     &dummy,dummy,dummy,dummy,dummy,dummy,                              &
     &dummy,dummy,dummy,dummy
 
      ntwin=1
      if(ilapa.ne.ifipa) ntwin=2
      if(ntwin.eq.1) then
        read(nfile,end=999) i11,ifipa,b,c,d,e,f,g,h,p
        p1=p
      else
        read(nfile,end=999) ia,ifipa,b,c,d,e,f,g,h,p,                   &
     &ilapa,b,c1,d1,e1,f1,g1,h1,p1
      endif
 
      read(nfile,end=999) i22
      ndiff=i22-i11
      rewind(nfile)
      comment='SUSSIX TREATED DATA'
 
      rewind(91)
      rewind(92)
      rewind(nfile)
!.....WRITES THE SAME HEADER AS THE INPUT FILE
      write(nfile) sixtit,comment,cdate,ctime,                          &
     &prgram,ifipa,ilapa,itopa,icode,num1,qwc(1),qwc(2),qwc(3),         &
     &clo(1),clop(1),clo(2),clop(2),clo(3),clop(3),                     &
     &di0x,dip0x,di0z,dip0z,dummy,dummy,                                &
     &ta(1,1),ta(1,2),ta(1,3),                                          &
     &ta(1,4),ta(1,5),ta(1,6),                                          &
     &ta(2,1),ta(2,2),ta(2,3),                                          &
     &ta(2,4),ta(2,5),ta(2,6),                                          &
     &ta(3,1),ta(3,2),ta(3,3),                                          &
     &ta(3,4),ta(3,5),ta(3,6),                                          &
     &ta(4,1),ta(4,2),ta(4,3),                                          &
     &ta(4,4),ta(4,5),ta(4,6),                                          &
     &ta(5,1),ta(5,2),ta(5,3),                                          &
     &ta(5,4),ta(5,5),ta(5,6),                                          &
     &ta(6,1),ta(6,2),ta(6,3),                                          &
     &ta(6,4),ta(6,5),ta(6,6),                                          &
     &dmmac,dnms,dizu0,dnumlr,                                          &
     &dummy,dummy,dummy,dummy,dummy,dummy,                              &
     &dummy,dummy,dummy,dummy,dummy,dummy,                              &
     &dummy,dummy,dummy,dummy,dummy,dummy,                              &
     &dummy,dummy,dummy,dummy,dummy,dummy,                              &
     &dummy,dummy,dummy,dummy,dummy,dummy,                              &
     &dummy,dummy,dummy,dummy,dummy,dummy,                              &
     &dummy,dummy,dummy,dummy,dummy,dummy,                              &
     &dummy,dummy,dummy,dummy
 
!---------------------------------------------------------------------
 
!....CONVERTS THE DATA FROM NORMALIZED COORDINATE TO PHISICAL COORDINATES
 
!--REPAIR PARTIAL OR TOTAL INCOMPLETE TRANSFER MATRIX
 
      if(abs(ta(1,1)).le.pieni.and.abs(ta(2,2)).le.pieni) then
        ta(1,1)=one
        ta(2,2)=one
      endif
      if(abs(ta(3,3)).le.pieni.and.abs(ta(4,4)).le.pieni) then
        ta(3,3)=one
        ta(4,4)=one
      endif
      if(abs(ta(5,5)).le.pieni.and.abs(ta(6,6)).le.pieni) then
        ta(5,5)=one
        ta(6,6)=one
      endif
      its6d=0
 
!--Convert if requested
 
      if(iconv.ne.1) then
        do i=1,6
          do j=1,6
            if(i.ne.j) then
              ta(i,j)=zero
            else
              ta(i,j)=one
            endif
          enddo
        enddo
      endif
 
!--TEST IF TRANSFER MATRIX HAS LONGITUDINAL PART
 
      do i=1,6
        tasum=tasum+abs(ta(i,5))+abs(ta(i,6))
      enddo
      do i=1,4
        tasum=tasum+abs(ta(5,i))+abs(ta(6,i))
      enddo
      tasum=tasum-two
      if(abs(tasum).ge.pieni) its6d=1
 
!--TRANSPOSE MATRIX
 
      do i=1,6
        do j=1,6
          t(i,j)=ta(j,i)
        enddo
      enddo
 
!--CHECK IF 1 OR 2 SETS OF COORDINATES ARE WRITTEN
 
      write(6,*)'WRITING SUSSIX TREATED DATA INTO THE BINARY FILE'
      do i=1,numl
        if(ntwin.eq.1) read(91,*,end=999)c,d,e,f,g,h
        if(ntwin.eq.2) then
          read(91,*,end=999) c,d,e,f,g,h
          read(92,*,end=999) c1,d1,e1,f1,g1,h1
        endif
 
!--TRANSFER FIRST SET OF DATA
        xyzv(1)=c
        xyzv(2)=d
        xyzv(3)=e
        xyzv(4)=f
        xyzv(5)=g
!--!!IMPORTANT THE MOMENTUM HAS A TRIVIAL SCALING!!
        xyzv(6)=h/1d3
 
        do iq=1,6
          txyz(iq)=zero
          do jq=1,6
            txyz(iq)=txyz(iq)+t(jq,iq)*xyzv(jq)
          enddo
        enddo
        c=txyz(1)
        d=txyz(2)
        e=txyz(3)
        f=txyz(4)
        g=txyz(5)
        h=txyz(6)
 
!--!! ADD DISPERSION ONLY IF THE LONGITUDINAL
!--PART OF THE TRANSFER MATRIX HAS NOT BEEN
!--CALCULATED BUT THERE ARE SYNCHROTRON OSCILLATIONS
 
        if(icode.ge.4.and.its6d.eq.0) then
          c=c+di0x*h
          d=d+dip0x*h
          e=e+di0z*h
          f=f+dip0z*h
        endif
 
!--ADD CLOSED ORBIT
 
        c=c+clo(1)
        d=d+clop(1)
        e=e+clo(2)
        f=f+clop(2)
        g=g+clo(3)
        h=h+clop(3)
 
        if(ntwin.eq.2) then
!--TRANSFER SECOND SET OF DATA IF THEY EXIST
          xyzv(1)=c1
          xyzv(2)=d1
          xyzv(3)=e1
          xyzv(4)=f1
          xyzv(5)=g1
!--!!IMPORTANT THE MOMENTUM HAS A TRIVIAL SCALING!!
          xyzv(6)=h1/1d3
 
          do iq=1,6
            txyz(iq)=zero
            do jq=1,6
              txyz(iq)=txyz(iq)+t(jq,iq)*xyzv(jq)
            enddo
          enddo
          c1=txyz(1)
          d1=txyz(2)
          e1=txyz(3)
          f1=txyz(4)
          g1=txyz(5)
          h1=txyz(6)
 
!--!! ADD DISPERSION ONLY IF THE LONGITUDINAL
!--PART OF THE TRANSFER MATRIX HAS NOT BEEN
!--CALCULATED BUT THERE ARE SYNCHROTRON OSCILLATIONS
 
          if(icode.ge.4.and.its6d.eq.0) then
            c1=c1+di0x*h1
            d1=d1+dip0x*h1
            e1=e1+di0z*h1
            f1=f1+dip0z*h1
          endif
 
!--ADD CLOSED ORBIT
 
          c1=c1+clo(1)
          d1=d1+clop(1)
          e1=e1+clo(2)
          f1=f1+clop(2)
          g1=g1+clo(3)
          h1=h1+clop(3)
 
        endif
 
!.....WRITES THE DATA IN BINARY FORMAT ON UNIT NFILE
 
        ia=(i-1)*ndiff
        if(ntwin.eq.1) then
          write(nfile) ia,ifipa,b,c,d,e,f,g,h,p
        elseif(ntwin.eq.2) then
          write(nfile) ia,ifipa,b,c,d,e,f,g,h,p,                        &
     &ilapa,b,c1,d1,e1,f1,g1,h1,p1
        endif
      enddo
 999  continue
      rewind(nfile)
      rewind(91)
      rewind(92)
 
      return
      end
      subroutine kerset(ercode,lgfile,limitm,limitr)
      implicit none
      integer i,kounte,l,lgfile,limitm,limitr,log,logf
      parameter(kounte = 27)
      character*6         ercode,   code(kounte)
      logical             mflag,    rflag
      integer             kntm(kounte),       kntr(kounte)
!-----------------------------------------------------------------------
      data      logf      /  0  /
      data      code(1), kntm(1), kntr(1)  / 'C204.1', 255, 255 /
      data      code(2), kntm(2), kntr(2)  / 'C204.2', 255, 255 /
      data      code(3), kntm(3), kntr(3)  / 'C204.3', 255, 255 /
      data      code(4), kntm(4), kntr(4)  / 'C205.1', 255, 255 /
      data      code(5), kntm(5), kntr(5)  / 'C205.2', 255, 255 /
      data      code(6), kntm(6), kntr(6)  / 'C305.1', 255, 255 /
      data      code(7), kntm(7), kntr(7)  / 'C308.1', 255, 255 /
      data      code(8), kntm(8), kntr(8)  / 'C312.1', 255, 255 /
      data      code(9), kntm(9), kntr(9)  / 'C313.1', 255, 255 /
      data      code(10),kntm(10),kntr(10) / 'C336.1', 255, 255 /
      data      code(11),kntm(11),kntr(11) / 'C337.1', 255, 255 /
      data      code(12),kntm(12),kntr(12) / 'C341.1', 255, 255 /
      data      code(13),kntm(13),kntr(13) / 'D103.1', 255, 255 /
      data      code(14),kntm(14),kntr(14) / 'D106.1', 255, 255 /
      data      code(15),kntm(15),kntr(15) / 'D209.1', 255, 255 /
      data      code(16),kntm(16),kntr(16) / 'D509.1', 255, 255 /
      data      code(17),kntm(17),kntr(17) / 'E100.1', 255, 255 /
      data      code(18),kntm(18),kntr(18) / 'E104.1', 255, 255 /
      data      code(19),kntm(19),kntr(19) / 'E105.1', 255, 255 /
      data      code(20),kntm(20),kntr(20) / 'E208.1', 255, 255 /
      data      code(21),kntm(21),kntr(21) / 'E208.2', 255, 255 /
      data      code(22),kntm(22),kntr(22) / 'F010.1', 255,   0 /
      data      code(23),kntm(23),kntr(23) / 'F011.1', 255,   0 /
      data      code(24),kntm(24),kntr(24) / 'F012.1', 255,   0 /
      data      code(25),kntm(25),kntr(25) / 'F406.1', 255,   0 /
      data      code(26),kntm(26),kntr(26) / 'G100.1', 255, 255 /
      data      code(27),kntm(27),kntr(27) / 'G100.2', 255, 255 /
      save
!-----------------------------------------------------------------------
      logf  =  lgfile
      l  =  0
      if(ercode .ne. ' ')  then
        do l = 1, 6
          if(ercode(1:l) .eq. ercode)  goto 12
        enddo
 12     continue
      endif
      do i  =  1, kounte
        if(l .eq. 0)  goto 13
        if(code(i)(1:l) .ne. ercode(1:l))  goto 14
 13     if(limitm.ge.0) kntm(i)  =  limitm
        if(limitr.ge.0) kntr(i)  =  limitr
 14     continue
      enddo
      return
      entry kermtr(ercode,log,mflag,rflag)
      log  =  logf
      do i  =  1, kounte
        if(ercode .eq. code(i))  goto 21
      enddo
      write(*,1000)  ercode
      call abend
      return
 21   rflag  =  kntr(i) .ge. 1
      if(rflag  .and.  (kntr(i) .lt. 255))  kntr(i)  =  kntr(i) - 1
      mflag  =  kntm(i) .ge. 1
      if(mflag  .and.  (kntm(i) .lt. 255))  kntm(i)  =  kntm(i) - 1
      if(.not. rflag)  then
        if(logf .lt. 1)  then
          write(*,1001)  code(i)
        else
          write(logf,1001)  code(i)
        endif
      endif
      if(mflag .and. rflag)  then
        if(logf .lt. 1)  then
          write(*,1002)  code(i)
        else
          write(logf,1002)  code(i)
        endif
      endif
      return
1000  format(' KERNLIB LIBRARY ERROR. ' /                               &
     &' ERROR CODE ',a6,' NOT RECOGNIZED BY KERMTR',                    &
     &' ERROR MONITOR. RUN ABORTED.')
1001  format(/' ***** RUN TERMINATED BY CERN LIBRARY ERROR ',           &
     &'CONDITION ',a6)
1002  format(/' ***** CERN LIBRARY ERROR CONDITION ',a6)
      end
      subroutine dinv(n,a,idim,ir,ifail)
!-----------------------------------------------------------------------
!
!     ******************************************************************
!
!     REPLACES A BY ITS INVERSE.
!
!     (PARAMETERS AS FOR DEQINV.)
!
!     CALLS ... DFACT, DFINV, F010PR, ABEND.
!
!     ******************************************************************
!-----------------------------------------------------------------------
      implicit none
      integer idim,ifail,jfail,k,kprnt,n
      integer ir
      real t1,t2,t3
      double precision a,det,temp,s,c11,c12,c13,c21,c22,c23,c31,c32,c33
      character*6 name
      dimension ir(n),a(idim,n)
      data name/'DINV'/,kprnt/0/
      save
!-----------------------------------------------------------------------
!
!  TEST FOR PARAMETER ERRORS.
!
      if((n.lt.1).or.(n.gt.idim)) goto 7
!
!  TEST FOR N.LE.3.
!
      if(n.gt.3) goto 6
      ifail=0
      if(n.lt.3) goto 4
!
!  N=3 CASE.
!
!     COMPUTE COFACTORS.
      c11=a(2,2)*a(3,3)-a(2,3)*a(3,2)
      c12=a(2,3)*a(3,1)-a(2,1)*a(3,3)
      c13=a(2,1)*a(3,2)-a(2,2)*a(3,1)
      c21=a(3,2)*a(1,3)-a(3,3)*a(1,2)
      c22=a(3,3)*a(1,1)-a(3,1)*a(1,3)
      c23=a(3,1)*a(1,2)-a(3,2)*a(1,1)
      c31=a(1,2)*a(2,3)-a(1,3)*a(2,2)
      c32=a(1,3)*a(2,1)-a(1,1)*a(2,3)
      c33=a(1,1)*a(2,2)-a(1,2)*a(2,1)
      t1=abs(sngl(a(1,1)))
      t2=abs(sngl(a(2,1)))
      t3=abs(sngl(a(3,1)))
!
!     (SET TEMP=PIVOT AND DET=PIVOT*DET.)
      if(t1.ge.t2) goto 1
      if(t3.ge.t2) goto 2
!        (PIVOT IS A21)
      temp=a(2,1)
      det=c13*c32-c12*c33
      goto 3
    1 if(t3.ge.t1) goto 2
!     (PIVOT IS A11)
      temp=a(1,1)
      det=c22*c33-c23*c32
      goto 3
!     (PIVOT IS A31)
 2    temp=a(3,1)
      det=c23*c12-c22*c13
!
!     SET ELEMENTS OF INVERSE IN A.
    3 if(det.eq.0d0) goto 8
      s=temp/det
      a(1,1)=s*c11
      a(1,2)=s*c21
      a(1,3)=s*c31
      a(2,1)=s*c12
      a(2,2)=s*c22
      a(2,3)=s*c32
      a(3,1)=s*c13
      a(3,2)=s*c23
      a(3,3)=s*c33
      return
!
    4 if(n.lt.2) goto 5
!
!  N=2 CASE BY CRAMERS RULE.
!
      det=a(1,1)*a(2,2)-a(1,2)*a(2,1)
      if(det.eq.0d0) goto 8
      s=1d0/det
      c11   =s*a(2,2)
      a(1,2)=-s*a(1,2)
      a(2,1)=-s*a(2,1)
      a(2,2)=s*a(1,1)
      a(1,1)=c11
      return
!
!  N=1 CASE.
!
    5 if(a(1,1).eq.0d0) goto 8
      a(1,1)=1d0/a(1,1)
      return
!
!  N.GT.3 CASES.  FACTORIZE MATRIX AND INVERT.
!
    6 call dfact(n,a,idim,ir,ifail,det,jfail)
      if(ifail.ne.0) return
      call dfinv(n,a,idim,ir)
      return
!
!  ERROR EXITS.
!
    7 ifail=+1
      call f010pr(name,n,idim,k,kprnt)
      return
!
    8 ifail=-1
      return
!
      end
      subroutine f010pr(name,n,idim,k,kprnt)
!     ******************************************************************
!
!     PRINT ROUTINE FOR PARAMETER ERRORS IN MATRIX SUBROUTINES $EQINV,
!     $EQN, $INV (WHERE $ IS A LETTER SPECIFYING THE ARITHMETIC TYPE).
!
!     NAME         (CHARACTER*6) NAME OF THE CALLING ROUTINE.
!
!     N,IDIM,K     PARAMETERS OF THE CALLING ROUTINE (WITH K=0 IF K IS
!                  NOT TO BE PRINTED).
!
!     KPRNT        PRINT FLAG FOR K (K IS NOT PRINTED IF KPRNT=0).
!
!     ******************************************************************
      implicit none
      integer idim,k,kprnt,lgfile,n
      character*6 name
      logical mflag,rflag
!-----------------------------------------------------------------------
      call kermtr('F010.1',lgfile,mflag,rflag)
      if(mflag) then
        if(lgfile.eq.0)  then
          if(kprnt.eq.0) write(*,2000) name,n,idim
          if(kprnt.ne.0) write(*,2001) name,n,idim,k
        else
          if(kprnt.eq.0) write(lgfile,2000) name,n,idim
          if(kprnt.ne.0) write(lgfile,2001) name,n,idim,k
        endif
      endif
      if(.not. rflag) call abend
      return
!
 2000 format( 7x, 11x,'subroutine ', a6, 14x,' ... parameter',          &
     &29x,' error (n.lt.1 or n.gt.idim).',                              &
     &6x, 3x,'n =', i4, 6x, 6x,'idim =', i4, 1x,'. ')
 2001 format( 7x, 11x,'subroutine ', a6, 14x,' ... parameter',          &
     &39x,' error (n.lt.1 or n.gt.idim or k.lt.1).',                    &
     &6x, 3x,'n =', i4, 6x, 6x,'idim =', i4, 6x, 3x,                    &
     &'k =', i4, 1x,'. ')
      end
      subroutine dfact(n,a,idim,ir,ifail,det,jfail)
      implicit none
      integer i,idim,ifail,imposs,ir,j,jfail,jm1,jover,jp1,jrange,      &
     &junder,k,l,n,normal,nxch
      real g1,g2,p,q,t
      double precision a,det,zero,one,s11,s12,tf
      character*6         hname
      dimension ir(*),a(idim,*)
      data      g1, g2              /  1.e-37,  1.e37  /
      data      hname               /  ' DFACT'  /
      data      zero, one           /  0.d0, 1.d0  /
      data      normal, imposs      /  0, -1  /
      data      jrange, jover, junder  /  0, +1, -1  /
!-----------------------------------------------------------------------
      if(idim .ge. n  .and.  n .gt. 0)  goto 110
      call tmprnt(hname,n,idim,0)
      return
 110  ifail  =  normal
      jfail  =  jrange
      nxch   =  0
      det    =  one
      do j  =  1, n
 120    k  =  j
        p  =  abs(sngl(a(j,j)))
        if(j .eq. n)  goto 122
        jp1  =  j+1
        do i  =  jp1, n
          q  =  abs(sngl(a(i,j)))
          if(q .le. p)  goto 121
          k  =  i
          p  =  q
 121      continue
        enddo
        if(k .ne. j)  goto 123
 122    if(p .gt. 0.)  goto 130
        det    =  zero
        ifail  =  imposs
        jfail  =  jrange
        return
 123    do l  =  1, n
          tf      =  a(j,l)
          a(j,l)  =  a(k,l)
          a(k,l)  =  tf
        enddo
        nxch      =  nxch + 1
        ir(nxch)  =  j*2**12 + k
 130    det     =  det * a(j,j)
        a(j,j)  =  one / a(j,j)
        t  =  abs(sngl(det))
        if(t .lt. g1)  then
          det    =  zero
          if(jfail .eq. jrange)  jfail  =  junder
        elseif(t .gt. g2)  then
          det    =  one
          if(jfail .eq. jrange)  jfail  =  jover
        endif
        if(j .eq. n)  goto 144
        jm1  =  j-1
        jp1  =  j+1
        do k  =  jp1, n
          s11  =  -a(j,k)
          s12  =  -a(k,j+1)
          if(j .eq. 1)  goto 142
          do i  =  1, jm1
            s11  =  a(i,k)*a(j,i)+s11
            s12  =  a(i,j+1)*a(k,i)+s12
          enddo
 142      a(j,k)    =  -s11 * a(j,j)
          a(k,j+1)  =  -(a(j,j+1)*a(k,j)+s12)
        enddo
 144    continue
      enddo
 150  if(mod(nxch,2) .ne. 0)  det  =  -det
      if(jfail .ne. jrange)   det  =  zero
      ir(n)  =  nxch
      return
      end
      subroutine dfinv(n,a,idim,ir)
      implicit none
      integer i,idim,ij,im2,ir,j,k,m,n,nm1,nmi,nxch
      double precision a,s31,s32,s33,s34,ti,zero
      character*6 hname
      dimension ir(*),a(idim,*)
      data      hname               /  ' DFINV'  /
      data      zero      /  0.d0  /
!-----------------------------------------------------------------------
      if(idim .ge. n  .and.  n .gt. 0)  goto 310
      call tmprnt(hname,n,idim,0)
      return
 310  if(n .eq. 1)  return
      a(2,1)  =  -a(2,2) * a(1,1)*a(2,1)
      a(1,2)  =  -a(1,2)
      if(n .eq. 2)  goto 330
      do i  =  3, n
        im2  =  i-2
        do j  =  1, im2
          s31  =  zero
          s32  =  a(j,i)
          do k  =  j, im2
            s31  =  a(k,j)*a(i,k)+s31
            s32  =  a(j,k+1)*a(k+1,i)+s32
          enddo
          a(i,j)  =  -a(i,i) * (a(i-1,j)*a(i,i-1)+s31)
          a(j,i)  =  -s32
        enddo
        a(i,i-1)  =  -a(i,i) * a(i-1,i-1)*a(i,i-1)
        a(i-1,i)  =  -a(i-1,i)
      enddo
 330  nm1  =  n-1
      do i  =  1, nm1
        nmi  =  n-i
        do j  =  1, i
          s33  =  a(i,j)
          do k  =  1, nmi
            s33  =  a(i+k,j)*a(i,i+k)+s33
          enddo
          a(i,j)  =  s33
        enddo
        do j  =  1, nmi
          s34  =  zero
          do k  =  j, nmi
            s34  =  a(i+k,i+j)*a(i,i+k)*s34
          enddo
          a(i,i+j)  =  s34
        enddo
      enddo
      nxch  =  ir(n)
      if(nxch .eq. 0)  return
      do m  =  1, nxch
        k   =  nxch - m+1
        ij  =  ir(k)
        i   =  ij / 4096
        j   =  mod(ij,4096)
        do k  =  1, n
          ti      =  a(k,i)
          a(k,i)  =  a(k,j)
          a(k,j)  =  ti
        enddo
      enddo
      return
      end
      subroutine tmprnt(name,n,idim,k)
      implicit none
      integer idim,k,lgfile,n
      character*6 name
      logical mflag,rflag
!-----------------------------------------------------------------------
      if(name(2:2) .eq. 'S') then
        call kermtr('F012.1',lgfile,mflag,rflag)
      else
        call kermtr('F011.1',lgfile,mflag,rflag)
      endif
      if(mflag) then
        if(lgfile .eq. 0) then
          if(name(3:6) .eq. 'FEQN') then
            write(*,1002) name, n, idim, k
          else
            write(*,1001) name, n, idim
          endif
        else
          if(name(3:6) .eq. 'FEQN') then
            write(lgfile,1002) name, n, idim, k
          else
            write(lgfile,1001) name, n, idim
          endif
        endif
      endif
      if(.not. rflag) call abend
      return
1001  format(7x, 31x,' parameter error in subroutine ', a6,             &
     &27x,' ... (n.lt.1 or idim.lt.n).',                                &
     &5x, 3x,'n =', i4, 5x, 6x,'idim =', i4, 1x,'. ')
1002  format(7x, 31x,' parameter error in subroutine ', a6,             &
     &37x,' ... (n.lt.1 or idim.lt.n or k.lt.1).',                      &
     &5x, 3x,'n =', i4, 5x, 6x,'idim =', i4, 5x, 3x,                    &
     &'k =', i4,1x,'.')
      end
      subroutine abend
      implicit none
      print*,'abend called ==> problem'
      return
      end
      SUBROUTINE CFFT(A,MSIGN)
      implicit none
      integer m,n,msign,nv2,nm1,i,j,k,le,le1,l,ip
      real c,s
      COMPLEX A,U,W,T
      integer maxiter
      parameter(maxiter=100000)
      dimension a(maxiter)
 
      IF(MSIGN.EQ.0) RETURN
      M=IABS(MSIGN)
      N=2**M
      NV2=N/2
      NM1=N-1
      J=1
      DO I=1,NM1
        IF(I.GE.J) GO TO 5
        T=A(J)
        A(J)=A(I)
        A(I)=T
 5      K=NV2
 6      IF(K.lt.J) then
          J=J-K
          K=K/2
          GO TO 6
        endif
        J=J+K
      enddo
      DO I=1,N,2
        T=A(I+1)
        A(I+1)=A(I)-T
        A(I )=A(I)+T
      enddo
      IF(M.EQ.1) RETURN
      C=0.
      S=ISIGN(1,MSIGN)
      LE=2
      DO L=2,M
        W=CMPLX(C,S)
        U=W
        C=SQRT(C*.5+.5)
        S=AIMAG(W)/(C+C)
        LE1=LE
        LE=LE1+LE1
        DO I=1,N,LE
          IP=I+LE1
          T=A(IP)
          A(IP)=A(I)-T
          A(I) =A(I)+T
        enddo
        DO J=2,LE1
          DO I=J,N,LE
            IP=I+LE1
            T=A(IP)*U
            A(IP)=A(I)-T
            A(I) =A(I)+T
          enddo
          U=U*W
        enddo
      enddo
      RETURN
      END
