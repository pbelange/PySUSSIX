      subroutine zfunr(z,maxn,tunea1,deltat,maxiter)
!=======================================================================
! Streamlined version of the original zfunr found in sussixNoO.f
!
! AUTHOR:     A. BAZZANI - BOLOGNA UNIVERSITY

! P. Belanger + chatGPT comments:        
! Initialization:
!     tunea1: Initial estimate of the tune.
!     deltat: Step size for tuning parameter updates.
!     zu: Complex unit.

! Iteration Loop (do ntest=1, 10):
!     tunea2: Update the tune estimate.
!     ztune2: Exponential term based on the updated tune.
!     calcr(ztune2, zf, z, maxn): Evaluate the function z(tune) with the updated tune.
!     calcr(ztune2, zfd, zd, maxn): Evaluate the derivative of z(tune) with respect to the tune.
!     dtunea2: Calculate the product of the function value and its derivative.
!     Check if the signs of dtunea1 and dtunea2 are suitable for the Newton method.

! Newton-like Iterative Loop (do ncont=1,100):
!     Update the tune estimate (tune3) based on the Newton method formula.
!     Evaluate the function (z(tune3)) and its derivative (z'(tune3)).
!     Adjust the tune estimates (tune1 and tune2) based on the comparison of signs.
!     Check for convergence (abs(tune2-tune1)) and exit the loop if satisfied.

! Result Recording:
!    ! Record the tuned estimates in the tunetest array.
!    ! The iterations involve updating the tune estimate based on the Newton method and checking 
!    !  for convergence. This process is part of the overall goal of refining the tune estimate in the zfunr subroutine.
!
!=======================================================================
      implicit none
      integer maxn,nc,ncont,nd,ntest,num
      double precision deltat,dtune1,dtune2,dtune3,dtunea1,dtunea2,     &
     &err,ratio,tune,tune1,tune2,tune3,tunea1,tunea2,tunetest,tuneval,  &
     &tunevmax
      double complex z,zd,zf,zfd,ztune,ztune1,ztune2,ztune3,zu,zw
      integer maxiter
      !parameter(maxiter=100000)
      integer mbpm
      double precision pieni,zero,one,two,pi,duepi,piph
      parameter(mbpm=1000,pieni=1d-17,zero=0d0,one=1d0,two=2d0)
      dimension z(maxiter),zd(maxiter),tunetest(10),tuneval(10)
    
      ! insertion P.Belanger
      !--------------------
      double precision tune_out
      double complex zw_out
      common/zfunr_out/tune_out(1),zw_out(1)


      
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

    !   write(6,*)'DFFT:',zf,zfd

      dtunea1=dble(zf)*dble(zfd)+dimag(zf)*dimag(zfd)
      num=1

    !   write(6,*)'BEFORE ALL:',tunea1, tunea2, dtunea1, dtunea2
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
        !   write(6,*)'test:',tunetest
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
      tune_out = tune
      zw_out   = zw
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
 