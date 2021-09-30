subroutine wrout(bsq, gsqrt, bsubu, bsubv, bsubs, br, bphi, bz, lu, lv)

      use stel_kinds, only: dp

      implicit none

      real(kind=dp), intent(in) :: bsq  (ns,*)
      real(kind=dp), intent(in) :: gsqrt(ns,*)
      real(kind=dp), intent(in) :: bsubu(ns,*)
      real(kind=dp), intent(in) :: bsubv(ns,*)
      real(kind=dp), intent(in) :: bsubs(ns,*)
      real(kind=dp), intent(in) :: br      (*)
      real(kind=dp), intent(in) :: bphi    (*)
      real(kind=dp), intent(in) :: bz      (*)
      real(kind=dp), intent(in) :: lu   (ns,*)
      real(kind=dp), intent(in) :: lv   (ns,*)

      real(kind=dp) :: bmod(nznt)
      real(kind=dp) :: phi(nsd)

      real(kind=dp) :: rmnc(mnmax)
      real(kind=dp) :: zmns(mnmax)
      real(kind=dp) :: lmns(mnmax)
      real(kind=dp) :: xm(mnmax)
      real(kind=dp) :: xn(mnmax)

      real(kind=dp) :: bmodmn(mnmax)
      real(kind=dp) :: bmodmn1(mnmax)

      ! set lwouttxt to false to get a binary ("unformatted") wout file
      logical, parameter :: lwouttxt = .true.

      ! THIS SUBROUTINE CREATES THE TEXT FILE 'WOUT'
      ! WHICH CONTAINS THE CYLINDRICAL COORDINATE
      ! SPECTRAL COEFFICIENTS RMNC,ZMNS,LMNS (FULL MESH)

      if (lwouttxt) then
        WRITE(8,701)GAM,NFP,NS,MPOL,NMAX,MNMAX,ITFSQ,NITER/100+1
      else
        WRITE(8    )GAM,NFP,NS,MPOL,NMAX,MNMAX,ITFSQ,NITER/100+1
      end if
  701 FORMAT(F10.3,7I6)

      do 20 js=1,ns

        call convert(rmnc,zmns,lmns,xm,xn,js,xc,xc(1+mns),
                  >  xc(1+2*mns),xc(1+3*mns),xc(1+4*mns),xc(1+5*mns))

        mn = 0
        do 25 lk=1,nznt
 25     bmod(lk)=sqrt(c2p0*abs(bsq(js,lk)/bscale**2-pres(js)))

        do 20 m = 0,mpol1
          nmin0 = -nmax
          if(m.eq.0)nmin0 = 0
            do 20 n = nmin0,nmax
              mn = mn+1
              dmult = c2p0*dnorm/(mscale(m)*nscale(iabs(n)))
              if(m.eq.0.and.n.eq.0)dmult=cp5*dmult
              gmn = czero
              bmn = czero
              bsubumn = czero
              bsubvmn = czero
              bsubsmn = czero
              bsupumn = czero
              bsupvmn = czero

              if (js.eq.1) &
                goto 22

              do 200 j = 1,ntheta2
                do 200 k = 1,nzeta
                  lk = k + nzeta*(j-1)
                  if(n.ge.0)then
                  tcosi = dmult*(cosmui(j,m)*cosnv(k,n)+sinmu(j,m) *sinnv(k,n))
                  tsini = dmult*(sinmu (j,m)*cosnv(k,n)-cosmui(j,m)*sinnv(k,n))
                  else
                  tcosi = dmult*(cosmui(j,m)*cosnv(k,-n)-sinmu(j,m) *sinnv(k,-n))
                  tsini = dmult*(sinmu (j,m)*cosnv(k,-n)+cosmui(j,m)*sinnv(k,-n))
                  endif
                  bmn = bmn + tcosi*bmod(lk)
                  gmn = gmn + tcosi*gsqrt(js,lk)
                  bsubumn = bsubumn + tcosi*bsubu(js,lk)
                  bsubvmn = bsubvmn + tcosi*bsubv(js,lk)
                  bsubsmn = bsubsmn + tsini*bsubs(js,lk)
                  bsupumn = bsupumn + tcosi*lv(js,lk)
                  bsupvmn = bsupvmn + tcosi*lu(js,lk)
 200              continue

              if(js.eq.ns/2)bmodmn (mn) = bmn
              if(js.eq.ns)  bmodmn1(mn) = bmn

              goto 20

 22           if (lwouttxt) then
                WRITE(8,702)NINT(XM(MN)),NINT(XN(MN))
              else
                WRITE(8    )NINT(XM(MN)),NINT(XN(MN))
              end if

 20           if (lwouttxt) then
                WRITE(8,703)RMNC(MN),ZMNS(MN),LMNS(MN),
     >          bmn,gmn,bsubumn/bscale,bsubvmn/bscale,bsubsmn/bscale,
     >          bsupumn/bscale,bsupvmn/bscale
              else
                WRITE(8    )RMNC(MN),ZMNS(MN),LMNS(MN),
     >          bmn,gmn,bsubumn/bscale,bsubvmn/bscale,bsubsmn/bscale,
     >          bsupumn/bscale,bsupvmn/bscale
              end if
      end do
 702  format(2i10)
 703  format(5e20.13)

      phi(1) = czero
      do 30 js = 2,ns
 30     phi(js) = twopi*hs*ssum(js-1,phips(2),1)
      end do

      fac = abs(bscale)**(gam-c2p0)

      if (lwouttxt) then
        WRITE(8,703)(IOTAS(JS),MASS(JS)*FAC,PRES(JS),PHIPS(JS),BUCO(JS),
     >  bvco(js),phi(js),vp(js),ju(js),jv(js),specw(js),js=2,ns)
      else
        WRITE(8    )(IOTAS(JS),MASS(JS)*FAC,PRES(JS),PHIPS(JS),BUCO(JS),
     >  bvco(js),phi(js),vp(js),ju(js),jv(js),specw(js),js=2,ns)
      end if

      if (lwouttxt) then
        WRITE(8,703)(FSQT(I),WDOT(I),I=1,100)
      else
        WRITE(8    )(FSQT(I),WDOT(I),I=1,100)
      end if

      if (nvac.eq.0) &
        return

      ! only free-boundary-related outputs below...

      write(3,60)rbtor/bscale,ctor/bscale,bscale
 60   format(/,' 2*PI*R*BTOR = ',1pe16.8,' NET TOROIDAL CURRENT = ',
     >1pe16.8,' BSCALE = ',1pe16.8)
      do 80 iprint=1,2
      if(iprint.eq.1)write(3,70)
      if(iprint.eq.2)write(3,75)
      ntskip=1+ntheta1/12
      nzskip=1+nzeta/6
      do 80 l=1,nzeta,nzskip
      zeta = real(360*(l-1))/real(nzeta)
      do 80 k=1,ntheta2,ntskip
      lk=l+nzeta*(k-1)
      nl = ns*lk
      if(iprint.eq.1)write(3,90)zeta,r1(nl)+r1(nl+nrzt),
     >z1(nl)+z1(nl+nrzt),(bsqsav(lk,n)/bscale**2,n=1,3),
     >bsqvac(lk)/bscale**2
      if(iprint.eq.2)write(3,95)zeta,r1(nl)+r1(nl+nrzt),
     >z1(nl)+z1(nl+nrzt),(1.5*br(nl) - 0.5*br(nl-1))/bscale,
     >(1.5*bphi(nl) - 0.5*bphi(nl-1))/bscale,
     >(1.5*bz(nl) - 0.5*bz(nl-1))/bscale,brv(lk)/bscale,
     >bphiv(lk)/bscale,bzv(lk)/bscale
 80   continue
 70   format(/,4x,'ZETA',8x,' Rb ',8x,' Zb ',6x,
     >'BSQMHDI',5x,'BSQVACI',5x,'BSQMHDF',5x,'BSQVACF',/)
 75   format(/,4x,'ZETA',8x,' Rb ',8x,' Zb ',6x,
     >'BR',8x,'BPHI',6x,'BZ',8x,'BRv',7x,'BPHIv',5x,'BZv',/)
 90   format(1pe10.2,1p6e12.4)
 95   format(1pe10.2,1p2e12.4,1p6e10.2)
      write(3,100)
 100  format(//,3x,'mb',2x,'nb',9x,'rbc',9x,'zbs',3x,'|B|(s=.5)',
     >3x,'|B|(s=1.)',6x,'mb',2x,'nb',9x,'rbc',9x,'zbs',3x,
     >'|B|(s=.5)',3x,'|B|(s=1.)'/)
      do 110 mn=1,mnmax,2
 110  write(3,115)nint(xm(mn)),nint(xn(mn)/nfp),rmnc(mn),zmns(mn),
     >bmodmn(mn),bmodmn1(mn),nint(xm(mn+1)),nint(xn(mn+1)/nfp),
     >rmnc(mn+1),zmns(mn+1),bmodmn(mn+1),bmodmn1(mn+1)
 115  format(i5,i4,1p4e12.4,3x,i5,i4,1p4e12.4)
      write(3,120)
 120  format(/,3x,'mf',2x,'nf',5x,'potvacs',6x,'mf',2x,'nf',5x,
     >'potvacs',6x,'mf',2x,'nf',5x,'potvacs'/)
      do 130 mn=1,mpmax,3
 130  write(3,135)nint(xmpot(mn)),nint(xnpot(mn)/nfp),
     >potvac(mn)/bscale,nint(xmpot(mn+1)),nint(xnpot(mn+1)/nfp),
     >potvac(mn+1)/bscale,nint(xmpot(mn+2)),nint(xnpot(mn+2)/nfp),
     >potvac(mn+2)/bscale
 135  format(i5,i4,1pe12.4,3x,i5,i4,1pe12.4,3x,i5,i4,1pe12.4)

      return
end
