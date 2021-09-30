subroutine wrout(bsq, gsqrt, bsubu, bsubv, bsubs, br, bphi, bz, lu, lv)

      use stel_kinds, only: dp
      use name0, only: c2p0
      use name1, only: nznt, mpol, nmax, mnmax
      use scalars, only: ns, itfsq, mns
      use inputdat, only: gam, nfp, niter
      use xstuff, only: xc
      use extfld, only: bscale
      use profs, only: pres

      implicit none

      real(kind=dp), intent(in) :: bsq  (ns,nznt)
      real(kind=dp), intent(in) :: gsqrt(ns,nznt)
      real(kind=dp), intent(in) :: bsubu(ns,nznt)
      real(kind=dp), intent(in) :: bsubv(ns,nznt)
      real(kind=dp), intent(in) :: bsubs(ns,nznt)
      real(kind=dp), intent(in) :: br      (nznt)
      real(kind=dp), intent(in) :: bphi    (nznt)
      real(kind=dp), intent(in) :: bz      (nznt)
      real(kind=dp), intent(in) :: lu   (ns,nznt)
      real(kind=dp), intent(in) :: lv   (ns,nznt)

      real(kind=dp) :: bmod(nznt)
      real(kind=dp) :: phi(nsd)

      real(kind=dp) :: rmnc(mnmax)
      real(kind=dp) :: zmns(mnmax)
      real(kind=dp) :: lmns(mnmax)
      real(kind=dp) :: xm(mnmax)
      real(kind=dp) :: xn(mnmax)

      real(kind=dp) :: bmodmn (mnmax) ! |B| at ns/2
      real(kind=dp) :: bmodmn1(mnmax) ! |B| at ns (LCFS)

      integer       :: js, mn, lk, m, n, nmin0

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

      do js =1, ns

        call convert(rmnc, zmns, lmns, xm, xn, js,          &
                     xc,          xc(1+  mns), xc(1+2*mns), &
                     xc(1+3*mns), xc(1+4*mns), xc(1+5*mns) )

        do lk = 1, nznt
          bmod(lk)=sqrt(c2p0*abs(bsq(js,lk)/bscale**2 - pres(js)))
        end do

        mn = 0
        do m = 0, mpol1

          nmin0 = -nmax
          if (m.eq.0) &
            nmin0 = 0

          do n = nmin0, nmax
            mn = mn+1

            dmult = c2p0*dnorm/(mscale(m)*nscale(iabs(n)))
            if (m.eq.0 .and. n.eq.0) &
              dmult=cp5*dmult

            gmn     = czero
            bmn     = czero
            bsubumn = czero
            bsubvmn = czero
            bsubsmn = czero
            bsupumn = czero
            bsupvmn = czero

            if (js.eq.1) then
              ! first flux surface: interlink mode number arrays with
              ! Fourier coefficient arrays
              if (lwouttxt) then
                WRITE(8,702) NINT(XM(MN)), NINT(XN(MN))
              else
                WRITE(8    ) NINT(XM(MN)), NINT(XN(MN))
              end if

            else ! js.gt.1

              do j = 1,ntheta2
                do k = 1,nzeta
                  lk = k + nzeta*(j-1)

                  if (n.ge.0) then
                    tcosi = dmult*(cosmui(j,m)*cosnv(k, n) + sinmu (j,m)*sinnv(k, n))
                    tsini = dmult*(sinmu (j,m)*cosnv(k, n) - cosmui(j,m)*sinnv(k, n))
                  else
                    tcosi = dmult*(cosmui(j,m)*cosnv(k,-n) - sinmu (j,m)*sinnv(k,-n))
                    tsini = dmult*(sinmu (j,m)*cosnv(k,-n) + cosmui(j,m)*sinnv(k,-n))
                  endif

                  bmn     = bmn     + tcosi*bmod    (lk)
                  gmn     = gmn     + tcosi*gsqrt(js,lk)
                  bsubumn = bsubumn + tcosi*bsubu(js,lk)
                  bsubvmn = bsubvmn + tcosi*bsubv(js,lk)
                  bsubsmn = bsubsmn + tsini*bsubs(js,lk)
                  bsupumn = bsupumn + tcosi*lv   (js,lk)
                  bsupvmn = bsupvmn + tcosi*lu   (js,lk)
                end do
              end do

              if (js.eq.ns/2) &
                bmodmn (mn) = bmn

              if (js.eq.ns) &
                bmodmn1(mn) = bmn

            end if

            if (lwouttxt) then
              WRITE(8,703) RMNC(MN), ZMNS(MN), LMNS(MN), bmn, gmn,         &
                           bsubumn/bscale, bsubvmn/bscale, bsubsmn/bscale, &
                           bsupumn/bscale, bsupvmn/bscale
            else
              WRITE(8    ) RMNC(MN), ZMNS(MN), LMNS(MN), bmn, gmn,         &
                           bsubumn/bscale, bsubvmn/bscale, bsubsmn/bscale, &
                           bsupumn/bscale, bsupvmn/bscale
            end if
          end do ! n
        end do ! m
      end do ! js
 702  format(2i10)
 703  format(5e20.13)

      phi(1) = czero
      do js = 2, ns
        phi(js) = twopi*hs * dsum(js-1, phips(2), 1)
      end do

      fac = abs(bscale)**(gam-c2p0)

      do js = 2, ns
        if (lwouttxt) then
          WRITE(8,703) IOTAS(JS), MASS(JS)*FAC, PRES(JS), PHIPS(JS), BUCO(JS), &
                       bvco(js), phi(js), vp(js), ju(js), jv(js),              &
                       specw(js)
        else
          WRITE(8    ) IOTAS(JS), MASS(JS)*FAC, PRES(JS), PHIPS(JS), BUCO(JS), &
                       bvco(js), phi(js), vp(js), ju(js), jv(js),              &
                       specw(js)
        end if
      end do

      do i = 1, 100
        if (lwouttxt) then
          WRITE(8,703) FSQT(I), WDOT(I)
        else
          WRITE(8    ) FSQT(I), WDOT(I)
        end if
      end do

      if (nvac.eq.0) &
        return

      ! only free-boundary-related outputs below...

      write(3,60) rbtor/bscale, ctor/bscale, bscale
 60   format(/,' 2*PI*R*BTOR = ',1pe16.8, &
               ' NET TOROIDAL CURRENT = ',1pe16.8,&
               ' BSCALE = ',1pe16.8)

      ntskip = 1 + ntheta1/12
      nzskip = 1 + nzeta/6

      write(3,70)
 70   format(/,4x,'ZETA',8x,' Rb ',8x,' Zb ',6x,
     >'BSQMHDI',5x,'BSQVACI',5x,'BSQMHDF',5x,'BSQVACF',/)
      do l = 1, nzeta, nzskip
          zeta = real(360*(l-1))/real(nzeta)

          do k=1,ntheta2,ntskip
            lk=l+nzeta*(k-1)
            nl = ns*lk

            if(iprint.eq.1)write(3,90)zeta,r1(nl)+r1(nl+nrzt),
     >      z1(nl)+z1(nl+nrzt),(bsqsav(lk,n)/bscale**2,n=1,3),
     >      bsqvac(lk)/bscale**2
        end do
      end do
 90   format(1pe10.2,1p6e12.4)

      write(3,75)
 75   format(/,4x,'ZETA',8x,' Rb ',8x,' Zb ',6x, &
                  'BR', 8x,'BPHI', 6x,'BZ',8x,   &
                  'BRv',7x,'BPHIv',5x,'BZv',/)
      do l = 1, nzeta, nzskip
        zeta = real(360*(l-1))/real(nzeta)
        do k=1,ntheta2,ntskip
          lk=l+nzeta*(k-1)
          nl = ns*lk

          write(3,95) zeta, r1(nl)+r1(nl+nrzt), z1(nl)+z1(nl+nrzt), &
                     (1.5*br(nl)   - 0.5*br(nl-1))/bscale,
                     (1.5*bphi(nl) - 0.5*bphi(nl-1))/bscale,
                     (1.5*bz(nl)   - 0.5*bz(nl-1))/bscale, &
                     brv(lk)/bscale, bphiv(lk)/bscale, bzv(lk)/bscale

          end do
        end do
      end do
 95   format(1pe10.2,1p2e12.4,1p6e10.2)

      write(3,100)
 100  format(//,3x,'mb',2x,'nb',9x,'rbc',9x,'zbs',3x,'|B|(s=.5)',
     >3x,'|B|(s=1.)',6x,'mb',2x,'nb',9x,'rbc',9x,'zbs',3x,
     >'|B|(s=.5)',3x,'|B|(s=1.)'/)
      do mn=1,mnmax,2
        write(3,115)nint(xm(mn)),nint(xn(mn)/nfp),rmnc(mn),zmns(mn),
     >  bmodmn(mn),bmodmn1(mn),nint(xm(mn+1)),nint(xn(mn+1)/nfp),
     >  rmnc(mn+1),zmns(mn+1),bmodmn(mn+1),bmodmn1(mn+1)
      end do
 115  format(i5,i4,1p4e12.4,3x,i5,i4,1p4e12.4)

      write(3,120)
 120  format(/,3x,'mf',2x,'nf',5x,'potvacs',6x,'mf',2x,'nf',5x,
     >'potvacs',6x,'mf',2x,'nf',5x,'potvacs'/)
      do mn=1,mpmax,3
        write(3,135)nint(xmpot(mn)),nint(xnpot(mn)/nfp),
     >  potvac(mn)/bscale,nint(xmpot(mn+1)),nint(xnpot(mn+1)/nfp),
     >  potvac(mn+1)/bscale,nint(xmpot(mn+2)),nint(xnpot(mn+2)/nfp),
     >  potvac(mn+2)/bscale
      end do
 135  format(i5,i4,1pe12.4,3x,i5,i4,1pe12.4,3x,i5,i4,1pe12.4)

      return
end
