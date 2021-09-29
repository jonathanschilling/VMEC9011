subroutine eqfor(bu, bv, bsq, rmag, zmag)

      use stel_kinds, only: dp

      implicit none

      real(kind=dp), intent(in) :: bu(ns,1)
      real(kind=dp), intent(in) :: bv(ns,1)
      real(kind=dp), intent(in) :: bsq(ns,1)
      real(kind=dp), intent(in) :: rmag(ns,0:nmax)
      real(kind=dp), intent(in) :: zmag(ns,0:nmax)

      write(3,5)
  5   format(/,'     S       EQUIF      PHI       CHI     <JTHETA>     IOTA     <JZETA>   <BSUBV>    P''/V''    V''/PHIP   <M>',/)

      betaxis =  c1p5*pres(2)/(dnorm * ddot(nznt, bsq(2,1), ns, wint(2), ns) - pres(2)) &
               -  cp5*pres(3)/(dnorm * ddot(nznt, bsq(3,1), ns, wint(2), ns) - pres(3))

      phi1 = czero
      chi1 = czero
      do i = 2,ns
        buco(i)  = dnorm/bscale * ddot(nznt, bu(i,1), ns, wint(2), ns)
        bvco(i)  = dnorm/bscale * ddot(nznt, bv(i,1), ns, wint(2), ns)

        pres(i)  = pres(i)/bscale**2
        phips(i) = phips(i)/bscale
      end do

      do js = 2, ns-1
        es = (js-1)/real(ns-1)

        t0     = cp5*(   vp(js) +    vp(js+1))
        aiotaf = cp5*(iotas(js) + iotas(js+1))

        jv(js) =  isigng*ohs*(buco(js+1) - buco(js))/t0
        ju(js) = -isigng*ohs*(bvco(js+1) - bvco(js))/t0

        t1 = jv(js)*aiotaf
        t2 = ohs*(pres(js+1) - pres(js))/t0
        t3 = cp5*(vp(js+1)/phips(js+1) + vp(js)/phips(js))

        phi1 = phi1 + phips(js)          /ohs
        chi1 = chi1 + phips(js)*iotas(js)/ohs

        equif(js) = (ju(js) - t1 - t2*t3) / (abs(ju(js)) + abs(t1) + abs(t2*t3))

        write(3,30) es, equif(js), &
                    isigng*twopi*phi1, &
                    isigng*twopi*chi1, &
                    ju(js), aiotaf, jv(js), bvco(js), t2, t3, specw(js)
      end do
 30   FORMAT(1P5E10.2,1p1E11.3,1P3E10.2,1P1E11.3,0PF7.3)

      volf = (twopi**2)*dsum(ns, vp, 1)/real(ns-1)
      write(3,40) voli, volf, betaxis
 40   format(/,'  INITIAL VOLUME  =  ',1pe20.10,'  FINAL VOLUME  =  ',1pe20.10,/, &
               '  BETA ON AXIS    =  ',2x,1pe12.4,/)

      write(3,50)
 50   format(2x,'MAGNETIC AXIS COEFFICIENTS',/,&
                '    n     rmag       zmag',/)
      do n = 0, nmax
        t1 = mscale(0)*nscale(n)

        ! TODO: why is the z coeff inverted ?
        write(3,70) n, t1*rmag(1,n), -t1*zmag(1,n)
      end do
 70   format(i5,1p2e12.4)

      return
end
