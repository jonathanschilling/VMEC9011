subroutine printout(i, w0, r00)

      use stel_kinds, only: dp
      use name0, only: czero, cp5, c1p0
      use name1, only: nznt, nvac
      use fsqu, only: wp, wb, mns
      use scalars, only: twopi, ns
      use inputdat, only: nfp
      use spectra, only: specw
      use xstuff, only: xc
      use profs, only: phips
      use extfld, only: ivac
      use scalefac, only: wint
      use magfield, only: dbsq, bsqsav
      use fsqu, only: fsqr, fsqz, fsql, fsqr1, fsqz1, fedge
      use time, only: delt

      implicit none

      integer,       intent(in) :: i, j
      real(kind=dp), intent(in) :: w0, r00

      real(kind=dp) :: betav, w, delbsq
      real(kind=dp) :: avm, den


      betav = wp/wb
      w = w0 * twopi*twopi/nfp

      ! avm: volume-averaged spectral width <M>
      avm = czero
      den = czero
      specw(1) = c1p0
      call spectrum(xc, xc(1+mns), xc(1+2*mns), xc(1+3*mns))
      do j = 2, ns
        den = den + phips(j)
        avm = avm + phips(j) * (specw(j) + specw(j-1))
      end do
      avm = cp5*avm/den

      if (ivac.ge.1) then
        delbsq = ddot(nznt, dbsq,        1, wint(2), ns) / &
                 ddot(nznt, bsqsav(1,3), 1, wint(2), ns)
      end if

      if (i.eq.1) then
        if (nvac.eq.0) then
          ! fixed-boundary header
          print   15
          write(3,20)
        else ! nvac != 0
          ! free-boundary header
          print   25
          write(3,30)
        end if
      end if
 15   format(/,' ITER    FSQR      FSQZ      FSQL      fsqr      fsqz     R00(0)     WMHD',/)
 20   format(/,' ITER    FSQR      FSQZ      FSQL      fsqr      fsqz      DELT     R00(0)        WMHD        BETA     <M>',/)
 25   format(/,' ITER    FSQR      FSQZ      FSQL      R00(0)     WMHD      DEL-BSQ',/)
 30   format(/,' ITER    FSQR      FSQZ      FSQL      fsqr      fsqz      DELT     R00(0)        WMHD        BETA     <M>   DEL-BSQ   FEDGE',/)

      if (nvac.eq.0) then
        ! fixed-boundary iteration line
        PRINT   35,I,FSQR,FSQZ,FSQL,FSQR1,FSQZ1,     R00,W                        ! matches format 15
        write(3,50)i,fsqr,fsqz,fsql,fsqr1,fsqz1,delt,r00,w,betav,avm              ! matches format 20
      else
        ! free-boundary iteration line
        PRINT   45,I,FSQR,FSQZ,FSQL,                 R00,W,          DELBSQ       ! matches format 25
        write(3,50)i,fsqr,fsqz,fsql,fsqr1,fsqz1,delt,r00,w,betav,avm,delbsq,fedge ! matches format 30
      end if
 35   format(i5,1p5e10.2,1pe10.3,1p2e11.4)
 45   format(i5,1p3e10.2,1pe10.3,1p2e11.4)
 50   format(i5,1p6e10.2,1pe11.4,1pe15.8,1pe9.2,0pf7.3,1p2e9.2)

      return
end
