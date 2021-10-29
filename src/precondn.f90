subroutine precondn(lu, bsq, gsqrt, r12, wint, &
                    xs, xu12, xue, xuo, xodd,  &
                    axm, axd, axp, bxm, bxd, cx)

      use stel_kinds, only: dp
      use name0, only: czero, cp25, cp5, c1p0, c1p5
      use name1, only: nznt, nsd1
      use scalars, only: iter2, ns, ohs
      use profs, only: pres
      use precond, only: sm, sp
      use scalefac, only: shalf

      implicit none

      ! general inputs
      real(kind=dp), intent(in)  :: lu   (ns,nznt) ! B^zeta = phip/sqrtg*[1 + d(lambda)/du]
      real(kind=dp), intent(in)  :: bsq  (ns,nznt) ! B^2/(2 \mu_0) + p
      real(kind=dp), intent(in)  :: gsqrt(ns,nznt) ! Jacobian on half-grid
      real(kind=dp), intent(in)  :: r12  (ns,nznt) ! R on half-grid
      real(kind=dp), intent(in)  :: wint (ns,nznt) ! weighting factor for surface integrals

      ! R/Z-specific inputs
      !                     force to precondition:    R   |  Z
      real(kind=dp), intent(in)  :: xs   (ns,nznt) !  Z_s |  R_s          half-grid
      real(kind=dp), intent(in)  :: xu12 (ns,nznt) !  Z_u |  R_u          half-grid
      real(kind=dp), intent(in)  :: xue  (ns,nznt) !  Z_u |  R_u (even-m) full-grid
      real(kind=dp), intent(in)  :: xuo  (ns,nznt) !  Z_u |  R_u (odd-m)  full-grid
      real(kind=dp), intent(in)  :: xodd (ns,nznt) !  Z   |  R   (odd-m)  full-grid

      ! preconditioning matrix elements outputs
      real(kind=dp), intent(out) :: axm  (nsd1,2)  ! aRm  | aZm
      real(kind=dp), intent(out) :: axd  (nsd1,2)  ! aRd  | aZd
      real(kind=dp), intent(out) :: axp  (nsd1,2)  ! aRp  | aZp
      real(kind=dp), intent(out) :: bxm  (nsd1,2)  ! bRm  | bZm
      real(kind=dp), intent(out) :: bxd  (nsd1,2)  ! bRd  | bZd
      real(kind=dp), intent(out) :: cx   (nsd1)    ! cR   | cZ

      ! internal temporary storage
      real(kind=dp) :: ax(nsd1,4)
      real(kind=dp) :: bx(nsd1,4) ! why 4? only 3 are used... (maybe to ease init loop?)
      real(kind=dp) :: ptau

      integer       :: i, js, lk
      real(kind=dp) :: t1, t2, t3
      real(kind=dp) :: axm_js_2, axp_js_2, axd_js_3, axd_js_4, D, r,g, b

      if (iter2.le.1) then
        ! setup interpolation magic
        do js = 2, ns
          sm(js) = sqrt( (js - c1p5)/(js - c1p0) ) ! sqrt(s_{j-1/2})/sqrt(s_j)
          sp(js) = sqrt( (js -  cp5)/(js - c1p0) ) ! sqrt(s_{j+1/2})/sqrt(s_j)
        end do
        sm(1) = czero
        sp(0) = czero
        sp(1) = sm(2)
      endif

      ! COMPUTE PRECONDITIONING MATRIX ELEMENTS FOR R,Z FORCE (ALL ARE MULTIPLIED BY 0.5).

      ! initialize ax, bx, cx to zero
      do i = 1, 4
        do js = 1, ns+1
          ax(js, i) = czero
          bx(js, i) = czero
        end do
      end do
      do js = 1, ns+1
         cx(js) = czero
      end do

      ! compute matrix elements on half-grid
      do js = 2, ns

        ! COMPUTE DOMINANT (1/DELTA-S)**2 PRECONDITIONING MATRIX ELEMENTS
        do lk = 1, nznt
          ptau = r12(js,lk)**2 * (bsq(js,lk) - pres(js)) * wint(js,lk)/gsqrt(js,lk)

          t3 = ohs  *  xu12(js, lk)

          ax(js,1) = ax(js,1) + ptau * t3 * t3

          t1 = cp5 * (xs(js,lk) + cp5*xodd(js,  lk)/shalf(js))
          t2 = cp5 * (xs(js,lk) + cp5*xodd(js-1,lk)/shalf(js))

          bx(js,1) = bx(js,1) + ptau*t1*t2
          bx(js,2) = bx(js,2) + ptau*t1**2
          bx(js,3) = bx(js,3) + ptau*t2**2

          cx(js) = cx(js) + cp25 * lu(js,lk)**2 * gsqrt(js,lk)*wint(js,lk)
        end do
      end do

      ! radial interpolation onto some other mesh ???
      ! averaging of neighboring half-grid points onto full grid
      do js = 1, ns
        ! even-m
        axm(js,1) =-ax(js,1)
        axd(js,1) = ax(js,1) + ax(js+1,1)
        axp(js,1) =          - ax(js+1,1)

        ! odd-m
        axm_js_2 = 0.0_dp
        axp_js_2 = 0.0_dp
        axd_js_3 = 0.0_dp
        axd_js_4 = 0.0_dp

        if (js .ge. 2) then

          r = 0.0_dp
          g = 0.0_dp

          do lk = 1, nznt
            ptau = r12(js,lk)**2 * (bsq(js,lk) - pres(js)) * wint(js,lk)/gsqrt(js,lk)

            D = ptau * xu12(js,lk)**2

            r = r + shalf(js)**2 * D * ohs**2
            g = g + shalf(js)**2 * D * ohs**2

            r = r + cp25*( ohs*(xue(js,lk)-xue(js-1,lk)) + shalf(js)*ohs*(xuo(js,lk)-xuo(js-1,lk)) ) * ptau *xu12(js,lk)
            g = g + cp5*ohs*xue(js,lk)*ptau*xu12(js,lk)
            g = g + cp5*ohs*xuo(js,lk)*shalf(js)*ptau*xu12(js,lk)

            r = r + xue(js-1,lk)*xue(js,lk)*ptau/(16.0_dp * shalf(js)**2)
            g = g + xue(js  ,lk)*xue(js,lk)*ptau/(16.0_dp * shalf(js)**2)

            r = r + ( xue(js-1,lk)*xuo(js,lk) + xuo(js-1,lk)*xue(js,lk) ) * ptau/(16.0_dp * shalf(js))
            g = g + xue(js,lk)*xuo(js,lk) * ptau/(8.0_dp * shalf(js))

            r = r + xuo(js-1,lk)*xuo(js,lk)*ptau/16.0_dp
            g = g + xuo(js  ,lk)*xuo(js,lk)*ptau/16.0_dp

          end do

          axm_js_2 = - r/shalf(js)**2
          axd_js_3 =   g/shalf(js)**2

        end if

        if (js .le. ns) then

          g = 0.0_dp
          b = 0.0_dp

          do lk = 1, nznt
            ptau = r12(js+1,lk)**2 * (bsq(js+1,lk) - pres(js+1)) * wint(js+1,lk)/gsqrt(js+1,lk)

            D = ptau * xu12(js+1,lk)**2

            b = b + shalf(js+1)**2 * D * ohs**2
            g = g + shalf(js+1)**2 * D * ohs**2

            b = b + cp25*( ohs*(xue(js+1,lk)-xue(js,lk)) + shalf(js+1)*ohs*(xuo(js+1,lk)-xuo(js,lk)) ) * ptau *xu12(js+1,lk)
            g = g - cp5*ohs*xue(js,lk)*ptau*xu12(js+1,lk)
            g = g - cp5*ohs*xuo(js,lk)*shalf(js)*ptau*xu12(js+1,lk)

            b = b - xue(js+1,lk)*xue(js,lk)*ptau/(16.0_dp * shalf(js+1)**2)
            g = g + xue(js  ,lk)*xue(js,lk)*ptau/(16.0_dp * shalf(js+1)**2)

            b = b - ( xue(js+1,lk)*xuo(js,lk) + xuo(js+1,lk)*xue(js,lk) ) * ptau/(16.0_dp * shalf(js+1))
            g = g + xue(js,lk)*xuo(js,lk) * ptau/(8.0_dp * shalf(js+1))

            b = b - xuo(js+1,lk)*xuo(js,lk)*ptau/16.0_dp
            g = g + xuo(js  ,lk)*xuo(js,lk)*ptau/16.0_dp

          end do

          axd_js_4 =  g/shalf(js+1)**2
          axp_js_2 = -b/shalf(js+1)**2
        end if

        axm(js,2) = axm_js_2 * sm(js) * sp(js-1)
        axp(js,2) = axp_js_2 * sm(js+1) * sp(js)
        axd(js,2) = axd_js_3 * sm(js) * sm(js) + axd_js_4 * sp(js) * sp(js)

        bxm(js,1) = bx(js,1)                                              ! off-diagonal, even-m, poloidal
        bxm(js,2) = bx(js,1) * sm(js) * sp(js-1)                          ! off-diagonal,  odd-m, poloidal
        bxd(js,1) = bx(js,2)                     + bx(js+1,3)             !     diagonal, even-m, poloidal
        bxd(js,2) = bx(js,2) * sm(js)**2         + bx(js+1,3) * sp(js)**2 !     diagonal,  odd-m, poloidal

        cx (js)   = cx(js)                       + cx(js+1)               !     diagonal,         toroidal
      end do

      return
end


