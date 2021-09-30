subroutine eqsolve(ns, intflag, ierflag)

      use stel_kinds, only: dp
      use name0,      only: c1pm13, &
                            c2pm8,  &
                            c1p0,   &
                            c100p
      use name1,      only: mnd,    &
                            nznt,   &
                            mnmax,  &
                            nsd
      use scalars,    only: hs,     &
                            ohs,    &
                            mns,    &
                            neqs,   &
                            nrzt,   &
                            iter1,  &
                            iter2,  &
                            irst,   &
                            ijacob, &
                            itfsq
      use inputdat,   only: niter,  &
                            gam,    &
                            ftol
      use time,       only: delt,   &
                            timer
      use xstuff,     only: xc
      use mnarray,    only: mscale, &
                            nscale
      use fsqu,       only: wb,     &
                            wp,     &
                            fsq,    &
                            fsqr,   &
                            fsqz,   &
                            fsql,   &
                            fsqt,   &
                            wdot
      use extfld,     only: ivac,   &
                            bscale

      implicit none

      integer, intent(in)  :: ns
      integer, intent(in)  :: intflag
      integer, intent(out) :: ierflag

      integer       :: miter
      real(kind=dp) :: w1   = 0.0_dp
      real(kind=dp) :: r01  = 0.0_dp
      real(kind=dp) :: res1 = 1.0_dp
      real(kind=dp) :: r00, w0
      real(kind=dp) :: r0dot, wdota
      real(kind=dp) :: res0

      !         INDEX OF LOCAL VARIABLES
      !
      ! hs      radial mesh size increment
      ! iequi   counter used to call EQFOR at end of run
      ! ijacob  counter for number of times jacobian changes sign
      ! irst    counter monitoring sign of jacobian; resets R, Z, and
      !         Lambda when jacobian changes sign and decreases time s
      ! isigng  sign of Jacobian : must be =1 (right-handed) or =-1 (left-handed)
      ! iterj   stores position in main iteration loop (j=1,2)
      ! itfsq   counter for storing FSQ into FSQT for plotting
      ! ivac    counts number of free-boundary iterations
      ! ndamp   number of iterations over which damping is averaged
      ! meven   parity selection label for even poloidal modes of R an
      ! modd    parity selection label for odd poloidal modes of R and
      ! rb      boundary coefficient array for R (rcc,rss)
      ! zb      boundary coefficient array for Z (zcs,zsc)
      ! gc      stacked array of R, Z, Lambda Spectral force coefficie
      ! xc      stacked array of scaled R, Z, Lambda Fourier coefficie
      ! STACKING ORDER:
      !   1:mns   => rmncc Fourier coefficients
      !   1+mns,2*mns => rmnss Fourier coefficients
      !   1+2*mns,3*mns => zmncs Fourier coefficients
      !   1+3*mns,4*mns => zmnsc Fourier coefficients
      !   1+4*mns,5*mns => lmncs Fourier coefficients
      !   1+5*mns,neqs => lmnsc Fourier coefficients

      ! INITIALIZE MESH-DEPENDENT SCALARS
      hs    = c1p0/real(ns-1)
      ohs   = c1p0/hs
      mns   = ns*mnd
      neqs  = 6*mns
      nrzt  = nznt*ns
      iter2 = 1
      irst  = 1

      ! sgg
      ! DELT = C1P1
      ! DELT=0.9
      DELT=0.6

      ijacob = -1
      print   *, 'NS = ',ns,' NO. FOURIER MODES = ',mnmax
      write(3,*) 'NS = ',ns,' NO. FOURIER MODES = ',mnmax

      ! COMPUTE INITIAL R, Z AND MAGNETIC FLUX PROFILES
      ! AND STORE XC, XCDOT FOR POSSIBLE RESTART

      call profil3d(xc, xc(1+2*mns), intflag)

      ! TODO: remove need for GOTO 10 below
 10   call restart

      ijacob = ijacob+1
      iter1  = iter2

      ! FORCE ITERATION LOOP
      do miter = iter1, niter

        ! offset in number of iterations
        iter2 = miter

        ! ADVANCE FOURIER AMPLITUDES OF R, Z, AND LAMBDA
        call evolve(ierflag)
        if (ierflag.ne.0) then
          ! evolve failed --> no need to try time step --> return
          return
        end if

        r00 = xc(1)*mscale(0)*nscale(0)
        w0  = wb + wp/(gam-c1p0) ! total energy
        r0dot = abs(r00 - r01)/r00
        wdota = abs(w0 - w1)/w0
        r01 = r00
        w1  = w0

        ! COMPUTE ABSOLUTE STOPPING CRITERION
        if (     (      fsqr.le.ftol               & ! { ftol satisfied for R
                  .and. fsqz.le.ftol               & !   ftol satisfied for Z
                  .and. fsql.le.ftol               & !   ftol satisfied for lambda
                  .and. miter-iter1.gt.10)         & !   at least 10 iterations      }
            .or. (      ns.lt.nsd                  & ! { not at final multi-grid step
                  .and. (fsqr+fsqz+fsql).le.c2pm8) & !   total force residual < 2e-8 }
            .or. ijacob.ge.100                     & ! # of Jacobian resets exceeded
            .or. miter.eq.niter ) then               ! niter iterations exceeded

          ! leave force iterations loop
          break

        endif

        if (ivac.eq.1) then
          print   110, miter, c1p0/bscale
          write(3,110) miter, c1p0/bscale
 110  format(/,'  VACUUM PRESSURE TURNED ON AT ',i4, ' ITERATIONS',/, &
               '  MAGNETIC FIELD IN PLASMA IS SCALED AT END BY ',1pe10.2,' (VAC. UNITS)',/)

          ivac = ivac + 1
        endif

        ! when using fine grid, update fsqt and wdot every 100 iterations
        if ( mod(miter, (niter/100 + 1)).eq.0 .and. (ns.eq.nsd) ) then
          itfsq       = itfsq+1
          fsqt(itfsq) = fsqr + fsqz
          wdot(itfsq) = max(wdota, c1pm13)
        end if

        ! TIME STEP CONTROL
        if (miter.eq.iter1) &
          res0 = fsq
        res0 = min(res0, fsq)

        if ( (fsq.le.res0) .and. ((iter2-iter1).gt.nsd/2) ) then
          ! forces reduced over some iterations --> save state
          call restart
        else if (fsq.gt.c100p*res1 .and. miter.gt.iter1) then
          ! forces much worse than before --> reset to saved state
          irst = 2
        else if (      fsq.gt.(c1p4*res1)      &
                 .and. mod(miter-iter1,5).eq.0 &
                 .and. (miter-iter1).gt.2*ns4  ) then
          ! bad convergence --> hard reset (?)
          irst = 3
        endif

        if (irst.ne.1) then
          ! convergence problem; need to restart from saved state
          goto 10 ! TODO: replace goto with second loop
        end if

        ! if we came here, save current value of total residual force
        res1 = res0

        ! PRINTOUT EVERY NSTEP ITERATIONS or in first iteration
        if (mod(miter,nstep).eq.0 .or. miter.eq.1) &
          call printout(iter2, w1, r00)

      enddo ! force iterations loop

      ! final printout
      call printout(iter2, w0, r00)

      ! note runtime
      call second(timer(0))

      if (ijacob.ge.100) then
        ! error: more than 100 Jacobian resets
        ierflag = 2
      endif

      if (ns.lt.nsd) then
        print   100, timer(0), ijacob
        write(3,100) timer(0), ijacob
      end if
 100  format(/,' COMPUTATIONAL TIME FOR INITIAL INTERPOLATION = ',1pe10.2,' SECONDS ', &
               ' JACOBIAN RESETS = ',i4,/)

      print   60, wdota, r0dot
      write(3,60) wdota, r0dot
 60   format(/,' d(ln W)/dt = ',1pe10.3,' d(ln R0)/dt = ',1pe10.3,/)

      return

end
