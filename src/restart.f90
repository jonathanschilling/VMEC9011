subroutine restart

      use name0,   only: czero, &
                         cp9,   &
                         cp96
      use scalars, only: irst,  &
                         neqs
      use xstuff,  only: xc,    &
                         xcdot, &
                         xstore
      use time,    only: delt

      implicit none

      integer :: l

      if (irst.eq.1) then
        ! copy current state vector into backup
        call scopy(neqs, xc, 1, xstore, 1)
      else
        ! zero velocity, restore state from backup
        do l = 1, neqs
          xcdot(l) = czero
          xc(l)    = xstore(l)
        enddo

        ! reduce time step;
        ! will arrive here with irst.eq.2 or irst.eq.3
        ! irst.eq.2 --> (irst-2) = 0, (3-irst) = 1
        !           --> delt *= 0.9
        ! irst.eq.3 --> (irst-2) = 3, (3-irst) = 0
        !           --> delt *= 0.96
        delt = delt*(cp96*(irst-2) + cp9*(3-irst))

        ! reset to normal state
        irst = 1

      endif

      return

end
