subroutine vsetup

      use name0,    only: czero, &
                          c1p0
      use name1,    only: nrztd
      use realsp,   only: rcon0, &
                          zcon0, &
                         setup_realsp
      use rforces, only: setup_rforces
      use time,     only: timer
      use scalars,  only: iequi, &
                          itfsq
      use extfld,   only: bscale, &
                          ivac
      use magfield, only: delbsq


      implicit none

      integer :: i
      integer :: l

      call setup_realsp
      call setup_rforces

      do l = 1, nrztd
        rcon0(l) = czero
        zcon0(l) = czero
      enddo

      do i = 0, 10
        timer(i) = czero
      enddo

      iequi  = 0
      itfsq  = 0
      bscale = c1p0
      delbsq = c1p0
      ivac   = 0

      return
end
