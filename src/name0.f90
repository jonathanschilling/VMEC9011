module name0
      use stel_kinds, only: dp
      implicit none

      real(kind=dp), parameter :: zerod  = 0.0e+00_dp
      real(kind=dp), parameter :: oned   = 1.0e+00_dp

      real(kind=dp), parameter :: czero  = zerod
      real(kind=dp), parameter :: c1pm13 = 1.0e-13_dp
      real(kind=dp), parameter :: c1pm10 = 1.0e-10_dp
      real(kind=dp), parameter :: c2pm8  = 2.0e-08_dp
      real(kind=dp), parameter :: cp15   = 1.5e-01_dp
      real(kind=dp), parameter :: cp25   = 2.5e-01_dp
      real(kind=dp), parameter :: cp5    = 5.0e-01_dp
      real(kind=dp), parameter :: cp707d = 7.071e-01_dp
      real(kind=dp), parameter :: cp9    = 9.0e-01_dp
      real(kind=dp), parameter :: cp96   = 9.6e-01_dp
      real(kind=dp), parameter :: c1p0   = oned
      real(kind=dp), parameter :: c1p1   = 1.1e+00_dp
      real(kind=dp), parameter :: c1p4   = 1.4e+00_dp
      real(kind=dp), parameter :: c1p5   = 1.5e+00_dp
      real(kind=dp), parameter :: c2p0   = 2.0e+00_dp
      real(kind=dp), parameter :: c3p0   = 3.0e+00_dp
      real(kind=dp), parameter :: c8p0   = 8.0e+00_dp
      real(kind=dp), parameter :: c100p  = 1.0e+02_dp
end module name0
