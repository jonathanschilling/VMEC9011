module name2
      use stel_kinds, only: dp
      implicit none
      real(kind=dp) :: iotas
      real(kind=dp) :: mass
      real(kind=dp) :: mscale
      real(kind=dp) :: nscale
      real(kind=dp) :: ju
      real(kind=dp) :: jv
end module name2

module current
      use stel_kinds, only: dp
      implicit none
      real(kind=dp) :: buco(nsd)
      real(kind=dp) :: bvco(nsd)
      real(kind=dp) :: ju(nsd)
      real(kind=dp) :: jv(nsd)
end module current

module extfld
      use stel_kinds, only: dp
      implicit none
      real(kind=dp) :: potvac(nznt)
      real(kind=dp) :: xmpot(nznt)
      real(kind=dp) :: xnpot(nznt)
      real(kind=dp) :: brv(nznt)
      real(kind=dp) :: bphiv(nznt)
      real(kind=dp) :: bzv(nznt)
      real(kind=dp) :: bscale
      integer       :: mpmax
      integer       :: ivac
      integer       :: ivac2
end module extfld

module fsqu
      use stel_kinds, only: dp
      implicit none
      real(kind=dp) :: fnorm
      real(kind=dp) :: fsqr
      real(kind=dp) :: fsqz
      real(kind=dp) :: fsql
      real(kind=dp) :: fnorm1
      real(kind=dp) :: fsqr1
      real(kind=dp) :: fsqz1
      real(kind=dp) :: fsql1
      real(kind=dp) :: fsq
      real(kind=dp) :: fedge
      real(kind=dp) :: wb
      real(kind=dp) :: wp
      real(kind=dp) :: fsqt(100)
      real(kind=dp) :: wdot(100)
      real(kind=dp) :: equif(nsd)
end module fsqu

module inputdat
      use stel_kinds, only: dp
      implicit none
      real(kind=dp) :: am(0:5)
      real(kind=dp) :: ai(0:5)
      real(kind=dp) :: raxis(0:nmax)
      real(kind=dp) :: zaxis(0:nmax)
      real(kind=dp) :: rb(0:nmax,0:mpol1,2)
      real(kind=dp) :: zb(0:nmax,0:mpol1,2)
      real(kind=dp) :: ftol
      real(kind=dp) :: gam
      real(kind=dp) :: ac(5)
      integer       :: ncurr
      integer       :: nfp
      integer       :: nstep
      integer       :: niter
      integer       :: nvacskip
end module inputdat

module mnarray
      use stel_kinds, only: dp
      implicit none
      real(kind=dp) :: xrz3(0:nmax,0:mpol1)
      real(kind=dp) :: xrz4(0:nmax,0:mpol1)
      real(kind=dp) :: xmpq(0:mpol1,3)
      integer       :: mscale(0:mpol1)
      integer       :: nscale(0:nmax1)
      integer       :: jmin1(0:mpol)
      integer       :: jmin2(0:mpol)
      integer       :: jlam(0:mpol)
end module mnarray

module profs
      use stel_kinds, only: dp
      implicit none
      real(kind=dp) :: iotas(nsd)
      real(kind=dp) :: mass(nsd)
      real(kind=dp) :: phips(nsd1)
      real(kind=dp) :: pres(nsd)
      real(kind=dp) :: vp(nsd)
      real(kind=dp) :: phip(nrztd)
end module profs

module scalars
      use stel_kinds, only: dp
      implicit none
      real(kind=dp) :: dnorm
      real(kind=dp) :: hs
      real(kind=dp) :: ohs
      real(kind=dp) :: twopi
      real(kind=dp) :: voli
      integer       :: ijacob
      integer       :: itfsq
      integer       :: iequi
      integer       :: irst
      integer       :: iter1
      integer       :: iter2
      integer       :: isigng
      integer       :: meven
      integer       :: modd
      integer       :: ndamp
      integer       :: ns
      integer       :: ns4
      integer       :: neqs
      integer       :: nrzt
      integer       :: mns
end module scalars

module scalefac
      use stel_kinds, only: dp
      implicit none
      real(kind=dp) :: faclam(mnd2*nsd)
      real(kind=dp) :: shalf(nrztd)
      real(kind=dp) :: sqrts(nrztd)
      real(kind=dp) :: scalxc(neq)
      real(kind=dp) :: wint(nrztd)
end module scalefac

module spectra
      use stel_kinds, only: dp
      implicit none
      real(kind=dp) :: faccon(0:nmax,0:mpol1)
      real(kind=dp) :: specw(nsd)
      real(kind=dp) :: tcon(nsd)
end module spectra

module time
      use stel_kinds, only: dp
      implicit none
      real(kind=dp) :: delt
      real(kind=dp) :: otav
      real(kind=dp) :: otau(15)
      real(kind=dp) :: timer(0:10)
end module time

module trignew
      use stel_kinds, only: dp
      implicit none
      real(kind=dp) :: cosmu  (ntheta2,0:mpol1)
      real(kind=dp) :: sinmu  (ntheta2,0:mpol1)
      real(kind=dp) :: cosmum (ntheta2,0:mpol1)
      real(kind=dp) :: sinmum (ntheta2,0:mpol1)
      real(kind=dp) :: cosmui (ntheta2,0:mpol1)
      real(kind=dp) :: cosmumi(ntheta2,0:mpol1)
      real(kind=dp) :: cosnv  (nzeta,  0:nmax)
      real(kind=dp) :: sinnv  (nzeta,  0:nmax)
      real(kind=dp) :: cosnvn (nzeta,  0:nmax)
      real(kind=dp) :: sinnvn (nzeta,  0:nmax)
end module trignew

module magfield
      use stel_kinds, only: dp
      implicit none
      real(kind=dp) :: bsqsav(nznt,3)
      real(kind=dp) :: dbsq(nznt)
      real(kind=dp) :: bsqvac(nznt)
      real(kind=dp) :: rbsq(nznt)
      real(kind=dp) :: curpol
      real(kind=dp) :: curtor
      real(kind=dp) :: rbtor
      real(kind=dp) :: ctor
      real(kind=dp) :: phiedge
      real(kind=dp) :: delbsq
end module magfield

module xstuff
      use stel_kinds, only: dp
      implicit none
      real(kind=dp) :: xcdot(neq)
      real(kind=dp) :: xstore(neq)
      real(kind=dp) :: xc(neq)
      real(kind=dp) :: gc(neq)
end module xstuff
