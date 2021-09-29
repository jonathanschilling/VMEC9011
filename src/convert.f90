subroutine convert(rmnc,zmns,lmns,xm,xn,js, &
                   rmncc,rmnss, zmncs,zmnsc, lmncs,lmnsc)

      include 'name1'
      include 'name0'
      include 'name2'

      real(kind=dp), intent(out) :: rmnc(*)
      real(kind=dp), intent(out) :: zmns(*)
      real(kind=dp), intent(out) :: lmns(*)
      real(kind=dp), intent(out) :: xm(*)
      real(kind=dp), intent(out) :: xn(*)
      real(kind=dp), intent(in)  :: rmncc(ns,0:nmax,0:mpol1)
      real(kind=dp), intent(in)  :: rmnss(ns,0:nmax,0:mpol1)
      real(kind=dp), intent(in)  :: zmncs(ns,0:nmax,0:mpol1)
      real(kind=dp), intent(in)  :: zmnsc(ns,0:nmax,0:mpol1)
      real(kind=dp), intent(in)  :: lmncs(ns,0:nmax,0:mpol1)
      real(kind=dp), intent(in)  :: lmnsc(ns,0:nmax,0:mpol1)

      ! CONVERTS INTERNAL MODE REPRESENTATION TO STANDARD
      ! FORM FOR OUTPUT (COEFFICIENTS OF cos(mu-nv), sin(mu-nv

      mn = 0
      do m = 0,mpol1

        nmin0 = -nmax
        if (m.eq.0) &
          nmin0 = 0

        do n = nmin0,nmax
          n1 = iabs(n)

          t1 = mscale(m)*nscale(n1)

          mn = mn + 1

          xm(mn) = real(m)
          xn(mn) = real(n*nfp)

          if( m.eq.0 )then
            rmnc(mn) = t1*rmncc(js,n,m)
            zmns(mn) =-t1*zmncs(js,n,m)
            lmns(mn) =-t1*lmncs(js,n,m)
          else if( n.eq.0 )then
            rmnc(mn) = t1*rmncc(js,n,m)
            zmns(mn) = t1*zmnsc(js,n,m)
            lmns(mn) = t1*lmnsc(js,n,m)
          else if( js.gt.1 )then
            sgn = real(n/n1)
            rmnc(mn) = cp5*t1*(rmncc(js,n1,m) + sgn*rmnss(js,n1,m))
            zmns(mn) = cp5*t1*(zmnsc(js,n1,m) - sgn*zmncs(js,n1,m))
            lmns(mn) = cp5*t1*(lmnsc(js,n1,m) - sgn*lmncs(js,n1,m))
          else if( js.eq.1 )then
            rmnc(mn) = czero
            zmns(mn) = czero
            lmns(mn) = czero
          endif

        enddo
      enddo

      return

end
