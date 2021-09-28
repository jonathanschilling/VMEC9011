subroutine profil3d(rmn, zmn, intflag)

      use stel_kinds, only: dp

      implicit none

      real(kind=dp), intent(out) :: rmn(ns,0:nmax,0:mpol1,2)
      real(kind=dp), intent(out) :: zmn(ns,0:nmax,0:mpol1,2)
      integer,       intent(in)  :: intflag

      pmass(x) = am(0)+x*(am(1)+x*(am(2)+x*(am(3)+x*(am(4)+x*am(5)))))
      piota(x) = ai(0)+x*(ai(1)+x*(ai(2)+x*(ai(3)+x*(ai(4)+x*ai(5)))))
      pcurr(x) =       x*(ac(1)+x*(ac(2)+x*(ac(3)+x*(ac(4)+x*ac(5)))))

      !          INDEX OF LOCAL VARIABLES
      !
      ! ai       vector of coefficients in phi-series for iota
      ! am       vector of coefficients in phi-series for mass
      ! ac       vector of coefficients in phi-series for toroidal current
      ! iotas    rotational transform , on half radial mesh
      ! jv    (-)toroidal current inside flux surface (vanishes like s)
      ! mass     mass profile on half-grid
      ! phip     radial derivative of phi/(2*pi) on half-grid
      ! phiedge  value of real toroidal flux at plasma edge (s=1)
      ! phips    same as phip , one-dimensional array
      ! pressure pressure profile on half-grid, mass/phip**gam
      ! shalf    sqrt(s), two-dimensional array on half-grid
      ! sqrts    sqrt(s), two-dimensional array on full-grid
      ! wint     two-dimensional array for normalizing angle integration

      ! VERIFY SIGN OF CURPOL CONSISTENT WITH TOROIDAL FLUX
      if (phiedge*curpol .lt. czero) then
        print   *,'CHANGING SIGN OF PHIEDGE '
        write(3,*)'CHANGING SIGN OF PHIEDGE '
        phiedge = -phiedge
      endif

      ! COMPUTE PHIP, MASS (PRESSURE), IOTA PROFILES ON HALF-G
      torcur = real(isigng)*curtor/(dnorm*twopi)
      do js = 2, ns
        phij      = hs*(js - c1p5)
        phips(js) = real(isigng) * phiedge / twopi
        iotas(js) = piota(phij)
        jv(js)    = -real(isigng) * torcur * pcurr(phij)
        mass(js)  = pmass(phij)*(abs(phips(js))*rb(0,0,1))**gam
      enddo

      do js = 1, ns
        do lk = 1, nznt
          loff = js + ns*(lk-1)

          shalf(loff) = sqrt(hs*abs(js - c1p5))
          sqrts(loff) = sqrt(hs*(js-1))
          phip(loff) = phips(js)
        enddo
      enddo

      lk = 0
      do lt = 1, ntheta2
        do lz = 1, nzeta
          do js = 2, ns
            wint(js+ns*lk) = cosmui(lt, 0)/mscale(0)
          enddo
          lk = lk+1
        enddo
      enddo

      do lk = 1, nznt
        shalf(1+ns+ns*(lk-1)) = c1p0
        wint (1   +ns*(lk-1)) = czero
      enddo

      ! SAVE COARSE-MESH, SCALED FOURIER COEFFICIENTS IN GC FOR INTERPOLATION
      ! and INITIALIZE XCDOT=0
      do l = 1, neqs
        gc(l)    = xc(l)
        xc(l)    = czero
        xcdot(l) = czero
      enddo
      do l = 1, 6*mns
        gc(l) = gc(l) * scalxc(l)
      enddo

      ! COMPUTE INITIAL R AND Z FOURIER COEFFICIENTS,
      ! FROM SCALED BOUNDARY VALUES, AND SCALXC ARRAY
      ! (1/SQRTS FACTOR FOR ODD M'S)
      do ntype = 1, 2
        do m = 0, mpol1
          do n = 0, nmax
            t1 = c1p0/(mscale(m)*nscale(n))

            do js = 1, ns
              l = js + ns*(n + nmax1*m) + (ntype-1)*mns

              if (mod(m,2).eq.1) then
                scalxc(l) = c1p0/max(sqrts(js), sqrts(2))
              else
                scalxc(l) = c1p0
              endif

              if (m.eq.0 .and. ntype.eq.1) then
                sm0 = c1p0 - (hs*(js-1))
                rmn(js,n,m,1) = (rb(n,m,1) + (raxis(n)-rb(n,m,1))*sm0)*t1
                zmn(js,n,m,1) = (zb(n,m,1) + (zaxis(n)-zb(n,m,1))*sm0)*t1
              else if (m.eq.0 .and. ntype.eq.2) then
                rmn(js,n,m,2) = czero
                zmn(js,n,m,2) = czero
              else if (m.ne.0) then
                facj = t1*sqrts(js)**m
                rmn(js,n,m,ntype) = rb(n,m,ntype)*facj
                zmn(js,n,m,ntype) = zb(n,m,ntype)*facj
              endif

            enddo
          enddo
        enddo
      enddo

      do l = 1, 4*mns
        scalxc(l+2*mns) = scalxc(l)
      enddo

      ! INTERPOLATE FROM COARSE TO FINE RADIAL GRID
      if (intflag.ne.0) then
        call interp(xc, gc, scalxc, intflag)
      endif

      return

end
