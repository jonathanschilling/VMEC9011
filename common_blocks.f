! In VMEC, use is made of a few common block data arragements
! to re-use already allocated memory as temporary workspace.
! This file lists the common blocks and where they were declared
! in the VMEC9011 version.

! bcovar
common/precond/ard(nsd1,2)
               arm(nsd1,2)
               brd(nsd1,2)
               brm(nsd1,2)
               cr(nsd1)
               azd(nsd1,2)
               azm(nsd1,2)
               bzd(nsd1,2)
               bzm(nsd1,2)

! block_data
common /bounds/mupper(3)
               mlower(3)

! name0
common/constant/czero
                cp15
                cp25
                cp5
                cp9
                cp96
                c1p0
                c1p1
                c1p4
                c1p5
                c2p0
                c3p0
                c8p0
                c100p
                c2pm8
                c1pm13

! name2
common/current/buco(nsd)
               bvco(nsd)
               ju(nsd)
               jv(nsd)
common/extfld/potvac(nznt)
              xmpot(nznt)
              xnpot(nznt)
              brv(nznt)
              bphiv(nznt)
              bzv(nznt)
              bscale
              mpmax
              ivac
              ivac2
common/fsqu/fnorm
            fsqr
            fsqz
            fsql
            fnorm1
            fsqr1
            fsqz1
            fsql1
            fsq
            fedge
            wb
            wp
            fsqt(100)
            wdot(100)
            equif(nsd)
common/inputdat/am(0:5)
                ai(0:5)
                raxis(0:nmax)
                zaxis(0:nmax)
                rb(0:nmax,0:mpol1,2)
                zb(0:nmax,0:mpol1,2)
                ftol
                gam
                ac(5)
                ncurr
                nfp
                nstep
                niter
                nvacskip
common/mnarray/xrz3(0:nmax,0:mpol1)
               xrz4(0:nmax,0:mpol1)
               xmpq(0:mpol1,3)
               mscale(0:mpol1)
               nscale(0:nmax1)
               jmin1(0:mpol)
               jmin2(0:mpol)
               jlam(0:mpol)
common/profs/iotas(nsd)
             mass(nsd)
             phips(nsd1)
             pres(nsd)
             vp(nsd)
             phip(nrztd)
common/scalars/dnorm
               hs
               ohs
               twopi
               voli
               ijacob
               itfsq
               iequi
               irst
               iter1
               iter2
               isigng
               meven
               modd
               ndamp
               ns
               ns4
               neqs
               nrzt
               mns
common/scalefac/faclam(mnd2*nsd)
                shalf(nrztd)
                sqrts(nrztd)
                scalxc(neq)
                wint(nrztd)
common/spectra/faccon(0:nmax,0:mpol1)
               specw(nsd)
               tcon(nsd)
common/time/delt
            otav
            otau(15)
            timer(0:10)
common/trignew/cosmu(ntheta2,0:mpol1)
               sinmu(ntheta2,0:mpol1)
               cosmum(ntheta2,0:mpol1)
               sinmum(ntheta2,0:mpol1)
               cosmui(ntheta2,0:mpol1)
               cosmumi(ntheta2,0:mpol1)
               cosnv(nzeta,0:nmax)
               sinnv(nzeta,0:nmax)
               cosnvn(nzeta,0:nmax)
               sinnvn(nzeta,0:nmax)
common/magfield/bsqsav(nznt,3)
                dbsq(nznt)
                bsqvac(nznt)
                rbsq(nznt)
                curpol
                curtor
                rbtor
                ctor
                phiedge
                delbsq
common/xstuff/xcdot(neq)
              xstore(neq)
              xc(neq)
              gc(neq)

! name3
common/realsp/r1(2*nrztd)
              ru(2*nrztd)
              rv(2*nrztd)
              z1(2*nrztd)
              zu(2*nrztd)
              zv(2*nrztd)   ! 12 up to here
              rcon(2*nrztd)
              zcon(2*nrztd) ! 16 up to here
              gcon(nrztd)
              ru0(nrztd)
              zu0(nrztd)
              rcon0(nrztd)
              zcon0(nrztd)

! name4
common/rforces/armn(2*nrztd)
               brmn(2*nrztd)
               crmn(2*nrztd)
               azmn(2*nrztd)
               bzmn(2*nrztd)
               czmn(2*nrztd)
               blmn(2*nrztd)
               clmn(2*nrztd)
