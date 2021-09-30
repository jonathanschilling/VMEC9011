        subroutine residue(gcr,gcz,gcl,bz0,work)
      include 'name1'
      include 'name0'
      include 'name2'
        common/precond/ard(nsd1,2),arm(nsd1,2),brd(nsd1,2),brm(nsd1,2),
     >        cr(nsd1),azd(nsd1,2),azm(nsd1,2),bzd(nsd1,2),bzm(nsd1,2)
        real gcr(ns,0:nmax,0:mpol1,2),gcz(ns,0:nmax,0:mpol1,2),
     >  gcl(ns,0:nmax,0:mpol1,2),work(mns,6)
*************
*                 IMPOSE M=1 MODE CONSTRAINT (ON Z CS, NTYPE = 1)
*************
        if(fsqz.lt.1.e-8.and.(iter1.ne.1))then
        do 10 n = 1,nmax
        do 10 js= 2,ns
 10     gcz(js,n,1,1) = czero
        endif
*************
*                 CONSTRUCT INVARIANT RESIDUALS
*************
        bz0 = c2p0*hs/bz0**2
        fsql = bz0*sdot(2*mns,gcl,1,gcl,1)
        call getfsq(gcr,gcz,fsqr,fsqz,fnorm,meven)
        fedge = fnorm*(sdot(mnd2,gcr(ns,0,0,1),ns,gcr(ns,0,0,1),ns)
     >        +        sdot(mnd2,gcz(ns,0,0,1),ns,gcz(ns,0,0,1),ns))
*************
*                 PERFORM PRECONDITIONING AND COMPUTE RESIDUES
*************
        call scalfor(gcr,arm,brm,ard,brd,cr,work,work(1,2),
     >  work(1,3),work(1,4),work(1,5))
        call scalfor(gcz,azm,bzm,azd,bzd,cr,work,work(1,2),
     >  work(1,3),work(1,4),work(1,5))
*************
*                 REDUCE R-Z FORCES IF INITIALLY VERY LARGE
*************
        fac = c1p0/(c1p0 + (fsqr+fsqz))
        if(iter1.gt.1)call sscal(4*mns,fac,gc,1)
        call getfsq(gcr,gcz,fsqr1,fsqz1,fnorm1,modd)
        if(iter1.eq.1)call sscal(4*mns,fac,gc,1)
*************
*                 CONSTRUCT PRECONDITIONED (SCALED) LAMBDA FORCES
*************
        do 30 l=1,2*mns
 30     gc(l+4*mns)=faclam(l)*gc(l+4*mns)
        fsql1 = hs*sdot(2*mns,gcl,1,gcl,1)
        return
        end
