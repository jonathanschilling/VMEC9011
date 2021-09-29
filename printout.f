
        subroutine printout(i,w0,r00)
      include 'name1'
      include 'name0'
      include 'name2'
        betav = wp/wb
        w = w0*twopi*twopi/nfp
        avm = czero
        den = czero
        specw(1) = c1p0
        call spectrum(xc,xc(1+mns),xc(1+2*mns),xc(1+3*mns))
        do 10 j=2,ns
        den = den+phips(j)
 10     avm = avm+phips(j)*(specw(j)+specw(j-1))
        avm = cp5*avm/den
        if( ivac.ge.1 )delbsq=sdot(nznt,dbsq,1,wint(2),ns)/
     >                 sdot(nznt,bsqsav(1,3),1,wint(2),ns)
        if(i.eq.1.and.nvac.ne.0)print 20
        if(i.eq.1.and.nvac.ne.0)write(3,15)
        if(i.eq.1.and.nvac.eq.0)print 30
        if(i.eq.1.and.nvac.eq.0)write(3,25)
 15     format(/,' ITER    FSQR      FSQZ      FSQL      fsqr      ',
     >  'fsqz      DELT     R00(0)        WMHD        BETA     ',
     >  '<M>   DEL-BSQ   FEDGE',/)
 20     format(/,' ITER    FSQR      FSQZ      FSQL      ',
     >  'R00(0)     WMHD      DEL-BSQ',/)
 25     format(/,' ITER    FSQR      FSQZ      FSQL      fsqr      ',
     >  'fsqz      DELT     R00(0)        WMHD        BETA     <M>',/)
 30     format(/,' ITER    FSQR      FSQZ      FSQL      fsqr      ',
     >  'fsqz     R00(0)     WMHD',/)
        if(nvac.ne.0)goto 35
        PRINT 45, I,FSQR,FSQZ,FSQL,FSQR1,FSQZ1,R00,W
        write(3,40)i,fsqr,fsqz,fsql,fsqr1,fsqz1,delt,r00,w,betav,avm
        return
 35     CONTINUE
        PRINT 50, I,FSQR,FSQZ,FSQL,R00,W,DELBSQ
        write(3,40)i,fsqr,fsqz,fsql,fsqr1,fsqz1,delt,r00,w,
     >  betav,avm,delbsq,fedge
 40     format(i5,1p6e10.2,1pe11.4,1pe15.8,1pe9.2,0pf7.3,1p2e9.2)
 45     format(i5,1p5e10.2,1pe10.3,1p2e11.4)
 50     format(i5,1p3e10.2,1pe10.3,1p2e11.4)
        return
        end
