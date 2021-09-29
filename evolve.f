
        subroutine evolve(ierflag)
      include 'name1'
      include 'name0'
      include 'name2'

        data ndamp1/10/
*************
*                 COMPUTE MHD FORCES
*************
        call funct3d
        if(iter2.ne.1.or.irst.ne.2)goto 10
        ierflag = 1
        return
*************
*                 COMPUTE DAMPING PARAMETER (DTAU) AND EVOLVE
*                 R, Z, AND LAMBDA ARRAYS IN FOURIER SPACE
*************
 10     if(iter2.ne.iter1)goto 15
        ndamp1 = min(ndamp,15)
        do 35 i = 1,ndamp1
 35     otau(i) = cp15/delt
 15     fsq1 = fsqr1 + fsqz1 + fsql1
        if(iter2.gt.iter1)dtau = min(abs(log(fsq/fsq1)),cp15)
        fsq = fsq1
        if(iter2.le.1)return
        do 25 i = 1,ndamp1-1
 25     otau(i) = otau(i+1)
        if(iter2.gt.iter1)otau(ndamp1) = dtau/delt
        otav = ssum(ndamp1,otau,1)/ndamp1
        dtau = delt*otav
        b1 = c1p0 - cp5*dtau
        fac = c1p0/(c1p0 + cp5*dtau)
        do 20 l = 1,neqs
        xcdot(l) = fac*(xcdot(l)*b1 + delt*gc(l))
 20     xc(l) = xc(l) + xcdot(l)*delt
        return
        end
