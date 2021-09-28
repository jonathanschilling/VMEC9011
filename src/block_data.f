************************************************************************
*               INITIALIZE COMMON BLOCK CONSTANTS
************************************************************************
        block data
      include 'name1'
      include 'name0'
      include 'name2'
        parameter(mup2=nmax+nmax1,mlo3=mup2+1,
     >  zerod=0.0e+00,oned=1.0e+00,cp707d=7.071e-01)
        common /bounds/ mupper(3),mlower(3)

        data mupper/nmax,mup2,mnd1/, mlower/0,nmax1,mlo3/
        DATA NDAMP/10/, NS4/25/, ISIGNG/-1/,
     >  meven/0/, modd/1/, jmin1/1,1,mpol1*2/, jmin2/1,2,mpol1*3/,
     >  jlam/2,3,mpol1*3/,phips(1)/zerod/,mscale/cp707d,mpol1*oned/,
     >  nscale/cp707d,nmax1*oned/,cp15/1.5e-01/,
     >  cp25/2.5e-01/,cp5/5.0e-01/,cp9/9.0e-01/,c2pm8/2.e-08/,
     >  cp96/9.6e-01/,c1p0/oned/,c1p1/1.1e+00/,c1p4/1.4e+00/,
     >  c1p5/1.5e+00/,c2p0/2.0e+00/,c1pm13/1.0e-13/,
     >  c3p0/3.0e+00/,c8p0/8.0e+00/,c100p/1.00e+02/,czero/zerod/
        end
