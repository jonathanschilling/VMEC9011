
        subroutine jacobian(r1,ru,z1,zu,zu12,ru12,zs,rs,gsqrt,r12,
     >  tau,wint,shalf,ohs,cp25,cp5,czero,nrzt,nznt,irst,meven,modd)
        implicit real*8 (a-h,o-z)
        real r1(nrzt,0:1),ru(nrzt,0:1),z1(nrzt,0:1),zu(nrzt,0:1),rs(*),
     >  zs(*),r12(*),gsqrt(*),ru12(*),zu12(*),shalf(*),tau(*),wint(*)
*************
*                 (RS, ZS)=(R, Z) SUB S, (RU12, ZU12)=(R, Z) SUB THETA(=
*                 AND GSQRT=SQRT(G) ARE DIFFERENCED ON HALF MESH
*************
        do 10 l=2,nrzt
        lm = l-1
        ru12(l) = cp5*(ru(l,meven) + ru(lm,meven) + shalf(l)*
     >                (ru(l,modd ) + ru(lm,modd )))
        zu12(l) = cp5*(zu(l,meven) + zu(lm,meven) + shalf(l)*
     >                (zu(l,modd ) + zu(lm,modd )))
        rs(l)   = ohs*(r1(l,meven) - r1(lm,meven) + shalf(l)*
     >                (r1(l,modd ) - r1(lm,modd )))
        zs(l)   = ohs*(z1(l,meven) - z1(lm,meven) + shalf(l)*
     >                (z1(l,modd ) - z1(lm,modd )))
        r12(l)  = cp5*(r1(l,meven) + r1(lm,meven) + shalf(l)*
     >                (r1(l,modd ) + r1(lm,modd )))
        gsqrt(l)= r12(l) * ( ru12(l)*zs(l) - rs(l)*zu12(l) + cp25*
     >  ( ru(l,modd )*z1(l,modd) + ru(lm,modd )*z1(lm,modd)
     >  - zu(l,modd )*r1(l,modd) - zu(lm,modd )*r1(lm,modd)
     >  +(ru(l,meven)*z1(l,modd) + ru(lm,meven)*z1(lm,modd)
     >  - zu(l,meven)*r1(l,modd) - zu(lm,meven)*r1(lm,modd))/shalf(l)) )
 10     tau(l) = wint(l)*gsqrt(l)
*************
*                 TEST FOR SIGN CHANGE IN JACOBIAN
*************
        taumax = czero
        taumin = czero
        do 20 l=2,nrzt
        taumax = max(tau(l),taumax)
 20     taumin = min(tau(l),taumin)
        if(taumax*taumin.lt.czero)irst=2
        return
        end