
        subroutine bss(r12,rs,zs,ru12,zu12,shalf,bsubs,lu,lv,
     >  br,bphi,bz,cp25,cp5,nrzt)
      include 'name1'
      include 'name3'
        real r12(*),rs(*),zs(*),ru12(*),zu12(*),shalf(*),
     >  bsubs(*),br(*),bphi(*),bz(*),lu(*),lv(*)
        do 10  l=2,nrzt
        lodd = l + nrzt
        rv12 =      cp5*(   rv(l)+   rv(l-1)
     >       +shalf(l)*(rv(lodd)+rv(lodd-1)))
        zv12 =      cp5*(   zv(l)+   zv(l-1)
     >       +shalf(l)*(zv(lodd)+zv(lodd-1)))
        gsu= rs(l)*ru12(l) + zs(l)*zu12(l)
     >        +cp25*(r1(lodd)*ru(lodd)+r1(lodd-1)*ru(lodd-1)
     >             +z1(lodd)*zu(lodd)+z1(lodd-1)*zu(lodd-1)
     >        +(    r1(lodd)*ru(l)   +r1(lodd-1)*ru(l-1)
     >        +     z1(lodd)*zu(l)   +z1(lodd-1)*zu(l-1))/shalf(l))
        gsv= rs(l)*rv12    + zs(l)*zv12
     >        +cp25*(r1(lodd)*rv(lodd)+r1(lodd-1)*rv(lodd-1)
     >             +z1(lodd)*zv(lodd)+z1(lodd-1)*zv(lodd-1)
     >        +(    r1(lodd)*rv(l)   +r1(lodd-1)*rv(l-1)
     >        +     z1(lodd)*zv(l)   +z1(lodd-1)*zv(l-1))/shalf(l))
        br(l)   = lv(l)*ru12(l) + lu(l)*rv12
        bphi(l) =                 lu(l)*r12(l)
        bz(l)   = lv(l)*zu12(l) + lu(l)*zv12
 10     bsubs(l)= lv(l)*gsu     + lu(l)*gsv
        return
        end
