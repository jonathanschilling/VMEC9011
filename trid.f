
        subroutine trid(a,d,b,c,gam,alf,jmin,nn,ns)
      include 'name1'
*************
*                 SOLVES B(I)*X(I-1)+D(I)*X(I)+A(I)*X(I+1)=C(I), I=IN,NN
*                 AND RETURNS ANSWER IN C(I)
*************
        real a(ns,0:mnd1),b(ns,0:mnd1),c(ns,0:mnd1),d(ns,0:mnd1),
     >  alf(0:mnd1,ns),gam(0:mnd1,ns)
        integer jmin(0:3)
        common /bounds/ mupper(3),mlower(3)
*************
*                 SEPARATE M=0 (IMODES=1), M=1 (IMODES=2), M>1 (IMODES=3
*************
        do 100 imodes = 1,3
        in = jmin(imodes-1)
        in1 = in + 1
        do 10 mn = mlower(imodes),mupper(imodes)
 10     gam(mn,in) = d(in,mn)
        do 20 i=in1,nn
        do 20 mn = mlower(imodes),mupper(imodes)
        alf(mn,i-1) = a(i-1,mn)/gam(mn,i-1)
 20     gam(mn,i)   = d(i,mn) - b(i,mn)*alf(mn,i-1)
        do 30 mn = mlower(imodes),mupper(imodes)
 30     c(in,mn) = c(in,mn)/gam(mn,in)
        do 40 i=in1,nn
        do 40 mn = mlower(imodes),mupper(imodes)
 40     c(i,mn) = (c(i,mn) - b(i,mn)*c(i-1,mn))/gam(mn,i)
        n2 = nn + in
        do 50 i=in1,nn
        i1 = n2 -i
        do 50 mn = mlower(imodes),mupper(imodes)
 50     c(i1,mn) = c(i1,mn) - alf(mn,i1)*c(i1+1,mn)
 100    continue
        return
        end
