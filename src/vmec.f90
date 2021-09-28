program vmec

      !             D * I * S * C * L * A * I * M * E * R
      !
      ! You are using a BETA version of the program VMEC, which is currently
      ! under development by S. P. Hirshman at the Fusion Energy Division of
      ! Oak Ridge National Laboratory.  Please report any problems or comments
      ! to him.  As a BETA version, this program is subject to change
      ! and improvement without notice.
      !
      ! THIS PROGRAM - VMEC (Variational Moments Equilibrium Code) -
      ! SOLVES THREE-DIMENSIONAL MHD EQUILIBRIUM EQUATIONS USING FOURIER SPECTRAL (MOMENTS) METHODS.
      ! A CYLINDRICAL COORDINATE REPRESENTATION IS USED (R-Z COORDINATES).
      ! THE POLOIDAL ANGLE VARIABLE IS RENORMALIZED THROUGH THE STREAM FUNCTION LAMBDA,
      ! WHICH IS SELF-CONSISTENTLY DETERMINED AND DIFFERENCED VARIATIONALLY ON THE FULL-RADIAL MESH.
      ! THE POLOIDAL ANGLE IS DETERMINED BY MINIMIZING <M> = m**2 S(m) , WHERE S(m) = Rm**2 + Zm**2 .
      ! AN EVEN-ODD DECOMPOSITION IN THE POLOIDAL MODE NO. OF R,Z, AND LAMDA
      ! IS USED TO IMPROVE RADIAL RESOLUTION.
      ! A FREE-BOUNDARY OPTION IS AVAILABLE (FOR nvac > 0),
      ! WITH A USER-SUPPLIED SUBROUTINE "VACUUM" NEEDED TO COMPUTE THE PLASMA BOUNDARY VALUE OF B**2.
      !
      ! Added features since last edition
      ! 1.  Implemented preconditioning algorithm for R,Z
      ! 2.  The physical (unpreconditioned) residuals are used
      !     to determine the level of convergence.
      ! 3.  The original (MOMCON) scaling of lambda is used, i.e.,
      !     Bsupu = phip*(iota - lamda[sub]v)/sqrt(g). This is needed to
      !     maintain consistency with the time-stepper for arbitrary PHI.
      !
      ! WRITTEN BY S. P. HIRSHMAN (8/28/85 - REVISED 3/1/86) BASED ON
      ! 1. S. P. Hirshman and J. C. Whitson, Phys. Fluids 26, 3553 (1983
      ! 2. S. P. Hirshman and H. K. Meier, Phys. Fluids 28, 1387 (1985).
      ! 3. S. P. Hirshman and D. K. Lee, Comp. Phys. Comm. 39, 161 (1986

      ! input/output files:
      ! unit | filename
      ! ---------------
      !   7  | indata
      !   3  | threed1
      !   8  | wout

      use name1, only: nsd
      implicit none

      character*80 :: werror(0:6)

      integer      :: ierflag
      integer      :: intflag
      integer      :: nsin

      !         PARAMETER STATEMENT DIMENSIONS AND RECOMMENDED VALUES
      !
      ! mpol    number of poloidal modes, starting at m=0 (note: mpol
      ! nmax    mode number of largest toroidal mode (normed to nfp),
      !         with -nmax <= n <= nmax
      ! ntheta  number of theta mesh points (at least 2*mpol+4)
      !         (note: in fft routines in VACUUM, ntheta and nzeta
      !         must be even with no prime factors greater than 5.)
      ! nzeta   number of zeta mesh points [at least 2*nmax+4 if nmax>
      !         or =1 if nmax = 0 (axisymmetric case)]
      ! nvac    turns on fixed (=0) or free(=1) boundary option

      !         INDEX OF LOCAL VARIABLES
      !
      ! ierflag specifies error condition if nonzero
      ! intflag = 0,    compute xc array from boundary data
      !         = nsin, compute xc array by interpolation

      OPEN(UNIT=7, FILE='INDATA',  STATUS='OLD',     ERR=901)
      OPEN(UNIT=3, FILE='THREED1', STATUS='NEW',     ERR=902)
      OPEN(UNIT=8, FILE='WOUT',    STATUS='UNKNOWN', ERR=903)
      goto 904
 901  print *,' INDATA FILE IS MISSING'
      stop
 902  print *,' THREED1 FILE ALREADY EXISTS: RENAME OR DELETE IT'
      stop
 903  print *,' WOUT FILE ALREADY EXIST: RENAME OR DELETE IT'
      stop
 904  continue

      print *,   'THIS IS THE PRECONDITIONED VMEC.FULL CODE: VMEC9011'
      WRITE(3,*) 'THIS IS THE PRECONDITIONED VMEC.FULL CODE: VMEC9011'

      ! INITIALIZE ALL VARIABLES SO VMEC CAN BE CALLED AS A SUBROUTINE.
      ierflag = 0
      intflag = 0
      call vsetup

      ! READ INPUT FILE INDATA AND STORE IN INPUTDAT COMMON BLOC
      call readin(nsin,ierflag)
      if (ierflag.ne.0) then
        goto 10
      endif

      ! COMPUTE INVARIANT ARRAYS
      ! IF VMEC IS USED AS A SUBROUTINE, CALL FIX_ARAY ONLY ONCE AT SETUP.
      call fixaray(ierflag)
      if (ierflag.ne.0) then
        goto 10
      endif

      ! MAKE INITIAL CALL TO EQUILIBRIUM SOLVER
      if (nsin.lt.2 .or. nsin.gt.nsd) then
        nsin = nsd
      endif
      call eqsolve(nsin, intflag, ierflag)

      ! PERFORM COARSE-TO-FINE MESH INTERPOLATION
      if (nsin.lt.nsd .and. ierflag.eq.0) then
        intflag = nsin
        call eqsolve(nsd, intflag, ierflag)
      end if

      ! WRITE OUTPUT TO FILES THREED1, WOUT
      if (ierflag.eq.0 .or. ierflag.eq.1 .or. ierflag.eq.5) then
        call output
      endif

 10   werror(0)=' EXECUTION TERMINATED NORMALLY'
      werror(1)=' INITIAL JACOBIAN CHANGED SIGN (NEED A BETTER GUESS)'
      werror(2)=' MORE THAN 100 JACOBIAN ITERATIONS (DECREASE DELT)'
      werror(3)=' m IN BOUNDARY ARRAY EXCEEDS mpol-1 (INCREASE MPOL)'
      werror(4)=' n IN BOUNDARY ARRAYS OUTSIDE nmax,(-nmax) RANGE'
      werror(5)=' ASSUMED SIGN OF JACOBIAN (ISIGNG) IS WRONG'
      print *,   werror(ierflag)
      write(3,*) werror(ierflag)

      call exit

end program vmec
