
all: vmec

vmec: crlib64.f \
      vac.f \
      name0 name1 name2 name3 name4 \
      alias.f       \
      bss.f         \
      eqfor.f       \
      extrap.f      \
      funct3d.f     \
      interp.f      \
      output.f      \
      printout.f    \
      residue.f     \
      spectrum.f    \
      trid.f        \
      vsetup.f      \
      bcovar.f      \
      convert.f     \
      eqsolve.f     \
      fixaray.f     \
      getfsq.f      \
      jacobian.f    \
      precondn.f    \
      profil3d.f    \
      restart.f     \
      tomnsp.f      \
      wrout.f       \
      block_data.f  \
      evolve.f      \
      forces.f      \
      getiota.f     \
      lamcal.f      \
      pressure.f    \
      readin.f      \
      scalfor.f     \
      totzsp.f      \
      vmec.f
	gfortran vac.f -c -o vac.f.o
	gfortran -std=legacy crlib64.f -c -o crlib64.f.o
	gfortran -fdefault-real-8 -std=legacy alias.f      -c -o alias.f.o
	gfortran -fdefault-real-8 -std=legacy bss.f        -c -o bss.f.o
	gfortran -fdefault-real-8 -std=legacy eqfor.f      -c -o eqfor.f.o
	gfortran -fdefault-real-8 -std=legacy extrap.f     -c -o extrap.f.o
	gfortran -fdefault-real-8 -std=legacy funct3d.f    -c -o funct3d.f.o
	gfortran -fdefault-real-8 -std=legacy interp.f     -c -o interp.f.o
	gfortran -fdefault-real-8 -std=legacy output.f     -c -o output.f.o
	gfortran -fdefault-real-8 -std=legacy printout.f   -c -o printout.f.o
	gfortran -fdefault-real-8 -std=legacy residue.f    -c -o residue.f.o
	gfortran -fdefault-real-8 -std=legacy spectrum.f   -c -o spectrum.f.o
	gfortran -fdefault-real-8 -std=legacy trid.f       -c -o trid.f.o
	gfortran -fdefault-real-8 -std=legacy vsetup.f     -c -o vsetup.f.o
	gfortran -fdefault-real-8 -std=legacy bcovar.f     -c -o bcovar.f.o
	gfortran -fdefault-real-8 -std=legacy convert.f    -c -o convert.f.o
	gfortran -fdefault-real-8 -std=legacy eqsolve.f    -c -o eqsolve.f.o
	gfortran -fdefault-real-8 -std=legacy fixaray.f    -c -o fixaray.f.o
	gfortran -fdefault-real-8 -std=legacy getfsq.f     -c -o getfsq.f.o
	gfortran -fdefault-real-8 -std=legacy jacobian.f   -c -o jacobian.f.o
	gfortran -fdefault-real-8 -std=legacy precondn.f   -c -o precondn.f.o
	gfortran -fdefault-real-8 -std=legacy profil3d.f   -c -o profil3d.f.o
	gfortran -fdefault-real-8 -std=legacy restart.f    -c -o restart.f.o
	gfortran -fdefault-real-8 -std=legacy tomnsp.f     -c -o tomnsp.f.o
	gfortran -fdefault-real-8 -std=legacy wrout.f      -c -o wrout.f.o
	gfortran -fdefault-real-8 -std=legacy block_data.f -c -o block_data.f.o
	gfortran -fdefault-real-8 -std=legacy evolve.f     -c -o evolve.f.o
	gfortran -fdefault-real-8 -std=legacy forces.f     -c -o forces.f.o
	gfortran -fdefault-real-8 -std=legacy getiota.f    -c -o getiota.f.o
	gfortran -fdefault-real-8 -std=legacy lamcal.f     -c -o lamcal.f.o
	gfortran -fdefault-real-8 -std=legacy pressure.f   -c -o pressure.f.o
	gfortran -fdefault-real-8 -std=legacy readin.f     -c -o readin.f.o
	gfortran -fdefault-real-8 -std=legacy scalfor.f    -c -o scalfor.f.o
	gfortran -fdefault-real-8 -std=legacy totzsp.f     -c -o totzsp.f.o
	gfortran -fdefault-real-8 -std=legacy vmec.f       -c -o vmec.f.o
	gfortran crlib64.f.o vac.f.o \
	         alias.f.o      \
	         bss.f.o        \
	         eqfor.f.o      \
	         extrap.f.o     \
	         funct3d.f.o    \
	         interp.f.o     \
	         output.f.o     \
	         printout.f.o   \
	         residue.f.o    \
	         spectrum.f.o   \
	         trid.f.o       \
	         vsetup.f.o     \
	         bcovar.f.o     \
	         convert.f.o    \
	         eqsolve.f.o    \
	         fixaray.f.o    \
	         getfsq.f.o     \
	         jacobian.f.o   \
	         precondn.f.o   \
	         profil3d.f.o   \
	         restart.f.o    \
	         tomnsp.f.o     \
	         wrout.f.o      \
	         block_data.f.o \
	         evolve.f.o     \
	         forces.f.o     \
	         getiota.f.o    \
	         lamcal.f.o     \
	         pressure.f.o   \
	         readin.f.o     \
	         scalfor.f.o    \
	         totzsp.f.o     \
	         vmec.f.o       \
	         -o vmec


