
list (APPEND vmec_sources
	${CMAKE_CURRENT_SOURCE_DIR}/alias.f90
	${CMAKE_CURRENT_SOURCE_DIR}/bcovar.f90
	${CMAKE_CURRENT_SOURCE_DIR}/bss.f90
	${CMAKE_CURRENT_SOURCE_DIR}/convert.f90
	${CMAKE_CURRENT_SOURCE_DIR}/eqfor.f90
	${CMAKE_CURRENT_SOURCE_DIR}/eqsolve.f90
	${CMAKE_CURRENT_SOURCE_DIR}/evolve.f90
	${CMAKE_CURRENT_SOURCE_DIR}/extrap.f90
	${CMAKE_CURRENT_SOURCE_DIR}/fixaray.f90
	${CMAKE_CURRENT_SOURCE_DIR}/forces.f90
	${CMAKE_CURRENT_SOURCE_DIR}/funct3d.f90
	${CMAKE_CURRENT_SOURCE_DIR}/getfsq.f90
	${CMAKE_CURRENT_SOURCE_DIR}/getiota.f90
	${CMAKE_CURRENT_SOURCE_DIR}/interp.f90
	${CMAKE_CURRENT_SOURCE_DIR}/jacobian.f90
	${CMAKE_CURRENT_SOURCE_DIR}/lamcal.f90
	${CMAKE_CURRENT_SOURCE_DIR}/name0.f90
	${CMAKE_CURRENT_SOURCE_DIR}/name1.f90
	${CMAKE_CURRENT_SOURCE_DIR}/name2.f90
	${CMAKE_CURRENT_SOURCE_DIR}/output.f90
	${CMAKE_CURRENT_SOURCE_DIR}/precondn.f90
	${CMAKE_CURRENT_SOURCE_DIR}/pressure.f90
	${CMAKE_CURRENT_SOURCE_DIR}/printout.f90
	${CMAKE_CURRENT_SOURCE_DIR}/profil3d.f90
	${CMAKE_CURRENT_SOURCE_DIR}/profile_functions.f90
	${CMAKE_CURRENT_SOURCE_DIR}/readin.f90
	${CMAKE_CURRENT_SOURCE_DIR}/realsp.f90
	${CMAKE_CURRENT_SOURCE_DIR}/residue.f90
	${CMAKE_CURRENT_SOURCE_DIR}/restart.f90
	${CMAKE_CURRENT_SOURCE_DIR}/rforces.f90
	${CMAKE_CURRENT_SOURCE_DIR}/scalfor.f90
	${CMAKE_CURRENT_SOURCE_DIR}/spectrum.f90
	${CMAKE_CURRENT_SOURCE_DIR}/stel_kinds.f90
	${CMAKE_CURRENT_SOURCE_DIR}/tomnsp.f90
	${CMAKE_CURRENT_SOURCE_DIR}/totzsp.f90
	${CMAKE_CURRENT_SOURCE_DIR}/trid.f90
	${CMAKE_CURRENT_SOURCE_DIR}/vacuum.f90
	${CMAKE_CURRENT_SOURCE_DIR}/vsetup.f90
	${CMAKE_CURRENT_SOURCE_DIR}/wrout.f90
)
set (vmec_sources "${vmec_sources}" PARENT_SCOPE)
