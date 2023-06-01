add_library(MolekylBibliotek
	"${PROJECT_SOURCE_DIR}/sources/Atom.cpp"
	"${PROJECT_SOURCE_DIR}/sources/AtomicNumber.cpp"
	"${PROJECT_SOURCE_DIR}/sources/Molecule.cpp"
	"${PROJECT_SOURCE_DIR}/sources/HF.cpp"
	"${PROJECT_SOURCE_DIR}/sources/DIIS.cpp"
	"${PROJECT_SOURCE_DIR}/sources/CI.cpp"
	"${PROJECT_SOURCE_DIR}/sources/TDHF.cpp"
)

target_include_directories(MolekylBibliotek PUBLIC
	"${PROJECT_SOURCE_DIR}/include"
	"${PROJECT_SOURCE_DIR}/sources/Eigen"
)
