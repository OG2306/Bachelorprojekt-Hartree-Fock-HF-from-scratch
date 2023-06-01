/********************************************************
I denne .h-fil er defineret klassen HF, der udregner
HF-energien af et molekyle gennem SCF-proceduren.
Et objekt af klassen konstrueres ud fra et molekyle og en filsti 
til den mappe, der indeholder integralerne nødvendige for dannelsen af Fock-matricen m.m.
------------BESKRIVELSE AF VIGTIGE DATAMEDLEMMER------------
smx: Matricen med alle overlapintegralerne
tmx: Matricen med alle én-elektron-integralerne for kinetisk energi
vmx: Matricen med alle én-elektron-integralerne for kernetiltrækningen
hcore: Matricen tmx + vmx
jmx: Matricen med Coulomb-integraler
kmx: Matricen med exchange-integraler
fmx: Fock-matricen: fmx = hcore + 2jmx - kmx
orthsmx: Orthogonaliseringsmatricen, der bruges i diagonaliseringen af Fock-matricen
cmx: Matricen med alle koefficienterne til LCAO-MO'erne
dmx: Densitetsmatricen, der udregnes fra cmx
eritensor: En tensor med alle to-elektron-integralerne. eritensor[0][1][2][3] = (1 2|3 4) (Mulliken-notation)
moenergy: En vektor med de udregnede molekylorbitalenergier (egenværdierne til Fock-matricen). Bestemmes i diagonalizefmx
nobasisfcn: Antallet af basisfunktioner, der er brugt i beregningen. Bruges mange steder i medlemsfunktionerne
------------BESKRIVELSE AF VIGTIGE MEDLEMSFUNKTIONER------------
HFSCF: Det er denne funktion, der laver normale HF SCF
HFSCFDIIS: Det er denne funktion, der laver DIIS-accelereret HF SCF
computeorthsmx: Funktionen danner orthogonaliseringsmatricen
diagonalizefmx: Funktionen diagonaliserer Fock-matricen og bestemmer koefficienterne til LCAO-MO'erne. Bestemmer samtidigt molekylorbitalenergierne
computedmx: Funktionen laver densitetsmatricen ud fra et sæt af basisfunktionkoefficienter
computeenergy: Funktionen udregner den elektroniske energi ud fra en Fock- og densitetsmatrix
computefmx: Funktionen laver en ny Fock-matrix.
MP2 og MP2Quick: Funktionerne udregner MP2-korrektionen til HF-energien med en langsom og hurtigere algoritme.
********************************************************/

#ifndef _HF
#define _HF

#include <fstream>
#include "Molecule.h"
#include "../sources/Eigen/Dense"
#include "../sources/Eigen/Eigenvalues"
#include "../sources/Eigen/Core"
using namespace std;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> DynMatrix;

class HF {
private:
    ////////////////////////////////////////////////////////////////////////////////
    //PRIVATE DATAMEDLEMMER
    ////////////////////////////////////////////////////////////////////////////////

    Molecule& mymolecule; //Molekylet hvis energi skal beregnes
    const char* const integraldir; //Filstien til mappen med integraldata
    char basissetname[20]; //Navnet på basissættet, der bruges til energiberegningen
    int nobasisfcn; //Antallet af basisfunktioner, der bruges i energiberegningen
    DynMatrix smx, tmx, vmx, hcore, jmx, kmx, fmx, orthsmx, cmx, dmx; //Relevante matricer
    double**** eritensor; //Tensoren med alle to-elektron-integralerne. Svarer til en 2D matrix, hvor hvert indeks er en 2D-matrix selv
    double* moenergy; //En vektor med molekylorbitalenergierne
    double elecenergy; //Den elektroniske energi
    bool energyconverged; //Er energien konvergeret?
    DynMatrix muxmx, muymx, muzmx; //Matricer til dipolintegralerne

    ////////////////////////////////////////////////////////////////////////////////
    //PRIVATE LOAD-FUNKTIONER
    ////////////////////////////////////////////////////////////////////////////////

    /*Load navnet på det anvendte basissæt og antallet af brugte basisfunktioner*/
    void loadbasissetinfo();

    /*Load én-elektron-integraler fra en datafil. Stien til mappen, hvor datafilen ligger, er gemt i integraldir, 
    og kun filnavnet skal bruges som parameter. Anden parameter er en pointer til den matrix, integralerne skal gemmes i.*/
    void loadintegrals(const char* const fname, DynMatrix* mx);

    /*Load Coulomb- og exchange-integralerne*/
    void loaderiintegrals();

    ////////////////////////////////////////////////////////////////////////////////
    //PRIVATE SET-FUNKTIONER
    ////////////////////////////////////////////////////////////////////////////////

    /*Lav Hcore ud fra T- og V-matricen. Kaldes i constructoren*/
    void sethcore() {hcore = tmx + vmx;}

    /*Angiv startgættet for Fock-matricen*/
    void setfmxstartguess(DynMatrix* startguess) {fmx = *startguess;}

    ////////////////////////////////////////////////////////////////////////////////
    //PRIVATE COMPUTE-FUNKTIONER
    ////////////////////////////////////////////////////////////////////////////////

    /*Lav S^(-1/2)-matricen (orthsmx) ud fra overlapmatricen. Kaldes i constructoren*/
    void computeorthsmx();

    /*Diagonalisér Fock-matricen og udregn koefficienterne i LCAO-MO'erne
    Funktionen tager et input af hensyn til DIIS*/
    void diagonalizefmx(DynMatrix* mx);

    /*Lav densitetsmatricen ud fra basisfunktionkoefficienterne*/
    void computedmx();

    /*Udregn den elektroniske energi ud fra Fock- og densitetsmatricen
    Funktionen tager et input (en Fock-matrix) af hensyn til DIIS*/
    void computeenergy(DynMatrix* mx);

    /*Lav Coulomb- og exchange-matricerne ud fra en densitetsmatrix og to-elektron-integralerne.
    Udregn herefter Fock-matricen i henhold til formlen fmx = hcore + 2jmx - kmx*/
    void computefmx();

    /*Udregn root-mean-squared-forskellen for to densitetsmatricer*/
    double computermsd(DynMatrix* mx1, DynMatrix* mx2);

public:
    /*Et objekt af HF-klassen konstrueres ud fra et Molecule-objekt og
    en streng, der angiver stien til den mappe, der indeholder al integraldata*/
    HF(Molecule& molecule, const char* const dir); //Constructor
    ~HF(); //Destructor

    ////////////////////////////////////////////////////////////////////////////////
    //PUBLIC PRINTFUNKTIONER
    ////////////////////////////////////////////////////////////////////////////////

    /*Print navnet på det anvendte basissæt til en fil*/
    void printbasissetname(ofstream& outfile) const {outfile << basissetname;}

    /*Print antallet af basisfunktioner, der bruges i beregningen, til en fil*/
    void printnobasisfcn(ofstream& outfile) const {outfile << nobasisfcn;}

    /*Print en matrix med dimensioner (nobasisfcn, nobasisfcn) til en fil*/
    void printmx(ofstream& outfile, DynMatrix* mx);

    ////////////////////////////////////////////////////////////////////////////////
    //PUBLIC GET-FUNKTIONER
    ////////////////////////////////////////////////////////////////////////////////

    /*Returnerer antallet af basisfunktioner brugt i beregningen*/
    int getnobasisfcn() const {return nobasisfcn;}

    /*Returnerer den elektroniske energy*/
    double getelecenergy() const {return elecenergy;}

    /*Returnerer addresser til diverse matricer*/
    DynMatrix* getsmx() {return &smx;}
    DynMatrix* gettmx() {return &tmx;}
    DynMatrix* getvmx() {return &vmx;}

    ////////////////////////////////////////////////////////////////////////////////
    //PUBLIC SET-FUNKTIONER
    ////////////////////////////////////////////////////////////////////////////////

    /*Angiv om energien er konvergeret eller ej*/
    void setenergyconverged(bool flag) {energyconverged = flag;}

    ////////////////////////////////////////////////////////////////////////////////
    //PUBLIC LOAD-FUNKTIONER
    ////////////////////////////////////////////////////////////////////////////////

    void loaddipoleintegrals(); //Indlæs dipolintegralerne til de rette matricer

    ////////////////////////////////////////////////////////////////////////////////
    //PUBLIC COMPUTE-FUNKTIONER
    ////////////////////////////////////////////////////////////////////////////////

    /*Udregn HF-energien af molekylet gennem SCF-proceduren.
    maxiter er det maksimale antal SCF-iterationer.
    Hvis root mean square-forskellen mellem alle elementerne i to på hinanden følgende 
    densitetsmatricer falder under delta1, antages energien at være konvergeret.
    Hvis energiforskellen mellem to på hinanden følgende elek. energier falder under delta2,
    antages energien også at være konvergeret.
    De elektroniske energier, kernefrastødningen og den endelige HF-energi printes til filen outfile.*/
    void HFSCF(int maxiter, double delta1, double delta2, ofstream& outfile);

    /*DIIS-accelereret HF SCF.*/
    void HFSCFDIIS(int maxiter, double delta1, double delta2, ofstream& outfile, int nosavedmcs);

    ////////////////////////////////////////////////////////////////////////////////
    //VENNER
    ////////////////////////////////////////////////////////////////////////////////

    /*Udregn MP2-korrektionen til HF-energien af et molekyle, hvis HF-energi er bestemt
    via et HF objekt og HFSCF-metoden. MP2 skalerer med antallet af basisfunktioner 
    opløftet til 8. MP2Quick skalerer med antallet af basisfunktioner opløftet til 5.*/
    friend void MP2(HF&, ofstream&);
    friend void MP2Quick(HF&, ofstream&);
    
    /*Udregn (det elektroniske) dipolmoment af et molekyle, hvis (konvergerede) densitetsmatrix er bestemt via et HF objekt
    og HFSCF-metoden. Funktionen kræver, at de nødvendige dipolintegraler er udregnet og
    udskrevet til filerne mux.dat, muy.dat og muz.dat. Filerne skal ligge i samme mappe som
    de integraler, der er brugt i HF-beregningen (filstien til mappen er gemt i integraldir,
    som er et datamedlem af HF-objektet.)*/
    friend void dipolemoment(HF&, ofstream&);

    /*DIIS-objektet skal tilgå private datamedlemmer fra HF-klassen som bl.a. fmx (Fock-matricen), dmx (densitetsmatricen) og smx (overlapmatricen)*/
    friend class DIIS;

    /*CI-objekter skal have adgang til information om basissætstørrelse, integraler m.m.*/
    friend class CI;

    /*TDHF-objekter skal have adgang til bl.a. cmx-matricen med LCAO-MO-koefficienter*/
    friend class TDHF;
};

#endif