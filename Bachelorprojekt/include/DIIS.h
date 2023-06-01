/********************************************************
I denne .h-fil er defineret klassen DIIS, der bruges til
at accelerere energikonvergeringen under HF SCF. Et objekt
af klassen laves ud fra en pointer til et HF-objekt, hvis
Fock-, densitets- og overlapmatricer skal tilgås. Desuden
tager constructoren en int som input, der angiver, hvor
mange fejl- og Fock-matricer, der skal gemmes maksimalt.
------------BESKRIVELSE AF VIGTIGE DATAMEDLEMMER------------
myhf: Det HF-objekt, hvis HF SCF-beregning skal accelereres
maxdimensions: Det maksimale antal fejl- og Fock-matricer, et DIIS-objekt gemmer
curdimension: Det antal fejl- og Fock-matricer, der er gemt i et DIIS-objekt på et givent tidspunkt (curdimensions <= maxdimensions)
emx: Et array af fejlvektorer (matricer).
fmxhistory: Et array af gemte Fock-matricer. Fock-matrix-gættet i DIIS skrives som en lineær kombination af tidligere Fock-matricer
bmx: En matrix, der dannes ud fra de gemte fejlmatricer og som bruges til at finde koefficienterne til den lineære kombination af Fock-matricer
lccoeff: Koefficienterne til den lineære kombination af Fock-matricer
lcfmx: Fock-matrix-gættet, der er lavet som en lineær kombination af tidligere Fock-matricer ("rigtige" Fock-matricer - ikke tidligere lineære kombinationer)
------------BESKRIVELSE AF VIGTIGE MEDLEMSFUNKTIONER------------
savefmx: Gemmer den Fock-matrix, der gives som første argument til funktionen, i den position i arrayet fmxhistory, der gives som andet argument
computeemx: Udregner en fejlmatrix ud fra Fock-, densitets- og overlapmatricerne i myhf-objektet. Fejlmatricen gemmes i den position i arrayet emx, der gives som argument til funktionen
computebmx: Udregner B-matricen ud fra de gemte fejlmatricer
solvebmx: Bestemmer ud fra B-matricen koefficienterne til den lineære kombination af Fock-matricer: Løser et ligningssystem af typen Bc = x
computelcfmx: Bestemmer Fock-matrix-gættet som en lineær kombination af tidligere Fock-matricer ("rigtige" Fock-matricer - ikke tidligere lineære kombinationer)
********************************************************/
#ifndef _DIIS
#define _DIIS

#include "HF.h"
#include "../sources/Eigen/Dense"
#include "../sources/Eigen/Eigenvalues"
#include "../sources/Eigen/Core"
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> DynMatrix;

class DIIS {
private:
    ////////////////////////////////////////////////////////////////////////////////
    //PRIVATE DATAMEDLEMMER
    ////////////////////////////////////////////////////////////////////////////////

    HF* myhf; //HF-objektet, hvis Fock-, densitets-, og overlapmatricer skal tilgås
    int maxdimensions; //Antallet af fejl- og Fock-matricer, der gemmes maksimalt
    int curdimensions; //Antallet af fejl- og Fock-matricer, der er gemt lige nu
    DynMatrix* emx; //Fejlmatricerne
    DynMatrix* fmxhistory; //De gemte Fock-matricer
    DynMatrix bmx; //B-matricen, der bruges i DIIS til at bestemme koefficienterne til den lineære kombination af Fock-matricer
    Eigen::VectorXd lccoeff; //En vektor, der gemmer på koefficienterne til den lineære kombination af Fock-matricer i DIIS
    DynMatrix lcfmx; //Fock-matrix-gættet i DIIS, der skrives som en lineær kombination af tidligere Fock-matricer

public:
    DIIS(HF* hf, int d);
    ~DIIS();

    ////////////////////////////////////////////////////////////////////////////////
    //PUBLIC FUNKTIONER
    ////////////////////////////////////////////////////////////////////////////////

    /*Gem en Fock-matrix*/
    void savefmx(DynMatrix* mx, int index);

    /*Udregn en fejlmatrix ud fra myhf->fmx, myhf->dmx og myhf->smx*/
    void computeemx(int index);

    /*Lav B-matricen*/
    void computebmx();

    /*Løs B-matricen og bestem koefficienterne til Fock-matricerne i DIIS*/
    void solvebmx();

    /*Bestem Fock-matrix-gættet i DIIS ud fra en lineær kombination af tidligere Fock-matricer*/
    void computelcfmx();

    ////////////////////////////////////////////////////////////////////////////////
    //VENNER
    ////////////////////////////////////////////////////////////////////////////////

    /*Private datamedlemmer fra et DIIS-objekt skal tilgås i HFSCFDIIS-metoden i HF-klassen*/
    friend class HF;
};

#endif