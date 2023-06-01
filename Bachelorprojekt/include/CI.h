/********************************************************
I denne .h-fil er defineret klassen CI, der med udgangspunkt
i en HF-beregning for et givent molekyle, udregner eksitationsenergier
for molekylet. CIS-metoden udregner alle CIS-eksitationsenergierne.
Et objekt af klassen defineres ud fra en pointer til et HF-objekt,
der gemmer på information om det undersøgte molekyle og integraldata.
------------BESKRIVELSE AF VIGTIGE DATAMEDLEMMER------------
myhf: HF-objektet
nosos: Det totale antal spinorbitaler, molekylet i HF-objektet beskrives med
nooccsos: Antallet af okkuperede spinorbitaler i molekylet
novrsos: Antallet af virtuelle spinorbitaler i molekylet
hmx: H-matricen, der skal gemme på forventningsværdier af den elektroniske Hamilton-operator for forskellige kombinationer af eksiterede Slater-determinanter
motmx: Matricen med integralerne for kinetisk energi i MO-basis (i rumorbital basis)
movmx: Matricen med integralerne for kernetiltrækningen i MO-basis (i rumorbital basis)
mofmx: Matricen med Fock-integralerne i MO-basis (i rumorbital basis)
moeritensor: Tensoren med to-elektron-integralerne i MO-basis (i rumorbital basis)
excienergies: Et array med eksitationsenergierne fra CIS-beregningen
------------BESKRIVELSE AF VIGTIGE MEDLEMSFUNKTIONER------------
CIS: Laver CIS-beregningen og bestemmer alle CIS-eksitationsenergierne for molekylet i myhf-objektet.
Udskriver resultaterne til en fil via outfile-variablen.
********************************************************/

#ifndef _CI
#define _CI

#include <fstream>
#include <time.h>
#include "HF.h"
#include "../sources/Eigen/Dense"
#include "../sources/Eigen/Eigenvalues"
#include "../sources/Eigen/Core"
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> DynMatrix;

class CI {
private:
    ////////////////////////////////////////////////////////////////////////////////
    //PRIVATE DATAMEDLEMMER
    ////////////////////////////////////////////////////////////////////////////////

    HF* myhf; //HF-objektet, der gemmer på information om det undersøgte molekyle og alle de relevante integraler
    int nosos, nooccsos, novrsos; //Det totale antal spinorbitaler for molekylet i myhf-objektet, antallet af okkuperede spinorbitaler og antalle af virtuelle spinorbitaler
    DynMatrix hmx; //CIS-matricen, der gemmer på forventningsværdier af den elektroniske Hamilton-operator for forskellige eksiterede Slater-determinanter
    double** motmx; //En matrix, der gemmer på integralerne for kinetisk energi i MO-basis (i rumorbital-basis)
    double** movmx; //En matrix, der gemmer på integralerne for kernetiltrækningen i MO-basis (i rumorbital-basis)
    double** mofmx; //En matrix, der gemmer på Fock-integralerne i MO-basis (i rumorbital-basis)
    double**** moeritensor; //En tensor, der gemmer på to-elektron-integralerne i MO-basis (i rumorbital-basis)
    double* excienergies; //En vektor, der gemmer på eksitationsenergierne fra CIS-beregningen
    ofstream& outfile; //Den fil resultaterne af en beregning udskrives til
    time_t time1; //En variabel til at holde styr på beregningstiden
    
    ////////////////////////////////////////////////////////////////////////////////
    //PRIVATE FUNKTIONER
    ////////////////////////////////////////////////////////////////////////////////

    /*Transformér to-elektron-integralerne fra AO- til MO-basis (rumorbital-basis)*/
    void computemoeritensor();

    /*Transformér integralerne for kinetisk energi fra AO- til MO-basis (rumorbital-basis)*/
    void computemotmx();

    /*Transformér integralerne for kernetiltrækningen fra AO- til MO-basis (rumorbital-basis)*/
    void computemovmx();

    /*Transformér Fock-integralerne fra AO- til MO-basis (rumorbital-basis)*/
    void computemofmx();

public:
    CI(HF* hf, ofstream& ofile);
    ~CI();

    ////////////////////////////////////////////////////////////////////////////////
    //PUBLIC FUNKTIONER
    ////////////////////////////////////////////////////////////////////////////////

    /*Udregn den fulde CIS-matrix til brug i CIS-metoden (og til brug i TDHF)*/
    void computehmx();

    /*Lav CIS-beregningen og bestem eksitationsenergierne for det undersøgte molekyle*/
    void CIS();

    /*Bestem kun singlet-eksitationsenergierne*/
    void CISSINGLETS();

    /*Bestem kun triplet-eksitationsenergierne*/
    void CISTRIPLETS();

    ////////////////////////////////////////////////////////////////////////////////
    //VENNER
    ////////////////////////////////////////////////////////////////////////////////

    /*TDHF-klassen skal have adgang til bl.a. H-matricen fra et CI-objekt, HF-objektet m.m.*/
    friend class TDHF;
};

#endif