/********************************************************
I denne .h-fil er defineret klassen TDHF, der udfører
molekylberegninger baseret på time-dependent Hartree-Fock.
Et objekt af klassen konstrueres ud fra en pointer til et
CI-objekt og en reference til et ofstream-objekt.
I CI-objektet er gemt oplysninger om det undersøgte molekyle
og relevante integraler, der bruges i TDHF-beregningerne.
------------BESKRIVELSE AF VIGTIGE DATAMEDLEMMER------------
amx: Matricen, der betegnes A i TDHF. Den er givet i ligning 37 i Toulouse og projekt 12 i Crawford. A-matricen svarer til CIS-matricen. Den hentes fra CI-objektet i TDHF-klassen
bmx: Matricen, der betegnes B i TDHF. Den er givet i ligning 38 i Toulouse og projekt 12 i Crawford. Den udregnes vha. de standard to-elektron-integraler, som er gemt i MO-basis i CI-objektet
abbamx: Matricen der laves ud fra A og B med A i øverste venstre kvadrant, B i øverste højre, B* i nederste venstre og A* i nederste højre.
Bemærk at elementerne i A og B ikke er komplekse (i hvert fald ikke i de beregninger, jeg foretager), og derfor er det ikke nødvendigt at kompleks-konjugere.
trsdipx, trsdipy og trsdipz: Vektorer med overgangsdipolintegralerne for de tre retninger. Integralerne er i SO-basis og sorteret først med ændring i virtuelle orbitaler
og så med ændring i okk. orbitaler. Dvs. [<a|µ|r>, <a|µ|s>, ..., <b|µ|r>, <b|µ|s>, ...], hvor a og b er okkuperede spinorbitaler og r og s er virtuelle spinorbitaler.
trsdip: Et array af tre pointere, der peger på trsdipx, trsdipy og trsdipz. Integralerne tilgås som hhv. trsdip[0], trsdip[1] og trsdip[2].
------------BESKRIVELSE AF VIGTIGE MEDLEMSFUNKTIONER------------
statpol: Udregner de statiske polarisabiliteter af det undersøgte molekyle
dynpol: Udregner de dynamiske polarisabiliteter af det undersøgte molekyle (med inversion). Metoden tager en frekvens som argument.
dynpol: Udregner de dynamiske polarisabiliteter af det undersøgte molekyle (uden inversion). Metoden tager en frekvens som argument.
excitations: Udregner det, der svarer til CIS-eksitationsenergierne, men vha. TDHF. Udregner både eksitationsenergier og deeksitationsenergier
excitationsquicker: Udregner eksitationsenergierne som excitations-metoden, men ikke deeksitationsenergierne. Matricen, der skal diagonaliseres, er derfor mindre
********************************************************/

#ifndef _TDHF
#define _TDHF

#include "CI.h"
#include "../sources/Eigen/Dense"
#include "../sources/Eigen/Eigenvalues"
#include "../sources/Eigen/Core"
#include <valarray>
#include <fstream>
#include <time.h>
using namespace std;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> DynMatrix;

class TDHF {
private:
    ////////////////////////////////////////////////////////////////////////////////
    //PRIVATE DATAMEDLEMMER
    ////////////////////////////////////////////////////////////////////////////////

    int nosinglex; //Antallet af mulige enkelteksitationer (givet ud fra antallet af elektroner i det undersøgte molekyle og basissætstørrelsen)
    CI* myci; //En pointer til et CI-objekt, der bl.a. gemmer på CIS-matricen
    DynMatrix& amx; //Matricen, vi kalder A i TDHF (svarer til CIS-matricen)
    DynMatrix bmx; //Matricen, vi kalder B i TDHF
    DynMatrix abbamx; //AB(B*)(A*)-matricen, der bruges i TDHF
    valarray<double> trsdipx, trsdipy, trsdipz; //Vektorerne med overgangsdipolintegraler. Vektorerne optræder i ligning 91 i Toulouse. OBS: SO-BASIS
    valarray<double>* trsdip[3]; //Et array af pointere, der skal pege på overgangsdipolintegralerne. trsdip[0] = trsdipx, trsdip[1] = trsdipy og trsdip[2] = trsdipz
    ofstream& outfile; //Den fil resultatet af beregningerne udskrives til
    time_t time1; //En variabel til at holde styr på beregningstiden

    ////////////////////////////////////////////////////////////////////////////////
    //PRIVATE FUNKTIONER
    ////////////////////////////////////////////////////////////////////////////////

    /*Lav B-matricen i TDHF*/
    void makebmx();

    /*Lav AB(B*)(A*)-matricen*/
    void makeabbamx();

    /*Lav AO- --> SO-transformationen af overgangsdipolintegralerne*/
    void computetrsdip();

public:
    TDHF(CI*, ofstream& ofile);

    ////////////////////////////////////////////////////////////////////////////////
    //PUBLIC FUNKTIONER
    ////////////////////////////////////////////////////////////////////////////////

    /*Udregn de statiske polarisabiliteter for det undersøgte molekyle*/
    void statpol();

    /*Udregn de dynamiske polarisabiliteter for det undersøgte molekyle*/
    void dynpol(double);

    /*Udregn de dynamiske polarisabiliteter for det undersøgte molekyle uden inversion*/
    void dynpol2(double);

    /*Udregn eksitationsenergierne for det undersøgte molekyle (inklusiv deeksitationsenergierne)*/
    void excitations();

    /*Udregn eksitationsenergerne for det undersøgte molekyle (EKSKLUSIV deeksitationsenergierne)*/
    void excitationsquicker();
};

#endif