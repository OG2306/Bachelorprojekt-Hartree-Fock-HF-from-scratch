/*
I denne .h-fil er defineret klassen Molecule.
Et objekt af klassen (et molekyle) konstrueres ud fra en .xyz-fil.
Datamedlemmet atom er en liste af objekter
af typen Atom, der indeholder oplysninger om alle
atomerne i molekylet (koordinater, atomsymboler m.m.).
Datamedlemmet noatoms angiver det totale antal atomer i molekylet.
Datamedlemmet noelec angiver det totale antal elektroner i molekylet.
Datamedlemmet nooccmos angiver antallet af dobbeltokkuperede molekylorbitaler for grundtilstanden af molekylet.
Datamedlemmet ncrepenergy angiver frastødningen mellem atomkernerne i molekylet (a.u.)
Datamedlemmerne ncdipolemomentx, -y og -z angiver x-, y- og z-komponenterne af molekylets kernedipolmoment
M.m.
*/

#ifndef _MOLECULE
#define _MOLECULE

#include "Atom.h"
#include <fstream>

class Molecule {
private:
    ////////////////////////////////////////////////////////////////////////////////
    //PRIVATE DATAMEDLEMMER
    ////////////////////////////////////////////////////////////////////////////////

    int noatoms; //Antal atomer i molekylet
    char comment[100]; //Kommentaren fra .xyz-filen. Max 100 tegn.
    int noelec; //Antal elektroner i molekylet
    int nooccmos; //Antal dobbeltokkuperede MO'er for gs af molekylet
    char** atomlabel; //Navne på atomerne i molekylet (C1, C2, H1 og lignende)
    double ncrepenergy; //Frastødningen mellem atomkernerne i molekylet (a.u.)
    double ncdipolemomentx; //x-komponenten af kernedipolmomentet i molekylet (a.u.)
    double ncdipolemomenty; //y-komponenten af kernedipolmomentet i molekylet (a.u.)
    double ncdipolemomentz; //z-komponenten af kernedipolmomentet i  molekylet (a.u.)
    static const int maxsize; //Det maksimale antal tilladte atomer i molekylet

public:
    Atom* atom; //Atomerne i molekylet. OBS er blevet gjort public pga. HF::setnobasisfcn()

    /*Et objekt af klassen konstrueres ud fra en .xyz-fil*/
    Molecule(const char* const xyzfilename);
    ~Molecule(); //Destructor

    ////////////////////////////////////////////////////////////////////////////////
    //PUBLIC PRINTFUNKTIONER
    ////////////////////////////////////////////////////////////////////////////////

    //Print det læste indhold fra .xyz-filen til terminalen
    void printxyz() const;

    //Print det læste indhold fra .xyz-filen til en fil
    void printxyz(ofstream&) const;

    //Print atommærkatet for atom n i molekylet til terminalen
    void printlabel(int n) const;

    //Print atommærkatet for atom n i molekylet til en fil
    void printlabel(int n, ofstream&) const;

    //Print antallet af elektroner i molekylet til terminalen
    void printnoelec() const {cout << noelec;}

    //Print antallet af elektroner i molekylet til en fil
    void printnoelec(ofstream& outfile) const {outfile << noelec;}

    //Print afstanden mellem alle atomer i molekylet i en tabel i terminalen
    void printdistances() const;

    //Print afstanden mellem alle atomer i molekylet i en tabel i en fil
    void printdistances(ofstream&) const;

    //Print kernefrastødningen mellem atomkernerne i molekylet til terminalen
    void printncrepenergy() const {cout << ncrepenergy;}

    //Print kernefrastødningen mellem atomkernerne i molekylet til en fil
    void printncrepenergy(ofstream& outfile) const {outfile << ncrepenergy;}

    //Print antallet af okkuperede molekylorbitaler i molekylet til terminalen
    void printnooccmos() const {cout << nooccmos;}

    //Print antallet af okkuperede molekylorbitaler i molekylet til en fil
    void printnooccmos(ofstream& outfile) const {outfile << nooccmos;}

    ////////////////////////////////////////////////////////////////////////////////
    //PUBLIC GET-FUNKTIONER
    ////////////////////////////////////////////////////////////////////////////////

    //Returnerer antallet af elektroner i molekylet
    int getnoelec() const {return noelec;}

    //Returnerner antallet af dobbeltokkuperede molekylorbitaler for gs af molekylet
    int getnooccmos() const {return nooccmos;}

    //Returnerer antallet af atomer i molekylet
    int getnoatoms() const {return noatoms;}
    
    //Returnerer kernefrastødningen mellem atomkernerne
    double getncrepenergy() const {return ncrepenergy;}

    //Returnerer x-komponenten af kernedipolmomentet
    double getncdipolemomentx() const {return ncdipolemomentx;}

    //Returnerer y-komponenten af kernedipolmomentet
    double getncdipolemomenty() const {return ncdipolemomenty;}

    //Returnerer z-komponenten af kernedipolmomentet
    double getncdipolemomentz() const {return ncdipolemomentz;}

    ////////////////////////////////////////////////////////////////////////////////
    //PUBLIC COMPUTE-FUNKTIONER
    ////////////////////////////////////////////////////////////////////////////////

    //Udregn frastødningen mellem atomkernene i molekylet
    void computencrepenergy();

    //Udregn atomkernedipolmomentet
    void computencdipolemoment();

    ////////////////////////////////////////////////////////////////////////////////
    //ANDRE
    ////////////////////////////////////////////////////////////////////////////////

    //Ændr atomkoordinater fra Ångstrom til Bohr
    void convertalltobohr();
};

#endif