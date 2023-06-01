/*
I denne .h-fil er defineret klassen Atom.
Et objekt af klassen (et atom) konstrueres med 
værdier for datamedlemmerne symbol (atomsymbol),
atomicnumber og x, y og z (koordinater)
Atomsymbolet kan angives med setsymbol()
Atomnummeret kan angives med setatomicnumber()
Atomet kan placeres med funktionen setposition()
M.m.
*/

#ifndef _ATOM
#define _ATOM

#include <iostream>
#include <fstream>
#include <cstring>
#include <math.h>
using namespace std;

class Atom {
private:
    ////////////////////////////////////////////////////////////////////////////////
    //PRIVATE DATAMEDLEMMER
    ////////////////////////////////////////////////////////////////////////////////

    char symbol[3]; //Symbolet for atomet
    int atomicnumber; //Atomnummeret Z for atomet
    double x, y, z; //Koordinaterne til atomet

public:
    //Constructor
    Atom(const int r = 1, const char* const s = "H", const double t = 0, const double u = 0, const double v = 0)
    : atomicnumber(r), x(t), y(u), z(v) {strncpy(symbol, s, sizeof(symbol)); symbol[sizeof(symbol)-1] = '\0';}

    ////////////////////////////////////////////////////////////////////////////////
    //PUBLIC SET-FUNKTIONER
    ////////////////////////////////////////////////////////////////////////////////

    //Angiv atomnummeret for atomet med et atomsymbol som input. Returnerer atomnummeret eller 0.
    int setatomicnumber(const char* const);

    //Angiv atomnummeret for atomet ud fra this->symbol
    int setatomicnumber();

    //Angiv symbolet for atomet
    void setsymbol(const char* const s) {strncpy(symbol, s, sizeof(symbol)); symbol[sizeof(symbol)-1] = '\0';}

    //Placér atomet
    void setposition(const double t, const double u, const double v) {x = t; y = u; z = v;}

    void setx(const double t) {x = t;} //Angiv x-koordinatet for atomet
    void sety(const double u) {y = u;} //Angiv y-koordinatet for atomet
    void setz(const double v) {z = v;} //Angiv z-koordinatet for atomet

    ////////////////////////////////////////////////////////////////////////////////
    //PUBLIC PRINTFUNKTIONER
    ////////////////////////////////////////////////////////////////////////////////

    //Print atomnummeret og koordinaterne for atomet til terminalen
    void print() const;

    //Print atomnummeret og koordinaterne for atomet til en fil
    void print(ofstream&) const;

    ////////////////////////////////////////////////////////////////////////////////
    //PUBLIC GET-FUNKTIONER
    ////////////////////////////////////////////////////////////////////////////////

    //Hent x-koordinatet for atomet
    double getx() const {return x;}

    //Hent y-koordinatet for atomet
    double gety() const {return y;}

    //Hent z-koordinatet for atomet
    double getz() const {return z;}

    //Hent atomnummeret for atomet
    int getatomicnumber() const {return atomicnumber;}

    ////////////////////////////////////////////////////////////////////////////////
    //ANDRE PUBLIC FUNKTIONER
    ////////////////////////////////////////////////////////////////////////////////

    //Forskyd atomet
    void displace(const double t, const double u, const double v) {x += t; y += u; z += v;};

    //Omregn koordinaterne til Bohr antaget de er opgivet i Ångstrom
    void converttobohr();

    //Udregn afstanden mellem to atomer
    friend double distance(Atom myatom1, Atom myatom2)
    {
        return sqrt( (myatom1.x - myatom2.x)*(myatom1.x - myatom2.x) + (myatom1.y - myatom2.y)*(myatom1.y - myatom2.y) + (myatom1.z - myatom2.z)*(myatom1.z - myatom2.z) );
    }
};

#endif