#include "../include/Atom.h"
#include "../include/AtomicNumber.h" //For findatomicnumber()

////////////////////////////////////////////////////////////////////////////////
//PUBLIC SET-FUNKTIONER
////////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------------

int Atom::setatomicnumber(const char* const s)
{
    atomicnumber = findatomicnumber(s);
    return atomicnumber;
}

//--------------------------------------------------------------------------------

int Atom::setatomicnumber()
{
    atomicnumber = findatomicnumber(symbol);
    return atomicnumber;
}

//--------------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
//PUBLIC PRINTFUNKTIONER
////////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------------

void Atom::print() const
{
    int w = 15; //Bredden (width)
    int p = 12; //Præcisionen
    cout.setf(ios_base::fixed, ios_base::floatfield);
    cout << symbol << ' ';
    cout.width(w);
    cout.precision(p);
    cout << x << ' ';
    cout.width(w);
    cout.precision(p);
    cout << y << ' ';
    cout.width(w);
    cout.precision(p);
    cout << z << '\n';
}

//--------------------------------------------------------------------------------

void Atom::print(ofstream& outfile) const
{
    int w = 15; //Bredden (width)
    int p = 12; //Præcisionen
    outfile.setf(ios_base::fixed, ios_base::floatfield);
    outfile << symbol << ' ';
    outfile.width(w);
    outfile.precision(p);
    outfile << x << ' ';
    outfile.width(w);
    outfile.precision(p);
    outfile << y << ' ';
    outfile.width(w);
    outfile.precision(p);
    outfile << z << '\n';
}

//--------------------------------------------------------------------------------

void Atom::converttobohr()
{
    double convf = 1.8897260;
    setx(x*convf);
    sety(y*convf);
    setz(z*convf);
}

//--------------------------------------------------------------------------------