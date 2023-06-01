#include "../include/Molecule.h"
#include "../include/Atom.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <cstring>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
//STATISKE DATAMEDLEMMER
////////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------------

const int Molecule::maxsize = 50;

//--------------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
//CONSTRUCTOR
////////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------------

Molecule::Molecule(const char* const xyzfilename)
{
    char line[100]; //Den linje eller det tekststykke, der læses i .xyz-filen
    noelec = 0; //Initiér noelec


    /*Lav en input-filstrøm til .xyz-filen*/
    ifstream infile(xyzfilename, ios_base::in);
    if (!infile) {
        cout << "Fejl: Filen \"" << xyzfilename << "\" kunne ikke åbnes.\n";
        exit(1);
    }


    /*Hent første linje (antal atomer i molekylet)*/
    infile.getline(line, sizeof(line), '\n');
    if (!(noatoms = atoi(line))) {
        cout << "Fejl: Første linje i filen \"" << xyzfilename << "\" angiver ikke antal atomer.\n";
        exit(1);
    }
    else if (noatoms > maxsize) {
        cout << "Fejl: Det angivne antal atomer, " << noatoms << ", i filen \"" << xyzfilename << "\" overstiger det maksimalt tilladte antal, " << maxsize << ".\n";
        exit(1);
    }


    /*Hent anden linje (kommentarlinjen)*/
    infile.getline(line, sizeof(line), '\n');
    strncpy(comment, line, sizeof(comment));
    comment[sizeof(comment)-1] = '\0';


    /*Gør plads til alle atomerne og atommærkaterne*/
    atom = new Atom [noatoms];
    atomlabel = new char* [noatoms];
    for (int i = 0; i < noatoms; i++) atomlabel[i] = new char [5];


    /*Hent alle atomsymbolerne og -koordinaterne ned*/
    int counter = 0;
    while (counter < noatoms)
    {
        /*Hent atomsymbolet og tjek for eof*/
        infile.getline(line, sizeof(line), ' ');
        if (infile.eof()) {
            cout << "Fejl: Der er uoverensstemmelse mellem det oplyste antal atomer, " << noatoms << ", og det faktiske antal, " << counter << ", i filen \"" << xyzfilename << "\".\n";
            exit(1);
        }


        /*Angiv atomnummeret og atomsymbolet for the givne atom. Hvis atomsymbolet ikke genkendes returneres 0.*/
        if (!atom[counter].setatomicnumber(line)) {
            cout << "Fejl: Atomsymbolet \"" << line << "\" i filen \"" << xyzfilename << "\" kunne ikke genkendes.\n";
            exit(1);
        }
        atom[counter].setsymbol(line);


        /*Læg antal elektroner for atomet til molekylets elektrontal og print atomsymbolet + indeks til atommærkatet*/
        noelec += atom[counter].getatomicnumber();
        snprintf(atomlabel[counter], sizeof(atomlabel[counter]), "%s%d", line, 1);


        /*Ændr indekset på atommærkatet (1, 2, 3 ...), hvis grundstoffet optræder flere gange */
        int newindex = 2; //Angiver det mulige nye indeks
        for (int i = 0; i < counter; i++) {
            if (strcmp(atomlabel[i], atomlabel[counter]) == 0) {
                snprintf(atomlabel[counter], sizeof(atomlabel[counter]), "%s%d", line, newindex); //Print det nye atommærkat
                newindex++; //Forøg newindex med 1 og fortsæt med at tjekke for ens atomlabels
            }
        }

        
        /*Hent atomkoordinaterne*/
        while (infile.peek() == ' ') {
            infile.ignore(); //Fjern eventuelle ekstra mellemrum mellem atomsymbolet og x-koordinatet
        }

        infile.getline(line, sizeof(line), ' '); //Hent x-koordinatet
        atom[counter].setx(atof(line));

        while (infile.peek() == ' ') {
            infile.ignore();
        }

        infile.getline(line, sizeof(line), ' '); //Hent y-koordinatet
        atom[counter].sety(atof(line));

        while (infile.peek() == ' ') {
            infile.ignore();
        }

        infile.getline(line, sizeof(line), '\n'); //Hent z-koordinatet
        atom[counter].setz(atof(line));


        counter++;
    }


    //Tjek for eof (ignorér tom linje)
    if (!infile.eof() && infile.peek() != -1) {
        cout << "Fejl: Det oplyste antal atomer, " << noatoms << ", er mindre end det faktiske antal atomer i filen \"" << xyzfilename << "\".\n";
        exit(1);
    }


    //Angiv antallet af dobbeltokkuperede mo'er for gs af molekylet
    nooccmos = noelec / 2;


    //Udregn kernefrastødningen
    computencrepenergy();


    //Udregn kernedipolmomentet
    computencdipolemoment();


    infile.close();
}

//--------------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
//DESTRUCTOR
////////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------------

Molecule::~Molecule()
{
    for (int i = 0; i < noatoms; i++) delete[] atomlabel[i];
    delete[] atomlabel;
    delete[] atom;
}

//--------------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
//PUBLIC PRINTFUNKTIONER
////////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------------

void Molecule::printxyz() const
{
    cout << noatoms << '\n';
    cout << comment << '\n';
    for (int i = 0; i < noatoms; i++) {
        atom[i].print();
    }
}

//--------------------------------------------------------------------------------

void Molecule::printxyz(ofstream& outfile) const
{
    outfile << noatoms << '\n';
    outfile << comment << '\n';
    for (int i = 0; i < noatoms; i++) {
        atom[i].print(outfile);
    }
}

//--------------------------------------------------------------------------------

void Molecule::printlabel(int n) const
{
    if (n >= noatoms) {
        cout << "Fejl: Argumentet til printlabel(), " << n << ", er \"out of range\".\nAntallet af atomer i molekylet er: " << noatoms << ".\n";
        exit(1);
    }
    cout << atomlabel[n];
}

//--------------------------------------------------------------------------------

void Molecule::printlabel(int n, ofstream& outfile) const
{
    if (n >= noatoms) {
        outfile << "Fejl: Argumentet til printlabel(), " << n << ", er \"out of range\".\nAntallet af atomer i molekylet er: " << noatoms << ".\n";
        exit(1);
    }
    outfile << atomlabel[n];
}

//--------------------------------------------------------------------------------

void Molecule::printdistances() const
{
    int w = 15; //Bredden (width)
    int p = 12; //Præcisionen
    cout.setf(ios_base::fixed, ios_base::floatfield);


    /*Print alle atommærkaterne i øverste linje af tabellen (  C1      C2      H1 osv.)*/
    cout << "  ";
    for (int i = 0; i < noatoms; i++) {
        cout.width(w);
        cout << atomlabel[i];
    }
    cout << '\n';


    /*Print alle kombinationer af afstande mellem atomerne i de resterende linjer*/
    for (int i = 0; i < noatoms; i++) {
        cout << atomlabel[i]; //Print mærkatet for atom i
        for (int j = 0; j < noatoms; j++) {
            cout.width(w);
            cout.precision(p);
            cout << distance(atom[i], atom[j]);
        }
        cout << '\n';
    }
}

//--------------------------------------------------------------------------------

void Molecule::printdistances(ofstream& outfile) const
{
    int w = 15; //Bredden (width)
    int p = 12; //Præcisionen
    outfile.setf(ios_base::fixed, ios_base::floatfield);


    /*Print alle atommærkaterne i øverste linje af tabellen (C1      C2      H1 osv.)*/
    outfile << "  ";
    for (int i = 0; i < noatoms; i++) {
        outfile.width(w);
        outfile << atomlabel[i];
    }
    outfile << '\n';


    /*Print alle kombinationer af afstande mellem atomerne i de resterende linjer*/
    for (int i = 0; i < noatoms; i++) {
        outfile << atomlabel[i]; //Print mærkatet for atom i
        for (int j = 0; j < noatoms; j++) {
            outfile.width(w);
            outfile.precision(p);
            outfile << distance(atom[i], atom[j]);
        }
        outfile << '\n';
    }
}

//--------------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
//PUBLIC COMPUTE-FUNKTIONER
////////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------------

void Molecule::computencrepenergy()
{
    ncrepenergy = 0;
    /*Atomare enheder bruges*/
    for (int i = 0; i < (noatoms - 1); i++) { //Bemærk: i < (noatoms - 1).
        for (int j = i + 1; j < noatoms; j++) { //Bemærk: j initieres til i + 1. Så får vi alle interaktioner én gang.
            ncrepenergy += atom[i].getatomicnumber() * atom[j].getatomicnumber() * (1.0 / distance(atom[i], atom[j]));
        }
    }
}

//--------------------------------------------------------------------------------

/*OBS: Vær opmærksom på, hvor origo er placeret. Origo skal være det samme for udregningen
af elektron- og kernedipolmomentet.*/
void Molecule::computencdipolemoment()
{
    ncdipolemomentx = 0; ncdipolemomenty = 0; ncdipolemomentz = 0;
    for (int i = 0; i < noatoms; i++) {
        int Z = atom[i].getatomicnumber();
        ncdipolemomentx += atom[i].getx() * Z;
        ncdipolemomenty += atom[i].gety() * Z;
        ncdipolemomentz += atom[i].getz() * Z;
    }
}

//--------------------------------------------------------------------------------

void Molecule::convertalltobohr()
{
    for (int i = 0; i < noatoms; i++) {
        atom[i].converttobohr();
    }
}