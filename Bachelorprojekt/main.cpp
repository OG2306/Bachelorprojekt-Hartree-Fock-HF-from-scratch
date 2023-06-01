#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include "include/TDHF.h"
#include "include/CI.h"
#include "include/DIIS.h"
#include "include/HF.h"
#include "include/Molecule.h"
#include "include/Atom.h"
using namespace std;

int main(int argc, char* argv[])
{
    if (argc != 2) {
        cout << "Fejl: Programmet forventer en filsti til mappen med molekyl- og integraldata som input.\n";
        return 0;
    }
    const char* const integraldir = argv[1];

    char xyzinputfile[200];
    snprintf(xyzinputfile, sizeof(xyzinputfile), "%s%s", integraldir, "/geom.xyz"); //Gem filstien til .xyz-filen
    
    Molecule mymolecule(xyzinputfile); //Lav molekylet
    HF myhf(mymolecule, integraldir); //Lav HF-objektet
    ofstream outfile("beregningsresultater.txt", ios_base::out);

    /*KONVERTERING FRA ÅNGSTROM TIL BOHR HVIS NØDVENDIGT*/
    //mymolecule.convertalltobohr();
    //mymolecule.computencdipolemoment(),
    //mymolecule.computencrepenergy();
    /*----------------------------------*/

    outfile << "-------------------------------------\n";
    outfile << "-----------BACHELORPROJEKT-----------\n";
    outfile << "-------------------------------------\n\n";
    outfile << "Indholdet fra den indlæste .xyz-fil \"" << xyzinputfile << "\":\n\n";
    outfile << "//////////////////////////////////////////////////\n";
    mymolecule.printxyz(outfile);
    outfile << "//////////////////////////////////////////////////\n\n";

    //Energiberegningen
    outfile << "-----------HF-ENERGIBEREGNINGEN-----------\n\n";
    outfile << "Det brugte basissæt er: ";
    myhf.printbasissetname(outfile);
    outfile << "\nAntallet af elektroner i molekylet er: ";
    mymolecule.printnoelec(outfile);
    outfile << "\nAntallet af basisfunktioner er: ";
    myhf.printnobasisfcn(outfile);
    outfile << "\nAntallet af dobbeltokkuperede MO'er er: ";
    mymolecule.printnooccmos(outfile);
    outfile << "\n\n";

    int maxiter = 100;
    double delta1 = 1.0e-12;
    double delta2 = 1.0e-12;
    
    myhf.HFSCF(maxiter, delta1, delta2, outfile); //Uden DIIS
    
    myhf.HFSCFDIIS(maxiter, delta1, delta2, outfile, 8); //Med DIIS

    //Dipolmoment
    
    outfile << "----------------DIPOLMOMENT----------------\n\n";

    dipolemoment(myhf, outfile);

    //MP2
    
    outfile << "--------------------MP2--------------------\n\n";
    MP2(myhf, outfile);
    
    outfile << "-----------------MP2Quick-----------------\n\n";
    MP2Quick(myhf, outfile);
    
    //CI

    outfile << "--------------------CIS--------------------\n\n";
    CI myCI(&myhf, outfile);
    myCI.CIS();

    //TDHF

    //outfile << "-------------------TDHF-------------------\n\n";
    TDHF myTDHF(&myCI, outfile);
    outfile << "--------STATISKE POLARISABILITETER--------\n\n";
    myTDHF.statpol();

    outfile << "-------DYNAMISKE POLARISABILITETER (MED INVERSION)-------\n\n";
    myTDHF.dynpol(0); //Vær opmærksom på at beregningstiden i dynpol, som den bliver printet til outputfilen, vil inkludere beregningstiden i statpol, fordi statpol bliver kaldt ovenfor

    outfile << "-------DYNAMISKE POLARISABILITETER (UDEN INVERSION)-------\n\n";
    myTDHF.dynpol2(0);

    outfile.close();
    return 0;
}