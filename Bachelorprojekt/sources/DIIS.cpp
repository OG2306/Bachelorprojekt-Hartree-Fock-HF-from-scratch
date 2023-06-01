#include <stdlib.h>
#include <iostream>
#include "../include/DIIS.h"
#include "../include/HF.h"
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"
using namespace std;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> DynMatrix;

////////////////////////////////////////////////////////////////////////////////
//CONSTRUCTOR
////////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------------

DIIS::DIIS(HF* hf, int d)
: myhf(hf), maxdimensions(d), curdimensions(0)
{
    emx = new DynMatrix [maxdimensions];
    fmxhistory = new DynMatrix [maxdimensions];
}

//--------------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
//DESTRUCTOR
////////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------------

DIIS::~DIIS()
{
    delete[] emx;
    delete[] fmxhistory;
}

//--------------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
//PUBLIC FUNKTIONER
////////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------------

void DIIS::savefmx(DynMatrix* mx, int index)
{
    if (index >= maxdimensions) {
        cout << "Fejl: Indekset i DIIS::savefmx er lig med eller overstiger det maksimale antal dimensioner.\n";
        exit(1);
    }
    fmxhistory[index] = *mx;
    if (curdimensions < maxdimensions) curdimensions++; //Opdater dimensionen
}

//--------------------------------------------------------------------------------

void DIIS::computeemx(int index)
{
    if (index >= maxdimensions) {
        cout << "Fejl: Indekset i DIIS::computeemx er lig med eller overstiger det maksimale antal dimensioner.\n";
        exit(1);
    }
    emx[index] = (myhf->fmx * myhf->dmx * myhf->smx) - (myhf->smx * myhf->dmx * myhf->fmx); //Formlen, som den er givet i projektet
}

//--------------------------------------------------------------------------------

void DIIS::computebmx()
{
    bmx = DynMatrix::Zero(curdimensions + 1, curdimensions + 1); //Giv B-matricen den rette størrelse
    for (int i = 0; i < curdimensions; i++) { //Sum over rækkerne - 1 i B-matricen
        for (int j = 0; j <= i; j++) { //Sum over kolonnerne - 1 i B-matricen, men kun ind til (og med) diagonalen
            double dotprod = 0;
            for (int k = 0; k < myhf->nobasisfcn; k++) { //Sum over rækkerne i fejlmatricen / -matricerne
                for (int l = 0; l < myhf->nobasisfcn; l++) { //Sum over kolonnerne i fejlmatricen / -matricerne
                    dotprod += emx[i](k, l) * emx[j](k, l); //Udregn prikproduktet af række k i fejlmatrix i med række k i fejlmatrix j
                }
            }
            bmx(i, j) = dotprod;
            bmx(j, i) = dotprod;
        }
    }
    for (int i = 0; i < curdimensions; i++) {
        //Elementerne i sidste række og sidste kolonne skal være -1 på nær allersidste element
        bmx(curdimensions, i) = -1;
        bmx(i, curdimensions) = -1;
    }
    bmx(curdimensions, curdimensions) = 0;
}

//--------------------------------------------------------------------------------

void DIIS::solvebmx()
{
    Eigen::VectorXd x = Eigen::VectorXd::Zero(curdimensions + 1); //Den vektor, vi løser B for i udtrykket Bc = x
    x(curdimensions) = -1; //Sidste entry i vektoren skal være -1
    lccoeff = bmx.lu().solve(x); //Bestem koefficienterne til den lineære kombination af Fock-matricer
}

//--------------------------------------------------------------------------------

void DIIS::computelcfmx()
{
    lcfmx = DynMatrix::Zero(myhf->nobasisfcn, myhf->nobasisfcn);
    for (int i = 0; i < curdimensions; i++) {
        lcfmx += lccoeff(i) * fmxhistory[i]; //Skriv Fock-matrix-gættet som en lineær kombination af tidligere Fock-matricer
    }
}

//--------------------------------------------------------------------------------