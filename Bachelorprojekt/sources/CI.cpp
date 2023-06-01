#include "../include/CI.h"
#include "../include/HF.h"
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"
#include <algorithm>
#include <math.h>
#include <time.h>
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> DynMatrix;

////////////////////////////////////////////////////////////////////////////////
//CONSTRUCTOR
////////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------------

CI::CI(HF* hf, ofstream& ofile)
: myhf(hf), outfile(ofile)
{
    time1 = time(NULL);
    nosos = myhf->nobasisfcn * 2;
    nooccsos = myhf->mymolecule.getnooccmos() * 2;
    novrsos = nosos - nooccsos;
    computemoeritensor();
    computemotmx();
    computemovmx();
    computemofmx();
}

//--------------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
//DESTRUCTOR
////////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------------

CI::~CI()
{
    for (int i = 0; i < myhf->nobasisfcn; i++) {delete[] motmx[i]; delete[] movmx[i]; delete[] mofmx[i];}
    delete[] motmx;
    delete[] movmx;
    delete[] mofmx;

    for (int i = 0; i < myhf->nobasisfcn; i++) { //Loop over rækker i den ydre 2D-matrix 
        for (int j = 0; j < myhf->nobasisfcn; j++) { //Loop over kolonner i den ydre 2D-matrix
            for (int k = 0; k < myhf->nobasisfcn; k++) { //Loop over rækker i de indre 2D-matricer
                delete[] moeritensor[i][j][k];
            }
            delete[] moeritensor[i][j];
        }
        delete[] moeritensor[i];
    }
    delete[] moeritensor;
}

//--------------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
//PRIVATE FUNKTIONER
////////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------------

void CI::computemoeritensor()
{
    //SAMME KODE SOM MP2Quick
    //Gør variablerne, der skal gemme to-elektron-integralerne i MO-basis / delvis MO-basis, klar
    /********************************************************************/
    moeritensor = new double*** [myhf->nobasisfcn]; //Variablen, der skal gemme to-elektron-integralerne i MO-basis
    double**** moeritensortmp = new double*** [myhf->nobasisfcn]; //En variabel, der skal gemme på midlertidige to-elektron-integraler under AO->MO-transformationen
    for (int i = 0; i < myhf->nobasisfcn; i++) { //Loop over rækker i den ydre 2D-matrix
        moeritensor[i] = new double** [myhf->nobasisfcn];
        moeritensortmp[i] = new double** [myhf->nobasisfcn];
        for (int j = 0; j < myhf->nobasisfcn; j++) { //Loop over kolonner i den ydre 2D-matrix
            moeritensor[i][j] = new double* [myhf->nobasisfcn];
            moeritensortmp[i][j] = new double* [myhf->nobasisfcn];
            for (int k = 0; k < myhf->nobasisfcn; k++) { //Loop over rækker i de indre 2D-matricer
                moeritensor[i][j][k] = new double [myhf->nobasisfcn];
                moeritensortmp[i][j][k] = new double [myhf->nobasisfcn];
                for (int l = 0; l < myhf->nobasisfcn; l++) {
                    moeritensor[i][j][k][l] = 0; //Initiér alle indekser til 0
                    moeritensortmp[i][j][k][l] = 0;
                }
            }
        }
    }
    /********************************************************************/

    //Transformer to-elektron-integralerne i AO-basis (eritensor) til MO-basis (moeritensor)
    /********************************************************************/
    //Transformer første indeks (ij|kl) -> (aj|kl)
    for (int a = 0; a < myhf->nobasisfcn; a++) { //Sum over MO'er (rumorbitaler) (okk. + vir.)
        for (int i = 0; i < myhf->nobasisfcn; i++) { //Sum over basisfunktioner i (ij|kl)
            for (int j = 0; j < myhf->nobasisfcn; j++) { //Sum over basisfunktioner i (ij|kl)
                for (int k = 0; k < myhf->nobasisfcn; k++) { //Sum over basisfunktioner i (ij|kl)
                    for (int l = 0; l < myhf->nobasisfcn; l++) { //Sum over basisfunktioner i (ij|kl)
                        moeritensortmp[a][j][k][l] += myhf->cmx(i, a) * myhf->eritensor[i][j][k][l];
                    }
                }
            }
        }
    }
    //Transformer andet indeks (aj|kl) -> (ab|kl)
    for (int b = 0; b < myhf->nobasisfcn; b++) { //Sum over MO'er (rumorbitaler) (okk. + vir.)
        for (int a = 0; a < myhf->nobasisfcn; a++) { //Sum over MO'er i (aj|kl)
            for (int j = 0; j < myhf->nobasisfcn; j++) { //Sum over basisfunktioner i (aj|kl)
                for (int k = 0; k < myhf->nobasisfcn; k++) { //Sum over basisfunktioner i (aj|kl)
                    for (int l = 0; l < myhf->nobasisfcn; l++) { //Sum over basisfunktioner i (aj|kl)
                        moeritensor[a][b][k][l] += myhf->cmx(j, b) * moeritensortmp[a][j][k][l];
                    }
                }
            }
        }
    }
    //Reset moeritensortmp
    for (int i = 0; i < myhf->nobasisfcn; i++) {
        for (int j = 0; j < myhf->nobasisfcn; j++) {
            for (int k = 0; k < myhf->nobasisfcn; k++) {
                for (int l = 0; l < myhf->nobasisfcn; l++) {
                    moeritensortmp[i][j][k][l] = 0;
                }
            }
        }
    }
    //Transformer tredje indeks (ab|kl) -> (ab|cl)
    for (int c = 0; c < myhf->nobasisfcn; c++) { //Sum over MO'er (rumorbitaler) (okk. + vir.)
        for (int a = 0; a < myhf->nobasisfcn; a++) { //Sum over MO'er i (ab|kl)
            for (int b = 0; b < myhf->nobasisfcn; b++) { //Sum over MO'er i (ab|kl)
                for (int k = 0; k < myhf->nobasisfcn; k++) { //Sum over basisfunktioner i (ab|kl)
                    for (int l = 0; l < myhf->nobasisfcn; l++) { //Sum over basisfunktioner i (ab|kl)
                        moeritensortmp[a][b][c][l] += myhf->cmx(k, c) * moeritensor[a][b][k][l];
                    }
                }
            }
        }
    }
    //Reset moeritensor
    for (int i = 0; i < myhf->nobasisfcn; i++) {
        for (int j = 0; j < myhf->nobasisfcn; j++) {
            for (int k = 0; k < myhf->nobasisfcn; k++) {
                for (int l = 0; l < myhf->nobasisfcn; l++) {
                    moeritensor[i][j][k][l] = 0;
                }
            }
        }
    }
    //Transformer fjerde indeks (ab|cl) -> (ab|cd)
    for (int d = 0; d < myhf->nobasisfcn; d++) { //Sum over MO'er (rumorbitaler) (okk. + vir.)
        for (int a = 0; a < myhf->nobasisfcn; a++) { //Sum over MO'er i (ab|cl)
            for (int b = 0; b < myhf->nobasisfcn; b++) { //Sum over MO'er i (ab|cl)
                for (int c = 0; c < myhf->nobasisfcn; c++) { //Sum over MO'er i (ab|cl)
                    for (int l = 0; l < myhf->nobasisfcn; l++) { //Sum over basisfunktioner i (ab|cl)
                        moeritensor[a][b][c][d] += myhf->cmx(l, d) * moeritensortmp[a][b][c][l];
                    }
                }
            }
        }
    }
    /********************************************************************/

    //Slet moeritensortmp
    /********************************************************************/
    for (int i = 0; i < myhf->nobasisfcn; i++) { //Loop over rækker i den ydre 2D-matrix 
        for (int j = 0; j < myhf->nobasisfcn; j++) { //Loop over kolonner i den ydre 2D-matrix
            for (int k = 0; k < myhf->nobasisfcn; k++) { //Loop over rækker i de indre 2D-matricer
                delete[] moeritensortmp[i][j][k];
            }
            delete[] moeritensortmp[i][j];
        }
        delete[] moeritensortmp[i];
    }
    delete[] moeritensortmp;
    /********************************************************************/
}

//--------------------------------------------------------------------------------

void CI::computemotmx()
{
    motmx = new double* [myhf->nobasisfcn];
    for (int i = 0; i < myhf->nobasisfcn; i++) {
        motmx[i] = new double [myhf->nobasisfcn];
        for (int j = 0; j < myhf->nobasisfcn; j++) {
            motmx[i][j] = 0; //Initiér indekser til 0
        }
    }

    for (int a = 0; a < myhf->nobasisfcn; a++) { //Sum over MO'er
        for (int b = 0; b < myhf->nobasisfcn; b++) { //Sum over MO'er
            for (int i = 0; i < myhf->nobasisfcn; i++) { //Sum over AO'er
                for (int j = 0; j < myhf->nobasisfcn; j++) { //Sum over AO'er
                    motmx[a][b] += myhf->cmx(i, a) * myhf->cmx(j, b) * myhf->tmx(i, j);
                }
            }
        }
    }
}

//--------------------------------------------------------------------------------

void CI::computemovmx()
{
    movmx = new double* [myhf->nobasisfcn];
    for (int i = 0; i < myhf->nobasisfcn; i++) {
        movmx[i] = new double [myhf->nobasisfcn];
        for (int j = 0; j < myhf->nobasisfcn; j++) {
            movmx[i][j] = 0; //Initiér indekser til 0
        }
    }

    for (int a = 0; a < myhf->nobasisfcn; a++) { //Sum over MO'er
        for (int b = 0; b < myhf->nobasisfcn; b++) { //Sum over MO'er
            for (int i = 0; i < myhf->nobasisfcn; i++) { //Sum over AO'er
                for (int j = 0; j < myhf->nobasisfcn; j++) { //Sum over AO'er
                    movmx[a][b] += myhf->cmx(i, a) * myhf->cmx(j, b) * myhf->vmx(i, j);
                }
            }
        }
    }
}

//--------------------------------------------------------------------------------

void CI::computemofmx()
{
    mofmx = new double* [myhf->nobasisfcn];
    for (int i = 0; i < myhf->nobasisfcn; i++) {
        mofmx[i] = new double [myhf->nobasisfcn];
        for (int j = 0; j < myhf->nobasisfcn; j++) {
            mofmx[i][j] = 0; //Initiér indekser til 0
        }
    }
    
    for (int i = 0; i < myhf->nobasisfcn; i++) {
        for (int j = 0; j < myhf->nobasisfcn; j++) {
            mofmx[i][j] += motmx[i][j] + movmx[i][j];
            for (int k = 0; k < myhf->mymolecule.getnooccmos(); k++) {
                mofmx[i][j] += 2 * moeritensor[i][j][k][k] - moeritensor[i][k][j][k];
            }
        }
    }
}

//--------------------------------------------------------------------------------

void CI::computehmx()
{
    hmx = DynMatrix::Zero(nooccsos * novrsos, nooccsos * novrsos); //Giv H-matricen de rette dimensioner (antallet af mulige enkelteksitationer ganget med antallet af mulige enkelteksitationer)
    int aspin, rspin, bspin, sspin; //Variabler, der holder styr på spinnet af elektronerne
    int index1 = -1, index2 = -1; //Variabler, der holder styr på hvilket entry i H-matricen, der skal tilgås

    //Dan H-matricen
    /********************************************************************/
    for (int a = 0; a < nooccsos; a++) { //Vi eksiterer en elektron fra orbital a til...
        for (int r = nooccsos; r < nosos; r++) { //... orbital r
            index1++;
            for (int b = 0; b < nooccsos; b++) { //Vi eksiterer en elektron fra orbital b til...
                for (int s = nooccsos; s < nosos; s++) { //... orbital s
                    index2++;
                    aspin = a % 2; rspin = r % 2; bspin = b % 2; sspin = s % 2; //Eksempel: a = 2 -> aspin = 0, b = 3 -> bspin = 1. Spinnet skifter mellem "0" og "1"

                    /*Bemærk: Når integralerne skal tilgås, divideres der med 2 (heltalsdivision), fordi integralerne er gemt i rumorbital-basis og ikke spinorbital-basis.
                    Når vi er interesserede i et integral, der vedrører spinorbitalerne i og j (f.eks. <i|f|j>), kan vi finde det i rumorbital-basis som integralet, der vedrører
                    rumorbitalerne i/2 og j/2. Variablerne aspin, bspin osv. tjekker om spinnet matcher mellem de spinorbitaler, vi gerne vil have forventningsværdien af.*/
                    if (a == b && rspin == sspin) hmx(index1, index2) += mofmx[r/2][s/2];
                    if (r == s && aspin == bspin) hmx(index1, index2) -= mofmx[a/2][b/2];
                    if (aspin == rspin && bspin == sspin) hmx(index1, index2) += moeritensor[a/2][r/2][b/2][s/2];
                    if (aspin == bspin && rspin == sspin) hmx(index1, index2) -= moeritensor[a/2][b/2][r/2][s/2];
                    if (fabs(hmx(index1, index2)) < 1.0e-12) hmx(index1, index2) = 0;
                }
            }
            index2 = -1;
        }
    }
    /********************************************************************/
}

////////////////////////////////////////////////////////////////////////////////
//PUBLIC FUNKTIONER
////////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------------

void CI::CIS()
{
    computehmx();

    //Diagonalisér H-matricen og bestem eksitationsenergierne
    /********************************************************************/
    Eigen::SelfAdjointEigenSolver<DynMatrix> solvedhmx(hmx);
    excienergies = new double [nooccsos * novrsos];
    for (int i = 0; i < nooccsos * novrsos; i++) excienergies[i] = solvedhmx.eigenvalues()(i);
    sort(excienergies, &(excienergies[nooccsos * novrsos]));
    /********************************************************************/
    time_t time2 = time(NULL);

    //Output resultaterne
    /********************************************************************/
    outfile << "Eksitationsenergier:\n\n";
    outfile.width(20);
    outfile << "Hartrees";
    outfile.width(18);
    outfile << "eV\n";
    for (int i = 0; i < nooccsos * novrsos; i++) {
        outfile.setf(ios_base::left, ios_base::adjustfield);
        outfile.width(3);
        outfile << i + 1;
        outfile.setf(ios_base::right, ios_base::adjustfield);
        outfile.precision(12);
        outfile.width(20);
        outfile << excienergies[i];
        outfile.width(20);
        outfile << excienergies[i] * 27.2114 << '\n'; //Konvertering til eV
    }
    outfile << "\n\n";
    outfile.precision(0);
    outfile << "**Beregningstid: Programmet var " << difftime(time2, time1) << " sekunder om at lave CIS-beregningen (inklusiv transformering af integralerne).\n\n";
    /********************************************************************/
    delete[] excienergies;
}

//--------------------------------------------------------------------------------

void CI::CISSINGLETS() //OBS: Eigen giver problemer
{
    int hmxdim = myhf->mymolecule.getnooccmos() * (myhf->nobasisfcn - myhf->mymolecule.getnooccmos());
    hmx = DynMatrix::Zero(hmxdim, hmxdim);
    int index1 = -1, index2 = -1;
    for (int a = 0; a < myhf->mymolecule.getnooccmos(); a++) {
        for (int i = myhf->mymolecule.getnooccmos(); i < myhf->nobasisfcn; i++) {
            index1++;
            for (int b = 0; b < myhf->mymolecule.getnooccmos(); b++) {
                for (int j = myhf->mymolecule.getnooccmos(); j < myhf->nobasisfcn; j++) {
                    index2++;
                    if (a == b) hmx(index1, index2) += mofmx[i][j];
                    if (i == j) hmx(index1, index2) -= mofmx[a][b];
                    hmx(index1, index2) += 2*moeritensor[a][i][j][b];
                    hmx(index1, index2) -= moeritensor[a][b][i][j];
                }
            }
            index2 = -1;
        }
    }
    //Diagonalisér H-matricen og bestem singleteksitationsenergierne
    /********************************************************************/
    Eigen::SelfAdjointEigenSolver<DynMatrix> solvedhmx(hmx);
    double* singletexcienergies = new double [hmxdim];
    for (int i = 0; i < hmxdim; i++) singletexcienergies[i] = solvedhmx.eigenvalues()(i);
    sort(singletexcienergies, &(singletexcienergies[hmxdim]));
    /********************************************************************/

    //Output resultaterne
    /********************************************************************/
    outfile << "Singlet Eksitationsenergier:\n\n";
    outfile.width(20);
    outfile << "Hartrees";
    outfile.width(18);
    outfile << "eV\n";
    for (int i = 0; i < hmxdim; i++) {
        outfile.setf(ios_base::left, ios_base::adjustfield);
        outfile.width(3);
        outfile << i + 1;
        outfile.setf(ios_base::right, ios_base::adjustfield);
        outfile.precision(12);
        outfile.width(20);
        outfile << singletexcienergies[i];
        outfile.width(20);
        outfile << singletexcienergies[i] * 27.2114 << '\n'; //Konvertering til eV
    }
    outfile << '\n';
    /********************************************************************/
    delete[] singletexcienergies;
}

//--------------------------------------------------------------------------------

void CI::CISTRIPLETS() //OBS: Eigen giver problemer
{
    int hmxdim = myhf->mymolecule.getnooccmos() * (myhf->nobasisfcn - myhf->mymolecule.getnooccmos());
    hmx = DynMatrix::Zero(hmxdim, hmxdim);
    int index1 = -1, index2 = -1;
    for (int a = 0; a < myhf->mymolecule.getnooccmos(); a++) {
        for (int i = myhf->mymolecule.getnooccmos(); i < myhf->nobasisfcn; i++) {
            index1++;
            for (int b = 0; b < myhf->mymolecule.getnooccmos(); b++) {
                for (int j = myhf->mymolecule.getnooccmos(); j < myhf->nobasisfcn; j++) {
                    index2++;
                    if (a == b) hmx(index1, index2) += mofmx[i][j];
                    if (i == j) hmx(index1, index2) -= mofmx[a][b];
                    hmx(index1, index2) -= moeritensor[a][b][i][j];
                }
            }
            index2 = -1;
        }
    }
    //Diagonalisér H-matricen og bestem tripleteksitationsenergierne
    /********************************************************************/
    Eigen::SelfAdjointEigenSolver<DynMatrix> solvedhmx(hmx);
    double* tripletexcienergies = new double [hmxdim];
    for (int i = 0; i < hmxdim; i++) tripletexcienergies[i] = solvedhmx.eigenvalues()(i);
    sort(tripletexcienergies, &(tripletexcienergies[hmxdim]));
    /********************************************************************/

    //Output resultaterne
    /********************************************************************/
    outfile << "Triplet Eksitationsenergier:\n\n";
    outfile.width(20);
    outfile << "Hartrees";
    outfile.width(18);
    outfile << "eV\n";
    for (int i = 0; i < hmxdim; i++) {
        outfile.setf(ios_base::left, ios_base::adjustfield);
        outfile.width(3);
        outfile << i + 1;
        outfile.setf(ios_base::right, ios_base::adjustfield);
        outfile.precision(12);
        outfile.width(20);
        outfile << tripletexcienergies[i];
        outfile.width(20);
        outfile << tripletexcienergies[i] * 27.2114 << '\n'; //Konvertering til eV
    }
    outfile << '\n';
    /********************************************************************/
    delete[] tripletexcienergies;
}

//--------------------------------------------------------------------------------