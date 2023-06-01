#include "../include/HF.h"
#include "../include/CI.h"
#include "../include/TDHF.h"
#include "../sources/Eigen/Dense"
#include "../sources/Eigen/Eigenvalues"
#include "../sources/Eigen/Core"
#include <math.h>
#include <fstream>
#include <valarray>
#include <time.h>
#include <iostream>
#include <fstream>
#include <assert.h>
using namespace std;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> DynMatrix;

////////////////////////////////////////////////////////////////////////////////
// CONSTRUCTOR
////////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------------

TDHF::TDHF(CI *ciobject, ofstream &ofile)
: myci(ciobject), outfile(ofile), amx(ciobject->hmx), nosinglex(ciobject->nooccsos * ciobject->novrsos), trsdipx(nosinglex), trsdipy(nosinglex), trsdipz(nosinglex)
{
    time1 = time(NULL);
    makebmx();
    makeabbamx();
    computetrsdip();
}

//--------------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
// PRIVATE FUNKTIONER
////////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------------

void TDHF::makebmx()
{
    bmx = DynMatrix::Zero(nosinglex, nosinglex); // Giv B-matricen de rette dimensioner
    int aspin, rspin, bspin, sspin;              // Variabler til at holde styr på spin
    int index1 = -1, index2 = -1;                // Variabler til at holde styr på indekserne i B-matricen
    for (int a = 0; a < myci->nooccsos; a++)
    { // Sum over okk. SO'er
        aspin = a % 2;
        for (int r = myci->nooccsos; r < myci->nosos; r++)
        { // Sum over vir. SO'er
            rspin = r % 2;
            index1++;
            for (int b = 0; b < myci->nooccsos; b++)
            { // Sum over okk. SO'er
                bspin = b % 2;
                for (int s = myci->nooccsos; s < myci->nosos; s++)
                { // Sum over vir. SO'er
                    sspin = s % 2;
                    index2++;
                    if (aspin == rspin && bspin == sspin)
                    {
                        bmx(index1, index2) += myci->moeritensor[a/2][r/2][b/2][s/2]; // OBS: integralerne er gemt i MO-basis!
                    }
                    if (aspin == sspin && bspin == rspin)
                    {
                        bmx(index1, index2) -= myci->moeritensor[a/2][s/2][b/2][r/2];
                    }
                }
            }
            index2 = -1;
        }
    }
}

//--------------------------------------------------------------------------------

void TDHF::makeabbamx()
{
    abbamx = DynMatrix::Zero(nosinglex * 2, nosinglex * 2); // Giv AB(B*)(A*)-matricen de rette dimensioner
    for (int i = 0; i < nosinglex * 2; i++)
    {
        for (int j = 0; j < nosinglex * 2; j++)
        {
            if (i < nosinglex && j < nosinglex)
            {                             // Hvis vi er i øverste venstre kvadrant
                abbamx(i, j) = amx(i, j); // Sæt lig med indekserne i A-matricen
            }
            else if (i < nosinglex && j >= nosinglex)
            {                                         // Hvis vi er i øverste højre kvadrant
                abbamx(i, j) = bmx(i, j % nosinglex); // Sæt lig med indekserne i B-matricen
            }
            else if (i >= nosinglex && j < nosinglex)
            {                                         // Hvis vi er i nederste venstre kvadrant
                abbamx(i, j) = bmx(i % nosinglex, j); // Vi arbejder med reelle basisfunktioner og B = B*
            }
            else
            {                                                     // Hvis vi er i nederste højre kvadrant
                abbamx(i, j) = amx(i % nosinglex, j % nosinglex); // Vi arbejder med reelle basisfunktioner og A = A*
            }
        }
    }
}

//--------------------------------------------------------------------------------

void TDHF::computetrsdip()
{
    myci->myhf->loaddipoleintegrals(); // Load dipolintegralerne i AO-basis fra datafilerne
    int aspin, rspin;                  // Variabler, der holder styr på spin
    int index = -1;                    // En variabel, der holder styr på indekset i trsdip-variablerne
    for (int a = 0; a < myci->nooccsos; a++)
    { // Sum over okkuperede SO'er
        aspin = a % 2;
        for (int r = myci->nooccsos; r < myci->nosos; r++)
        { // Sum over virtuelle SO'er
            rspin = r % 2;
            index++;
            for (int i = 0; i < myci->myhf->nobasisfcn; i++)
            { // Sum over AO'er
                for (int j = 0; j < myci->myhf->nobasisfcn; j++)
                { // Sum over AO'er
                    if (aspin == rspin)
                    {
                        trsdipx[index] += myci->myhf->cmx(i, a / 2) * myci->myhf->cmx(j, r / 2) * myci->myhf->muxmx(i, j);
                        trsdipy[index] += myci->myhf->cmx(i, a / 2) * myci->myhf->cmx(j, r / 2) * myci->myhf->muymx(i, j);
                        trsdipz[index] += myci->myhf->cmx(i, a / 2) * myci->myhf->cmx(j, r / 2) * myci->myhf->muzmx(i, j);
                    }
                }
            }
        }
    }
    /*Gør at integralerne kan tilgås via. indekserne 0, 1 og 2 i trsdip*/
    trsdip[0] = &trsdipx;
    trsdip[1] = &trsdipy;
    trsdip[2] = &trsdipz;
}

//--------------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
// PUBLIC FUNKTIONER
////////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------------

void TDHF::statpol()
{
    valarray<double> statpolarizabilities(9); // En variabel, der skal gemme på de statiske polarisabiliteter
    DynMatrix apbinv = (amx + bmx).inverse(); //(A + B)-invers (Toulouse ligning 44)
    valarray<double> tmp(nosinglex);          // En midlertidig variabel, der skal gemme på vektor-matrixproduktet af overgangsdipolintegralerne og (A + B)-invers (Toulouse ligning 44)
    for (int I = 0; I < 3; I++)
    { // Loop over x, y og z for overgangsdipolintegralerne i x-, y- og z-retningerne
        for (int J = 0; J < 3; J++)
        {
            for (int ar = 0; ar < nosinglex; ar++)
            { // Loop over kolonnerne i (A + B)-invers
                for (int bs = 0; bs < nosinglex; bs++)
                {                                                 // Loop over elementerne i dipolintegralvektorerne og rækkerne i (A + B)-invers
                    tmp[ar] += (*trsdip[I])[bs] * apbinv(bs, ar); // Vektor-matrixmultiplikationen i Toulouse ligning 44 (der ganges med 2 nedenfor)
                }
                statpolarizabilities[3 * I + J] += 2 * tmp[ar] * (*trsdip[J])[ar]; // Prikproduktet i Toulouse ligning 44
                tmp[ar] = 0;                                                       // Reset tmp
            }
        }
    }
    time_t time2 = time(NULL);
    // Output resultaterne
    /********************************************************************/
    outfile.precision(5);
    outfile << "xx ";
    outfile.width(10);
    outfile << statpolarizabilities[0] << '\n';
    outfile << "xy ";
    outfile.width(10);
    outfile << statpolarizabilities[1] << '\n';
    outfile << "xz ";
    outfile.width(10);
    outfile << statpolarizabilities[2] << '\n';
    outfile << "yx ";
    outfile.width(10);
    outfile << statpolarizabilities[3] << '\n';
    outfile << "yy ";
    outfile.width(10);
    outfile << statpolarizabilities[4] << '\n';
    outfile << "yz ";
    outfile.width(10);
    outfile << statpolarizabilities[5] << '\n';
    outfile << "zx ";
    outfile.width(10);
    outfile << statpolarizabilities[6] << '\n';
    outfile << "zy ";
    outfile.width(10);
    outfile << statpolarizabilities[7] << '\n';
    outfile << "zz ";
    outfile.width(10);
    outfile << statpolarizabilities[8] << "\n--\n";
    outfile << "Den isotrope værdi: " << (statpolarizabilities[0] + statpolarizabilities[4] + statpolarizabilities[8]) / 3.0 << "\n\n";
    outfile.precision(0);
    outfile << "**Beregningstid: Programmet var " << difftime(time2, time1) << " sekunder om at lave polarisabilitetsberegningen (inklusiv transformeringen af overgangsdipolintegralerne, men med A-matricen genbrugt fra CIS-beregningen).\n\n";
    /********************************************************************/
}

//--------------------------------------------------------------------------------

void TDHF::dynpol(double frequency)
{
    valarray<double> dynpolarizabilities(9); // Variablen, der skal gemme de dynamiske polarisabiliteter
    DynMatrix tdhessianinv = abbamx;         // Gør klar til den inverse af Hessian-matricen. OBS: Det jeg kalder Hessian-matricen her er den "korrigerede" ABB*A*-matrix (se ligning 85 i Toulouse)
    for (int i = 0; i < nosinglex; i++)
    {                                                            // Loop over første halvdel af diagonalen i Hessian-matricen
        tdhessianinv(i, i) -= frequency;                         // Første halvdel af diagonalen ændres
        tdhessianinv(i + nosinglex, i + nosinglex) += frequency; // Sidste halvdel af diagonalen ændres
    }
    tdhessianinv = tdhessianinv.inverse(); //Invertér

    valarray<double> tmp(2 * nosinglex); // En midlertidig variabel, der skal gemme på vektor-matrix-produktet af (trsdip[I] trsdip[I]*) og Hessianen, som det fremgår i ligning 91 i Toulouse (bemærk at elementerne er reelle)
    for (int I = 0; I < 3; I++)
    { // Loop over x, y og z for overgangsdipolintegralerne i x-, y- og z-retningerne
        for (int J = 0; J < 3; J++)
        {
            for (int ar = 0; ar < nosinglex; ar++)
            { // Loop over halvdelen af kolonnerne i Hessian'en
                for (int bs = 0; bs < nosinglex; bs++)
                {                                                                                           // Loop over elementerne i dipolintegralvektorerne og halvdelen af rækkerne i Hessian'en
                    tmp[ar] += (*trsdip[I])[bs] * tdhessianinv(bs, ar);                                     // Leddene fra øverste venstre halvdel af Hessian'en
                    tmp[ar] += (*trsdip[I])[bs] * tdhessianinv(bs + nosinglex, ar);                         // Leddene fra nederste venstre halvdel af Hessian'en
                    tmp[ar + nosinglex] += (*trsdip[I])[bs] * tdhessianinv(bs, ar + nosinglex);             // Leddene fra øverste højre halvdel af Hessian'en
                    tmp[ar + nosinglex] += (*trsdip[I])[bs] * tdhessianinv(bs + nosinglex, ar + nosinglex); // Leddene fra nederste højre halvdel af Hessian'en
                }
                dynpolarizabilities[3 * I + J] += tmp[ar] * (*trsdip[J])[ar];             // Prikproduktet mellem vektor-matrix-produktet og (trsdip[J] trsdip[J]*) (bemærk igen at elementerne er reelle)
                dynpolarizabilities[3 * I + J] += tmp[ar + nosinglex] * (*trsdip[J])[ar]; // Fordi vi kun summer over halvdelen af kolonnerne i Hessian'en / halvdelen af elementerne i (trsdip[J] trsdip[J]*) har vi to led her
                tmp[ar] = 0;                                                              // Reset
                tmp[ar + nosinglex] = 0;                                                  // Reset
            }
        }
    }
    time_t time2 = time(NULL);
    // Output resultaterne
    /********************************************************************/
    outfile.precision(5);
    outfile << "omega = " << frequency << "\n\n";
    outfile << "xx ";
    outfile.width(10);
    outfile << dynpolarizabilities[0] << '\n';
    outfile << "xy ";
    outfile.width(10);
    outfile << dynpolarizabilities[1] << '\n';
    outfile << "xz ";
    outfile.width(10);
    outfile << dynpolarizabilities[2] << '\n';
    outfile << "yx ";
    outfile.width(10);
    outfile << dynpolarizabilities[3] << '\n';
    outfile << "yy ";
    outfile.width(10);
    outfile << dynpolarizabilities[4] << '\n';
    outfile << "yz ";
    outfile.width(10);
    outfile << dynpolarizabilities[5] << '\n';
    outfile << "zx ";
    outfile.width(10);
    outfile << dynpolarizabilities[6] << '\n';
    outfile << "zy ";
    outfile.width(10);
    outfile << dynpolarizabilities[7] << '\n';
    outfile << "zz ";
    outfile.width(10);
    outfile << dynpolarizabilities[8] << "\n--\n";
    outfile << "Den isotrope værdi: " << (dynpolarizabilities[0] + dynpolarizabilities[4] + dynpolarizabilities[8]) / 3.0 << "\n\n";
    outfile.precision(0);
    outfile << "**Beregningstid: Programmet var " << difftime(time2, time1) << " sekunder om at lave polarisabilitetsberegningen i dynpol (inklusiv transformeringen af overgangsdipolintegralerne, men med A-matricen genbrugt fra CIS-beregningen).\n\n";
    /********************************************************************/
}

//--------------------------------------------------------------------------------

void TDHF::dynpol2(double frequency)
{
    valarray<double> dynpolarizabilities(9); // Variablen, der skal gemme de dynamiske polarisabiliteter
    DynMatrix tdhessian = abbamx;            // Gør klar til Hessian-matricen. OBS: Det jeg kalder Hessian-matricen her er den "korrigerede" ABB*A*-matrix (se ligning 85 i Toulouse)
    Eigen::VectorXd XY(2 * nosinglex);       // Vektoren (X Y) som vi bestemmer ved at løse ligning 85 i Toulouse
    for (int i = 0; i < nosinglex; i++)
    {                                                         // Loop over første halvdel af diagonalen i Hessian-matricen
        tdhessian(i, i) -= frequency;                         // Første halvdel af diagonalen ændres
        tdhessian(i + nosinglex, i + nosinglex) += frequency; // Sidste halvdel af diagonalen ændres
    }

    // Vektorer, der gemmer på overgangsdipolintegralerne - dobbelt: (mu_i mu_i*) (Bemærk at mu_i = mu_i* i min implementering)
    // Vektorerne svarer til vektoren (V W) (ligning 85 i Toulouse)
    Eigen::VectorXd dipxdipx(2 * nosinglex);
    Eigen::VectorXd dipydipy(2 * nosinglex);
    Eigen::VectorXd dipzdipz(2 * nosinglex);
    for (int i = 0; i < nosinglex; i++)
    {
        dipxdipx(i) = trsdipx[i];
        dipxdipx(i + nosinglex) = trsdipx[i];
        dipydipy(i) = trsdipy[i];
        dipydipy(i + nosinglex) = trsdipy[i];
        dipzdipz(i) = trsdipz[i];
        dipzdipz(i + nosinglex) = trsdipz[i];
    }
    Eigen::VectorXd *dipidipi[3] = {&dipxdipx, &dipydipy, &dipzdipz};

    for (int J = 0; J < 3; J++)
    {                                              // Loop over x, y, z
        XY = tdhessian.lu().solve((*dipidipi[J])); // Bestem (X Y)
        for (int I = 0; I < 3; I++)
        {
            // Polarisabiliteten IJ findes som prikproduktet mellem (mu_i mu_i*) og XY
            for (int ar = 0; ar < nosinglex; ar++)
            {
                dynpolarizabilities[3 * I + J] += (*dipidipi[I])(ar)*XY(ar);
                dynpolarizabilities[3 * I + J] += (*dipidipi[I])(ar + nosinglex) * XY(ar + nosinglex);
            }
        }
    }
    double isotpol = (dynpolarizabilities[0] + dynpolarizabilities[4] + dynpolarizabilities[8]) / 3.0;
    time_t time2 = time(NULL);
    // Output resultaterne
    /********************************************************************/
    outfile.precision(5);
    outfile << "omega = " << frequency << "\n\n";
    outfile << "xx ";
    outfile.width(10);
    outfile << dynpolarizabilities[0] << '\n';
    outfile << "xy ";
    outfile.width(10);
    outfile << dynpolarizabilities[1] << '\n';
    outfile << "xz ";
    outfile.width(10);
    outfile << dynpolarizabilities[2] << '\n';
    outfile << "yx ";
    outfile.width(10);
    outfile << dynpolarizabilities[3] << '\n';
    outfile << "yy ";
    outfile.width(10);
    outfile << dynpolarizabilities[4] << '\n';
    outfile << "yz ";
    outfile.width(10);
    outfile << dynpolarizabilities[5] << '\n';
    outfile << "zx ";
    outfile.width(10);
    outfile << dynpolarizabilities[6] << '\n';
    outfile << "zy ";
    outfile.width(10);
    outfile << dynpolarizabilities[7] << '\n';
    outfile << "zz ";
    outfile.width(10);
    outfile << dynpolarizabilities[8] << "\n--\n";
    outfile << "Den isotrope værdi: " << isotpol << "\n\n";
    outfile.precision(0);
    outfile << "**Beregningstid: Programmet var " << difftime(time2, time1) << " sekunder om at lave polarisabilitetsberegningen i dynpol2 (inklusiv transformeringen af overgangsdipolintegralerne, men med A-matricen genbrugt fra CIS-beregningen).\n\n";
    /********************************************************************/
}

//--------------------------------------------------------------------------------

void TDHF::excitations()
{
    DynMatrix abmbmamx = abbamx; // Gør klar til AB(-B)(-A)-matricen (Crawford, projekt 12)
    for (int i = nosinglex; i < nosinglex * 2; i++)
    {
        for (int j = 0; j < nosinglex * 2; j++)
        {
            abmbmamx(i, j) = -abmbmamx(i, j); // Lav ændringen
        }
    }
    Eigen::EigenSolver<DynMatrix> tdhfsolved(abmbmamx); // Diagonalisér matricen
    double *tdhfexcienergies = new double[nosinglex * 2];
    for (int i = 0; i < nosinglex * 2; i++)
    { // Sum over antallet af egenværdier
        tdhfexcienergies[i] = tdhfsolved.eigenvalues()(i).real();
    }
    sort(tdhfexcienergies, (tdhfexcienergies + nosinglex * 2));

    // Output resultaterne
    /********************************************************************/
    outfile << "Eksitationsenergier:\n\n";
    outfile.width(20);
    outfile << "Hartrees";
    outfile.width(18);
    outfile << "eV\n";
    for (int i = 0; i < 2 * nosinglex; i++)
    {
        outfile.setf(ios_base::left, ios_base::adjustfield);
        outfile.width(3);
        outfile << i + 1;
        outfile.setf(ios_base::right, ios_base::adjustfield);
        outfile.precision(12);
        outfile.width(20);
        outfile << tdhfexcienergies[i];
        outfile.width(20);
        outfile << tdhfexcienergies[i] * 27.2114 << '\n'; // Konvertering til eV
    }
    outfile << '\n';
    /********************************************************************/
    delete[] tdhfexcienergies;
}

//--------------------------------------------------------------------------------

void TDHF::excitationsquicker()
{
    DynMatrix targetmx = (amx + bmx) * (amx - bmx);   // Lav (A+B)*(A-B)-matricen (Crawford, projekt 12)
    Eigen::EigenSolver<DynMatrix> solvedmx(targetmx); // Diagonalisér matricen
    double *tdhfquickexcienergies = new double[nosinglex];
    for (int i = 0; i < nosinglex; i++)
    {                                                                      // Sum over antallet af egenværdier
        tdhfquickexcienergies[i] = sqrt(solvedmx.eigenvalues()(i).real()); // Husk kvadratroden
    }
    sort(tdhfquickexcienergies, (tdhfquickexcienergies + nosinglex));
    // Output resultaterne
    /********************************************************************/
    outfile << "Eksitationsenergier:\n\n";
    outfile.width(20);
    outfile << "Hartrees";
    outfile.width(18);
    outfile << "eV\n";
    for (int i = 0; i < nosinglex; i++)
    {
        outfile.setf(ios_base::left, ios_base::adjustfield);
        outfile.width(3);
        outfile << i + 1;
        outfile.setf(ios_base::right, ios_base::adjustfield);
        outfile.precision(12);
        outfile.width(20);
        outfile << tdhfquickexcienergies[i];
        outfile.width(20);
        outfile << tdhfquickexcienergies[i] * 27.2114 << '\n'; // Konvertering til eV
    }
    outfile << '\n';
    /********************************************************************/
    delete[] tdhfquickexcienergies;
}

//--------------------------------------------------------------------------------