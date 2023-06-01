#include <iostream>
#include <fstream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <time.h>
#include "../include/HF.h"
#include "../include/DIIS.h"
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"
using namespace std;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> DynMatrix;

////////////////////////////////////////////////////////////////////////////////
//CONSTRUCTOR
////////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------------

HF::HF(Molecule& molecule, const char* const dir) : mymolecule(molecule), integraldir(dir), energyconverged(false)
{
    loadbasissetinfo(); //Load basissætnavnet og antallet af basisfunktioner
    /*Load integralerne*/
    loadintegrals("/s.dat", &smx);
    loadintegrals("/v.dat", &vmx);
    loadintegrals("/t.dat", &tmx);
    loaderiintegrals();
    sethcore(); //Dan hcore
    computeorthsmx(); //Udregn orthogonaliseringsmatricen S^(-1/2) fra overlapmatricen
    moenergy = new double [nobasisfcn]; //Gør klar til vektoren med MO-energierne
}

//--------------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
//DESTRUCTOR
////////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------------

HF::~HF()
{
    /*Slet eritensor*/
    for (int i = 0; i < nobasisfcn; i++) {
        for (int j = 0; j < nobasisfcn; j++) {
            for (int k = 0; k < nobasisfcn; k++) {
                delete[] eritensor[i][j][k]; //Slet alle kolonner i eritensor[i][j] 2D-matricen
            }
            delete[] eritensor[i][j]; //Slet eritensor[i][j] 2D-matricen
        }
        delete[] eritensor[i]; //Slet alle kolonner i eritensor
    }
    delete[] eritensor;
    delete[] moenergy;
}

//--------------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
//PRIVATE LOAD-FUNKTIONER
////////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------------

void HF::loadbasissetinfo()
{
    char basisfile[200]; //Filstien til filen med navnet på det anvendte basissæt
    snprintf(basisfile, sizeof(basisfile), "%s%s", integraldir, "/basis.txt"); //Print filstien til variablen

    char noaofile[200]; //Filstien til filen med antallet af brugte basisfunktioner
    snprintf(noaofile, sizeof(noaofile), "%s%s", integraldir, "/noao.txt"); // Print filstien til variablen

    ifstream infile1(basisfile, ios_base::in);
    ifstream infile2(noaofile, ios_base::in);

    if (!infile1) {
        cout << "Fejl: Programmet kunne ikke finde filen \"" << basisfile << "\".\n";
        cout << "Bemærk at programmet forventer filnavnet \"basis.txt\" for filen med navnet på basissættet.\n";
        exit(1);
    }
    else if (!infile2) {
        cout << "Fejl: Programmet kunne ikke finde filen \"" << noaofile << "\".\n";
        cout << "Bemærk at programmet forventer filnavnet \"noao.txt\" for filen med antallet af basisfunktioner.\n";
        exit(1);
    }

    infile1.getline(basissetname, sizeof(basissetname), '\n'); //Henter basissætnavnet

    char noao[10];
    infile2.getline(noao, sizeof(noao), '\n'); //Henter antallet af basisfunktioner
    nobasisfcn = atoi(noao);
    if (nobasisfcn <= 0) {
        cout << "Fejl: Antallet af basisfunktioner er indlæst til " << noao << ".\n";
        exit(1);
    }

    infile1.close();
    infile2.close();
}

//--------------------------------------------------------------------------------

void HF::loadintegrals(const char* const fname, DynMatrix* mx)
{
    char integralfile[200]; //Filstien til filen med integraldata
    snprintf(integralfile, sizeof(integralfile), "%s%s", integraldir, fname); //Print filstien til variablen


    ifstream infile(integralfile, ios_base::in); //Indlæs filen
    if (!infile) {
        cout << "Fejl: Programmet kunne ikke finde filen \"" << integralfile << "\" i HF::loadintegrals.\n";
        exit(1);
    }

    
    char line[30]; //Den linje / det tekststykke, der læses i .dat-filen
    mx->resize(nobasisfcn, nobasisfcn); //Giv matricen den rette størrelse


    /*OBS: Indlæsningen af én-elektron-integralerne er designet efter formateringen af datafilerne fra projektet*/
    for (int i = 0; i < nobasisfcn; i++) { //Loop over alle rækkerne i matricen
        for (int j = 0; j <= i; j++) { //Loop over kolonnerne, men kun ind til (og med) diagonalen
            /*Ignorer mellemrum og informationen om indekserne*/
            while (infile.peek() == ' ') {infile.ignore();} //Ignorér mellemrum
            infile.getline(line, sizeof(line), ' '); //Ignorér første indeks
            while (infile.peek() == ' ') {infile.ignore();} //Ignorér mellemrum igen
            infile.getline(line, sizeof(line), ' '); //Ignorér andet indeks
            while (infile.peek() == ' ') {infile.ignore();}


            infile.getline(line, sizeof(line), '\n'); //Hent integralet
            (*mx)(i, j) = atof(line); //Læs værdien til matricen
            (*mx)(j, i) = (*mx)(i, j); //Pga. hermiticitet
        }
    }


    //Tjek for eof (ignorér tom linje)
    if (!infile.eof() && infile.peek() != -1) {
        cout << "Advarsel: Der er uoverensstemmelse mellem antallet af basisfunktioner og den forventede længde af filen \"" << integralfile << "\".\n";
    }


    infile.close();
}

//--------------------------------------------------------------------------------

void HF::loaderiintegrals()
{
    char eriintegralfile[200]; //Filstien til filen med integraldata for to-elektron-integralerne
    snprintf(eriintegralfile, sizeof(eriintegralfile), "%s%s", integraldir, "/eri.dat"); //Printer filstien til variablen
    

    ifstream infile(eriintegralfile, ios_base::in);
    if (!infile) {
        cout << "Fejl: Programmet kunne ikke finde filen \"" << eriintegralfile << "\".\n";
        cout << "Bemærk at programmet forventer filnavnet \"eri.dat\" for filen med integraldata for J- og K-integralerne.\n";
        exit(1);
    }


    char line[30]; //Den linje / det tekststykke, der læses i eri.dat-filen
    int index[4]; //Indekserne i eri.dat-filen
    eritensor = new double*** [nobasisfcn]; //Lav rækkerne i eritensor
    for (int i = 0; i < nobasisfcn; i++) {
        eritensor[i] = new double** [nobasisfcn]; //Lav kolonnerne i eritensor (hvert element er en 2D-matrix)
        for (int j = 0; j < nobasisfcn; j++) {
            eritensor[i][j] = new double* [nobasisfcn]; //Lav rækkerne i eritensor[i][j] 2D-matricen
            for (int k = 0; k < nobasisfcn; k++) {
                eritensor[i][j][k] = new double [nobasisfcn]; //Lav kolonnerne i eritensor[i][j] 2D-matricen
                for (int l = 0; l < nobasisfcn; l++) {
                    eritensor[i][j][k][l] = 0; //Initiér alle elementerne til 0
                }
            }
        }
    }


    /*OBS: Denne rækkefølge af loops sikrer, at alle unikke integraler indlæses.
    Indlæsningen af integralerne er designet efter, at symmetriidentiske integraler
    kun optræder en gang i .dat-filen*/
    for (int i = 0; i < nobasisfcn; i++) { //Loop over alle rækker i eritensor
        for (int j = 0; j <= i; j++) { //Loop over kolonnerne i eritensor, men kun ind til (og med) diagonalen
            for (int k = 0; k < nobasisfcn; k++) { //Loop over alle rækker eritensor[i][j] 2D-matricen
                for (int l = 0; l <= k; l++) { //Loop over kolonner i eritensor[i][j] 2D-matricen, men kun ind til (og med) diagonalen
                    //Hent indekserne
                    for (int m = 0; m < 4; m++) {
                        while (infile.peek() == ' ') infile.ignore(); //Ignorér mellemrum
                        infile.getline(line, sizeof(line), ' '); //Hent indekset
                        if (! (index[m] = atoi(line))) {
                            cout << "Fejl: Programmet kunne ikke indlæse indekset \"" << line << "\" i filen \"" << eriintegralfile << "\".\n";
                            exit(1);
                        }
                        else if (index[m] > nobasisfcn) {
                            cout << "Fejl: Det indlæste indeks \"" << line << "\" i filen \"" << eriintegralfile << "\" er større end antallet af basisfunktioner.";
                            exit(1);
                        }
                        index[m] -= 1; //VIGTIG: Integral 1234 bliver til 0123 i repræsentationen i programmet (tilsvarende for alle andre integraler)
                    }


                    //Hent integralet
                    while (infile.peek() == ' ') infile.ignore(); //Ignorér mellemrum mellem sidste indeks og integralet
                    infile.getline(line, sizeof(line), '\n'); //Hent integralet
                    double eriintegral = atof(line); //Integralet gemmes i en midlertidig variabel


                    //Gem integralet i eritensor. Symmetrien af integralerne bruges.
                    eritensor[index[0]][index[1]][index[2]][index[3]] = eriintegral;
                    eritensor[index[1]][index[0]][index[2]][index[3]] = eriintegral;
                    eritensor[index[0]][index[1]][index[3]][index[2]] = eriintegral;
                    eritensor[index[1]][index[0]][index[3]][index[2]] = eriintegral;
                    eritensor[index[2]][index[3]][index[0]][index[1]] = eriintegral;
                    eritensor[index[3]][index[2]][index[0]][index[1]] = eriintegral;
                    eritensor[index[2]][index[3]][index[1]][index[0]] = eriintegral;
                    eritensor[index[3]][index[2]][index[1]][index[0]] = eriintegral;
                    

                    //Afslut hvis eof eller en tom linje mødes
                    if (infile.eof() || infile.peek() == -1) {
                        infile.close();
                        return;
                    }
                }
            }
        }
    }


    infile.close();
}

//--------------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
//PRIVATE COMPUTE-FUNKTIONER
////////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------------

void HF::computeorthsmx()
{
    Eigen::SelfAdjointEigenSolver<DynMatrix> solvedsmx(smx); //Diagonalisér overlapmatricen (find egenvektorer og egenværdier)
    DynMatrix eigvalmxinvsqrt = DynMatrix::Zero(nobasisfcn, nobasisfcn); //Gør klar til den diagonale matrix med egenværdierne opløftet til -1/2
    for (int i = 0; i < nobasisfcn; i++) {
        double eigenvalue = solvedsmx.eigenvalues()(i); //Gem egenværdien i en midlertidig variabel
        if (eigenvalue != 0 && eigenvalue > 0) {
            eigvalmxinvsqrt(i, i) = 1.0 / sqrt(eigenvalue); //Sæt de diagonale elementer lig med egenværdierne opløftet til -1/2
        }
        else {
            cout << "Fejl: En af egenværdierne til overlapmatricen er indlæst til 0 eller et negativt tal.\n";
            exit(1);
        }
    }
    DynMatrix eigvecmx = solvedsmx.eigenvectors(); //Gem egenvektorerne til overlapmatricen i eigvecmx
    orthsmx = eigvecmx * eigvalmxinvsqrt * eigvecmx.transpose(); //Dan orthogonaliseringsmatricen
}

//--------------------------------------------------------------------------------

void HF::diagonalizefmx(DynMatrix* mx)
{

    DynMatrix fmxprime = orthsmx * (*mx) * orthsmx; //Transformér Fock-matricen (mx) til en orthonormal basis
    Eigen::SelfAdjointEigenSolver<DynMatrix> solvedfmxprime(fmxprime); //Diagonalisér den transformerede Fock-matrix (find egenvektorer og egenværdier)
    cmx = orthsmx * solvedfmxprime.eigenvectors(); //Transformér de orthonormale egenvektorer til den rigtige basis
    for (int i = 0; i < nobasisfcn; i++) moenergy[i] = solvedfmxprime.eigenvalues()(i); //Gem MO-energierne (egenværdierne)
    sort(moenergy, &(moenergy[nobasisfcn])); //Sortér MO-energierne i voksende orden (i tilfælde af, at de ikke er sorterede i forvejen)
}

//--------------------------------------------------------------------------------

void HF::computedmx()
{
    dmx = DynMatrix::Zero(nobasisfcn, nobasisfcn); //Initiér alle indekser til 0
    for (int i = 0; i < nobasisfcn; i++) { //En sum over alle rækker i densitetsmatricen
        for (int j = 0; j <= i; j++) { //En sum over alle kolonner i densitetsmatricen men kun ind til diagonalen
            for (int k = 0; k < mymolecule.getnooccmos(); k++) { //En sum over alle okkuperede molekylorbitaler
                dmx(i, j) += cmx(i, k) * cmx(j, k);
            }
            dmx(j, i) = dmx(i, j); //Pga. symmetrien D(ij) = D(ji)
        }
    }
}

//--------------------------------------------------------------------------------

void HF::computeenergy(DynMatrix* mx)
{
    elecenergy = 0; //Initiér den elektroniske energi til 0
    for (int i = 0; i < nobasisfcn; i++) {  //En sum over rækkerne i Fock-matricen
        for (int j = 0; j < nobasisfcn; j++) { //En sum over kolonnerne i Fock-matricen
            elecenergy += dmx(i, j) * (hcore(i, j) + (*mx)(i, j));
        }
    }
}

//--------------------------------------------------------------------------------

void HF::computefmx()
{
    /*Kontrahér to-elektron-integralerne*/
    jmx = DynMatrix::Zero(nobasisfcn, nobasisfcn); //Initiér elementerne i Coulomb-matricen til 0
    kmx = DynMatrix::Zero(nobasisfcn, nobasisfcn); //Initiér elementerne i exchange-matricen til 0
    for (int i = 0; i < nobasisfcn; i++) { //En sum over alle rækker i Coulomb- og exchange-matricerne
        for (int j = 0; j < nobasisfcn; j++) { //En sum over alle kolonner i Coulomb- og exchange-matricerne
            for (int k = 0; k < nobasisfcn; k++) { //En sum over alle rækker i densitetsmatricen og rækkerne i 2D-matricen eritensor[i][j] 
                for (int l = 0; l < nobasisfcn; l++) { //En sum over alle kolonner i densitetsmatricen og kolonnerne i 2D-matricen eritensor[i][j]
                    jmx(i, j) += dmx(k, l) * eritensor[i][j][k][l];
                    kmx(i, j) += dmx(k, l) * eritensor[i][k][j][l];
                }
            }
        }
    }
    fmx = hcore + (2*jmx) - kmx;
}

//--------------------------------------------------------------------------------

double HF::computermsd(DynMatrix* mx1, DynMatrix* mx2)
{
    double rmsd = 0;
    for (int i = 0; i < nobasisfcn; i++) {
        for (int j = 0; j < nobasisfcn; j++) {
            rmsd += ( ((*mx1)(i, j) - (*mx2)(i, j)) * ((*mx1)(i, j) - (*mx2)(i, j)));
        }
    }
    return sqrt(rmsd);
}

//--------------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
//PUBLIC PRINTFUNKTIONER
////////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------------

void HF::printmx(ofstream& outfile, DynMatrix* mx)
{
    int w = 18; //Bredden
    int p = 12; //Præcisionen
    outfile.setf(ios_base::fixed, ios_base::floatfield);
    for (int i = 0; i < nobasisfcn; i++) {
        for (int j = 0; j < nobasisfcn; j++) {
            outfile.width(w);
            outfile.precision(p);
            outfile << (*mx)(i, j);
            outfile << ' ';
        }
        outfile << '\n';
    }
}

//--------------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
//PUBLIC LOAD-FUNKTIONER
////////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------------

void HF::loaddipoleintegrals()
{
    loadintegrals("/mux.dat", &muxmx);
    loadintegrals("/muy.dat", &muymx);
    loadintegrals("/muz.dat", &muzmx);
}

//--------------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
//PUBLIC COMPUTE-FUNKTIONER
////////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------------

void HF::HFSCF(int maxiter, double delta1, double delta2, ofstream& outfile)
{
    time_t time1 = time(NULL);
    //Start med at bestemme den 0. elektroniske energi ud fra startgættet
    /********************************************************************/
    setfmxstartguess(&hcore); //Angiv startgættet for Fock-matricen
    diagonalizefmx(&fmx); //Diagonalisér Fock-matricen (hcore)
    computedmx(); //Udregn densitetsmatricen
    computeenergy(&fmx); //Udregn den 0. elektroniske energi
    /********************************************************************/

    //Formatering
    /********************************************************************/
    outfile << "Uden DIIS\n";
    outfile << "------------------------------------------------------------------------------------------------------------------------------------\n";
    outfile << "Iter.";
    outfile.width(28);
    outfile << "Elek. energi / Hartrees";
    outfile.width(23);
    outfile << "DeltaE / Hartrees";
    outfile.width(21);
    outfile << "RMSD\n";
    outfile.width(3);
    outfile << "0";
    outfile.width(26);
    outfile << elecenergy << '\n';
    /********************************************************************/

    //HF-SCF
    /********************************************************************/
    int counter = 1; //Antallet af SCF-iterationer
    double rmsd = 0; //root-mean-square-forskellen mellem elementerne i to på hinanden følgende densitetsmatricer
    double deltae = 0; //Forskelllen mellem to på hinanden følgende elektroniske energier
    double tempenergy = elecenergy; //En variabel, der gemmer på den elektroniske energi
    DynMatrix tempdmx = dmx; //En variabel, der gemmer på densitetsmatricen
    do {
        computefmx(); //Udregn den nye Fock-matrix
        diagonalizefmx(&fmx); //Diagonalisér Fock-matricen
        computedmx(); //Udregn den nye densitetsmatrix
        computeenergy(&fmx); //Udregn den nye elektroniske energy
        deltae = elecenergy - tempenergy; //Udregn energiforskellen
        rmsd = computermsd(&dmx, &tempdmx); //Udregn rmsd
        tempenergy = elecenergy;
        tempdmx = dmx;

        //Output
        /********************************************************************/
        outfile.width(3);
        outfile << counter;
        outfile.width(26);
        outfile << elecenergy;
        outfile.width(26);
        outfile << deltae;
        outfile.width(26);
        outfile << rmsd << '\n';
        /********************************************************************/
        counter++;
    } while (rmsd > delta1 && fabs(deltae) > delta2 && counter <= maxiter);
    /********************************************************************/
    time_t time2 = time(NULL);

    //Output de endelige resultater
    /********************************************************************/
    outfile << "------------------------------------------------------------------------------------------------------------------------------------\n\n";
    if (rmsd < delta1) {
        outfile << "Succes! Energien er konvergeret. RMSD er under den angivne grænse: " << delta1 << "\n";
        energyconverged = true;
    }
    else if (fabs(deltae) < delta2) {
        outfile << "Succes! Energien er konvergeret. Energiforskellen er under den angivne grænse: " << delta2 << "\n";
        energyconverged = true;
    }
    else {
        outfile << "Det maksimale antal iterationer er nået. Energiforskellen eller RMSD falder ikke under de angivne grænser.\n";
        energyconverged = false;
    }
    outfile.precision(0);
    outfile << "**Beregningstid: Programmet var " << difftime(time2, time1) << " sekunder om at udføre den ikke-accelererede HF SCF-beregning.\n\n";
    outfile.precision(12);
    outfile << "Den endelige elektroniske energi / Hartrees:\n";
    outfile << elecenergy << "\n\n";
    outfile << "Kernefrastødningen / Hartrees:\n";
    outfile << mymolecule.getncrepenergy() << "\n\n";
    outfile << "Den endelige HF-energi / Hartrees:\n";
    outfile << elecenergy + mymolecule.getncrepenergy() << "\n\n";
    /********************************************************************/
}

//--------------------------------------------------------------------------------

void HF::HFSCFDIIS(int maxiter, double delta1, double delta2, ofstream& outfile, int nosavedmcs)
{
    time_t time1 = time(NULL);
    DIIS mydiis(this, nosavedmcs); //Gør klar til DIIS. nosavedmcs er antallet af fejl- og Fock-matricer, der skal gemmes

    //Start med at bestemme den 0. elektroniske energi ud fra startgættet
    /*******************************************************************/
    setfmxstartguess(&hcore); //Angiv startgættet for Fock-matricen
    diagonalizefmx(&fmx); //Diagonalisér Fock-matricen (hcore) --> find koefficienterne til basisfunktionerne
    computedmx(); //Udregn densitetsmatricen
    computeenergy(&fmx); //Udregn den 0. elektroniske energi
    /*******************************************************************/

    //Formatering
    /*******************************************************************/
    outfile << "MED DIIS. Max " << nosavedmcs << " gemte fejlvektorer.\n";
    outfile << "------------------------------------------------------------------------------------------------------------------------------------\n";
    outfile << "Iter.";
    outfile.width(28);
    outfile << "Elek. energi / Hartrees";
    outfile.width(23);
    outfile << "DeltaE / Hartrees";
    outfile.width(21);
    outfile << "RMSD\n";
    outfile.width(3);
    outfile << "0";
    outfile.width(26);
    outfile << elecenergy << '\n';
    /*******************************************************************/
    
    //HF-SCF
    /*******************************************************************/
    int counter = 1; //Antallet af SCF-iterationer
    double rmsd = 0; //Root mean square for forskellen mellem elementerne i to på hinanden følgende densitetsmatricer
    double deltae = 0; //Forskelllen mellem to på hinanden følgende elektroniske energier
    double tempenergy = elecenergy; //En variabel, der gemmer på den elektroniske energi
    DynMatrix tempdmx = dmx; //En variabel, der gemmer på densitetsmatricen
    do {
        //Vi er klar til at lave DIIS, når counter > 1
        if (counter > 1) {
            computefmx(); //Udregn den nye Fock-matrix ud fra densitetsmatricen
            mydiis.savefmx(&fmx, (counter - 1) % mydiis.maxdimensions); //Gem den nye Fock-matrix i DIIS-objektet
            mydiis.computeemx((counter - 1) % mydiis.maxdimensions); //Lav en ny fejlmatrix ud fra den netop gemte Fock-matrix
            mydiis.computebmx(); //Udregn B-matricen i DIIS
            mydiis.solvebmx(); //Find koefficienterne til den lineære kombination af Fock-matricer
            mydiis.computelcfmx(); //Bestem Fock-matrix-gættet
            diagonalizefmx(&mydiis.lcfmx); //Find LCAO-MO-koefficienterne ud fra Fock-matrix-gættet
            computedmx(); //Udregn den nye densitetsmatrix
            computeenergy(&mydiis.lcfmx); //Udregn den nye elektroniske energy fra densitetsmatricen og Fock-matrix-gættet
            deltae = elecenergy - tempenergy; //Udregn energiforskellen
            rmsd = computermsd(&dmx, &tempdmx); //Udreng rmsd
            tempenergy = elecenergy;
            tempdmx = dmx;

            //Output
            /*******************************************************************/
            outfile.width(3);
            outfile << counter;
            outfile.width(26);
            outfile << elecenergy;
            outfile.width(26);
            outfile << deltae;
            outfile.width(26);
            outfile << rmsd << '\n';
            /*******************************************************************/
            counter++;
        }
        
        //Hvis counter == 1 kan vi ikke lave DIIS endnu
        else {
            computefmx(); //Udregn den nye Fock-matrix
            mydiis.savefmx(&fmx, (counter - 1) % mydiis.maxdimensions); //Gem Fock-matricen i DIIS-objektet
            mydiis.computeemx((counter - 1) % mydiis.maxdimensions); //Lav en fejlmatrix ud fra den netop gemte Fock-matrix
            diagonalizefmx(&fmx); //Diagonalisér Fock-matricen --> find koefficienterne til basisfunktionerne
            computedmx(); //Udregn den nye densitetsmatrix
            computeenergy(&fmx); //Udregn den nye elektroniske energy
            deltae = elecenergy - tempenergy; //Udregn energiforskellen
            rmsd = computermsd(&dmx, &tempdmx); //Udregn rmsd
            tempenergy = elecenergy;
            tempdmx = dmx;

            //Output
            /*******************************************************************/
            outfile.width(3);
            outfile << counter;
            outfile.width(26);
            outfile << elecenergy;
            outfile.width(26);
            outfile << deltae;
            outfile.width(26);
            outfile << rmsd << '\n';
            /*******************************************************************/
            counter++;
        }
    } while (rmsd > delta1 && fabs(deltae) > delta2 && counter <= maxiter);
    /*******************************************************************/
    time_t time2 = time(NULL);

    //Output de endelige resultater
    /*******************************************************************/
    outfile << "------------------------------------------------------------------------------------------------------------------------------------\n\n";
    if (rmsd < delta1) {
        outfile << "Succes! Energien er konvergeret. RMSD er under den angivne grænse: " << delta1 << "\n";
        energyconverged = true;
    }
    else if (fabs(deltae) < delta2) {
        outfile << "Succes! Energien er konvergeret. Energiforskellen er under den angivne grænse: " << delta2 << "\n";
        energyconverged = true;
    }
    else {
        outfile << "Det maksimale antal iterationer er nået. Energiforskellen eller RMSD falder ikke under de angivne grænser.\n";
        energyconverged = false;
    }
    outfile.precision(0);
    outfile << "**Beregningstid: Programmet var " << difftime(time2, time1) << " sekunder om at udføre den DIIS-accelererede HF SCF-beregning.\n\n";
    outfile.precision(12);
    outfile << "Den endelige elektroniske energi / Hartrees:\n";
    outfile << elecenergy << "\n\n";
    outfile << "Kernefrastødningen / Hartrees:\n";
    outfile << mymolecule.getncrepenergy() << "\n\n";
    outfile << "Den endelige HF-energi / Hartrees:\n";
    outfile << elecenergy + mymolecule.getncrepenergy() << "\n\n";
    /*******************************************************************/
}

////////////////////////////////////////////////////////////////////////////////
//VENNER
////////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------------

void MP2(HF& myhf, ofstream& outfile)
{
    if (!myhf.energyconverged) {
        cout << "Fejl: MP2-energiberegningen kræver, at HF-SCF-energien er konvergeret.\n";
        exit(1);
    }
    time_t t1 = time(NULL);

    //Gør variablen, der skal gemme to-elektron-integralerne i MO-basis, klar
    /********************************************************************/
    double**** moeritensor = new double*** [myhf.nobasisfcn]; //Variablen, der skal gemme to-elektron-integralerne i MO-basis
    for (int i = 0; i < myhf.nobasisfcn; i++) { //Loop over rækker i den ydre 2D-matrix
        moeritensor[i] = new double** [myhf.nobasisfcn]; //Lav rækkerne i den ydre matrix
        for (int j = 0; j < myhf.nobasisfcn; j++) { //Loop over kolonner i den ydre 2D-matrix
            moeritensor[i][j] = new double* [myhf.nobasisfcn]; //Gør klar til rækkerne i de indre matricer
            for (int k = 0; k < myhf.nobasisfcn; k++) { //Loop over rækkerne i de indre 2D-matricer
                moeritensor[i][j][k] = new double [myhf.nobasisfcn]; //Lav rækkerne i de indre matricer
                for (int l = 0; l < myhf.nobasisfcn; l++) {
                    moeritensor[i][j][k][l] = 0; //Initiér alle indekser til 0
                }
            }
        }
    }
    /********************************************************************/

    //Transformer to-elektron-integralerne i AO-basis (eritensor) til MO-basis (moeritensor)
    /********************************************************************/
    //Sum over MO'er (okk. + vir.)
    for (int a = 0; a < myhf.nobasisfcn; a++) {
        for (int b = 0; b < myhf.nobasisfcn; b++) {
            for (int c = 0; c < myhf.nobasisfcn; c++) {
                for (int d = 0; d < myhf.nobasisfcn; d++) {
                    //Sum over AO'er
                    for (int i = 0; i < myhf.nobasisfcn; i++) {
                        for (int j = 0; j < myhf.nobasisfcn; j++) {
                            for (int k = 0; k < myhf.nobasisfcn; k++) {
                                for (int l = 0; l < myhf.nobasisfcn; l++) {
                                    moeritensor[a][b][c][d] += myhf.cmx(i,a)*myhf.cmx(j,b)*myhf.cmx(k,c)*myhf.cmx(l,d) * myhf.eritensor[i][j][k][l];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    /********************************************************************/

    //Udregn andenordenskorrektionen til HF-energien, EMP2
    /********************************************************************/
    double EMP2 = 0;
    for (int a = 0; a < myhf.mymolecule.getnooccmos(); a++) { //Sum over okkuperede molekylorbitaler (rumorbitaler)
        for (int b = 0; b < myhf.mymolecule.getnooccmos(); b++) { //Sum over okkuperede molekylorbitaler (rumorbitaler)
            for (int r = myhf.mymolecule.getnooccmos(); r < myhf.nobasisfcn; r++) { //Sum over virtuelle molekylorbitaler (rumorbitaler)
                for (int s = myhf.mymolecule.getnooccmos(); s < myhf.nobasisfcn; s++) { //Sum over virtuelle molekylorbitaler (rumorbitaler)
                    EMP2 += (moeritensor[a][r][b][s]*(2*(moeritensor[a][r][b][s]) - moeritensor[a][s][b][r]))
                    / (myhf.moenergy[a] + myhf.moenergy[b] - myhf.moenergy[r] - myhf.moenergy[s]);
                }
            }
        }
    }
    /********************************************************************/

    //Slet moeritensor
    /********************************************************************/
    for (int i = 0; i < myhf.nobasisfcn; i++) { //Loop over rækker i den ydre 2D-matrix
        for (int j = 0; j < myhf.nobasisfcn; j++) { //Loop over kolonner i den ydre 2D-matrix
            for (int k = 0; k < myhf.nobasisfcn; k++) { //Loop over rækker i de indre 2D-matricer
                delete[] moeritensor[i][j][k]; //Slet rækker i de indre matricer
            }
            delete[] moeritensor[i][j]; //Slet de indre matricer
        }
        delete[] moeritensor[i]; //Slet rækker i den ydre matrix
    }
    delete[] moeritensor;
    /********************************************************************/
    
    time_t t2 = time(NULL);

    //Output resultatet
    /********************************************************************/
    outfile.precision(12);
    outfile << "E(MP2) / Hartrees:\n" << EMP2 << "\n\n";
    outfile << "Den korrigerede elektroniske energi / Hartrees:\n" << myhf.elecenergy + EMP2 << "\n\n";
    outfile << "Den korrigerede HF-energi / Hartrees:\n" << myhf.elecenergy + myhf.mymolecule.getncrepenergy() + EMP2 << "\n\n";

    outfile.precision(0);
    outfile << "**Beregningstid: Programmet var " << difftime(t2, t1) << " sekunder om at beregne MP2-energien med den langsomme algoritme.\n\n";
    /********************************************************************/
}

//--------------------------------------------------------------------------------

void MP2Quick(HF& myhf, ofstream& outfile)
{
    if (!myhf.energyconverged) {
        cout << "Fejl: MP2-energiberegningen kræver, at HF-SCF-energien er konvergeret.\n";
        exit(1);
    }
    time_t t1 = time(NULL);

    //Gør variablerne, der skal gemme to-elektron-integralerne i MO-basis / delvis MO-basis, klar
    /********************************************************************/
    double**** moeritensor = new double*** [myhf.nobasisfcn]; //Variablen, der skal gemme to-elektron-integralerne i MO-basis
    double**** moeritensortmp = new double*** [myhf.nobasisfcn]; //En variabel, der skal gemme på midlertidige to-elektron-integraler under AO->MO-transformationen
    for (int i = 0; i < myhf.nobasisfcn; i++) { //Loop over rækker i den ydre 2D-matrix
        moeritensor[i] = new double** [myhf.nobasisfcn];
        moeritensortmp[i] = new double** [myhf.nobasisfcn];
        for (int j = 0; j < myhf.nobasisfcn; j++) { //Loop over kolonner i den ydre 2D-matrix
            moeritensor[i][j] = new double* [myhf.nobasisfcn];
            moeritensortmp[i][j] = new double* [myhf.nobasisfcn];
            for (int k = 0; k < myhf.nobasisfcn; k++) { //Loop over rækker i de indre 2D-matricer
                moeritensor[i][j][k] = new double [myhf.nobasisfcn];
                moeritensortmp[i][j][k] = new double [myhf.nobasisfcn];
                for (int l = 0; l < myhf.nobasisfcn; l++) {
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
    for (int a = 0; a < myhf.nobasisfcn; a++) { //Sum over MO'er (rumorbitaler) (okk. + vir.)
        for (int i = 0; i < myhf.nobasisfcn; i++) { //Sum over basisfunktioner i (ij|kl)
            for (int j = 0; j < myhf.nobasisfcn; j++) { //Sum over basisfunktioner i (ij|kl)
                for (int k = 0; k < myhf.nobasisfcn; k++) { //Sum over basisfunktioner i (ij|kl)
                    for (int l = 0; l < myhf.nobasisfcn; l++) { //Sum over basisfunktioner i (ij|kl)
                        moeritensortmp[a][j][k][l] += myhf.cmx(i, a) * myhf.eritensor[i][j][k][l];
                    }
                }
            }
        }
    }
    //Transformer andet indeks (aj|kl) -> (ab|kl)
    for (int b = 0; b < myhf.nobasisfcn; b++) { //Sum over MO'er (rumorbitaler) (okk. + vir.)
        for (int a = 0; a < myhf.nobasisfcn; a++) { //Sum over MO'er i (aj|kl)
            for (int j = 0; j < myhf.nobasisfcn; j++) { //Sum over basisfunktioner i (aj|kl)
                for (int k = 0; k < myhf.nobasisfcn; k++) { //Sum over basisfunktioner i (aj|kl)
                    for (int l = 0; l < myhf.nobasisfcn; l++) { //Sum over basisfunktioner i (aj|kl)
                        moeritensor[a][b][k][l] += myhf.cmx(j, b) * moeritensortmp[a][j][k][l];
                    }
                }
            }
        }
    }
    //Reset moeritensortmp
    for (int i = 0; i < myhf.nobasisfcn; i++) {
        for (int j = 0; j < myhf.nobasisfcn; j++) {
            for (int k = 0; k < myhf.nobasisfcn; k++) {
                for (int l = 0; l < myhf.nobasisfcn; l++) {
                    moeritensortmp[i][j][k][l] = 0;
                }
            }
        }
    }
    //Transformer tredje indeks (ab|kl) -> (ab|cl)
    for (int c = 0; c < myhf.nobasisfcn; c++) { //Sum over MO'er (rumorbitaler) (okk. + vir.)
        for (int a = 0; a < myhf.nobasisfcn; a++) { //Sum over MO'er i (ab|kl)
            for (int b = 0; b < myhf.nobasisfcn; b++) { //Sum over MO'er i (ab|kl)
                for (int k = 0; k < myhf.nobasisfcn; k++) { //Sum over basisfunktioner i (ab|kl)
                    for (int l = 0; l < myhf.nobasisfcn; l++) { //Sum over basisfunktioner i (ab|kl)
                        moeritensortmp[a][b][c][l] += myhf.cmx(k, c) * moeritensor[a][b][k][l];
                    }
                }
            }
        }
    }
    //Reset moeritensor
    for (int i = 0; i < myhf.nobasisfcn; i++) {
        for (int j = 0; j < myhf.nobasisfcn; j++) {
            for (int k = 0; k < myhf.nobasisfcn; k++) {
                for (int l = 0; l < myhf.nobasisfcn; l++) {
                    moeritensor[i][j][k][l] = 0;
                }
            }
        }
    }
    //Transformer fjerde indeks (ab|cl) -> (ab|cd)
    for (int d = 0; d < myhf.nobasisfcn; d++) { //Sum over MO'er (rumorbitaler) (okk. + vir.)
        for (int a = 0; a < myhf.nobasisfcn; a++) { //Sum over MO'er i (ab|cl)
            for (int b = 0; b < myhf.nobasisfcn; b++) { //Sum over MO'er i (ab|cl)
                for (int c = 0; c < myhf.nobasisfcn; c++) { //Sum over MO'er i (ab|cl)
                    for (int l = 0; l < myhf.nobasisfcn; l++) { //Sum over basisfunktioner i (ab|cl)
                        moeritensor[a][b][c][d] += myhf.cmx(l, d) * moeritensortmp[a][b][c][l];
                    }
                }
            }
        }
    }
    /********************************************************************/

    //Udregn andenordenskorrektionen til HF-energien, EMP2
    /********************************************************************/
    double EMP2 = 0;
    for (int a = 0; a < myhf.mymolecule.getnooccmos(); a++) { //Sum over okkuperede molekylorbitaler (rumorbitaler)
        for (int b = 0; b < myhf.mymolecule.getnooccmos(); b++) { //Sum over okkuperede molekylorbitaler (rumorbitaler)
            for (int r = myhf.mymolecule.getnooccmos(); r < myhf.nobasisfcn; r++) { //Sum over virtuelle molekylorbitaler (rumorbitaler)
                for (int s = myhf.mymolecule.getnooccmos(); s < myhf.nobasisfcn; s++) { //Sum over virtuelle molekylorbitaler (rumorbitaler)
                    EMP2 += (moeritensor[a][r][b][s]*(2*(moeritensor[a][r][b][s]) - moeritensor[a][s][b][r]))
                    / (myhf.moenergy[a] + myhf.moenergy[b] - myhf.moenergy[r] - myhf.moenergy[s]);
                }
            }
        }
    }
    /********************************************************************/

    //Slet moeritensor og moeritensortmp
    /********************************************************************/
    for (int i = 0; i < myhf.nobasisfcn; i++) { //Loop over rækker i den ydre 2D-matrix
        for (int j = 0; j < myhf.nobasisfcn; j++) { //Loop over kolonner i den ydre 2D-matrix
            for (int k = 0; k < myhf.nobasisfcn; k++) { //Loop over rækker i de indre 2D-matricer
                delete[] moeritensor[i][j][k];
                delete[] moeritensortmp[i][j][k];
            }
            delete[] moeritensor[i][j];
            delete[] moeritensortmp[i][j];
        }
        delete[] moeritensor[i];
        delete[] moeritensortmp[i];
    }
    delete[] moeritensor;
    delete[] moeritensortmp;
    /********************************************************************/
    
    time_t t2 = time(NULL);

    //Output resultatet
    /********************************************************************/
    outfile.precision(12);
    outfile << "E(MP2) / Hartrees:\n" << EMP2 << "\n\n";
    outfile << "Den korrigerede elektroniske energi / Hartrees:\n" << myhf.elecenergy + EMP2 << "\n\n";
    outfile << "Den korrigerede HF-energi / Hartrees:\n" << myhf.elecenergy + myhf.mymolecule.getncrepenergy() + EMP2 << "\n\n";

    outfile.precision(0);
    outfile << "**Beregningstid: Programmet var " << difftime(t2, t1) << " sekunder om at beregne MP2-energien med den hurtige algoritme.\n\n";
    /********************************************************************/
}

//--------------------------------------------------------------------------------

void dipolemoment(HF& myhf, ofstream& outfile)
{
    //Load dipolintegralerne
    myhf.loaddipoleintegrals();


    //Udregn x-, y- og z-komponenterne af dipolvektoren og det totale dipolmoment for molekylet
    //Dipolmomentet udregnes ud fra det elektroniske dipolmoment, som bestemmes her, og kernedipolmomentet som er gemt i mymolecule-datamedlemmet
    /********************************************************************/
    double dipolex = 0, dipoley = 0, dipolez = 0; //x-, y- og z-komponenterne af dipolvektoren
    //Først udregnes x-, y- og z-komponenterne af den elektroniske dipolvektor
    for (int i = 0; i < myhf.nobasisfcn; i++) {
        for (int j = 0; j < myhf.nobasisfcn; j++) {
            dipolex += myhf.dmx(i, j) * myhf.muxmx(i, j);
            dipoley += myhf.dmx(i, j) * myhf.muymx(i, j);
            dipolez += myhf.dmx(i, j) * myhf.muzmx(i, j);
        }
    }
    dipolex *= 2; dipoley *= 2; dipolez *= 2;
    //Tilføj kernedipolmomentet
    dipolex += myhf.mymolecule.getncdipolemomentx();
    dipoley += myhf.mymolecule.getncdipolemomenty();
    dipolez += myhf.mymolecule.getncdipolemomentz();
    /********************************************************************/


    //Output resultaterne
    /********************************************************************/
    outfile << "µ_x / au: " <<  dipolex << '\n';
    outfile << "µ_y / au: " <<  dipoley << '\n';
    outfile << "µ_z / au: " << dipolez << '\n';
    outfile << "\n";
    outfile << "µ_tot / au: " << sqrt(dipolex*dipolex + dipoley*dipoley + dipolez*dipolez) << '\n';
    outfile << "\n\n";
    /********************************************************************/
}

//--------------------------------------------------------------------------------