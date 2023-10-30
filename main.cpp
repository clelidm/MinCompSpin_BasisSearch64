//g++ -std=c++11 main.cpp ReadDataFile.cpp Tool_BitOperations.cpp BestBasis_Init_SetOp.cpp

#include <iostream>
#include <sstream>

#include <map>
#include <set>
#include <list>
#include <vector>

#include <ctime> // for chrono
#include <ratio> // for chrono
#include <chrono> // for chrono

#include "data.h"

using namespace std;

/********************************************************************/
/**************************    PARAMETERS    ************************/
/********************************************************************/
// number of binary (spin) variables:
unsigned int n = 50;

// INPUT DATA FILES (optional):  
// the input datafile can also be specified directly in the main() function, as an argument of the function "read_datafile()":
std::string input_datafile = "INPUT/Big5-IPC1_VS3_Ne5.dat"; //Big5PT.sorted_Ne5";  //"INPUT/Big5PT.sorted";  //"INPUT/SCOTUS_n9_N895_Data.dat"; //

/******************************************************************************/
/**************************     READ FILE    **********************************/
/******************************************************************************/
/**************    READ DATA and STORE them in Nset    ************************/
vector<pair<uint64_t, unsigned int>> read_datafile64_vect(string datafilename, unsigned int *N, unsigned int r);

/******************************************************************************/
/**************************     Basis  Tools  *********************************/
/******************************************************************************/
bool Is_Basis(vector<Operator64> Basis, unsigned int n);

void PrintTerm_OpBasis(vector<Operator64> OpVect_Basis, unsigned int n, unsigned int N);
void PrintFile_OpBasis(vector<Operator64> OpVect_Basis, unsigned int n, unsigned int N, string filename);

map<unsigned int, unsigned int> Histo_BasisOpOrder(vector<Operator64> Basis);

/******************************************************************************/
/************************   BASIS SEARCH TOOLS    *****************************/
/******************************************************************************/
// Exhaustive Search:
vector<Operator64> BestBasis_ExhaustiveSearch(vector<pair<uint64_t, unsigned int>> Nvect, unsigned int n, unsigned int N, bool bool_print = false);

// Fixed Representation up to order `k_max``:
vector<Operator64> BestBasisSearch_FixedRepresentation(vector<pair<uint64_t, unsigned int>> Nvect, unsigned int n, unsigned int N, unsigned int k_max, unsigned int B_it, bool bool_print = false);

// Changing representation up to order `k_max``:
vector<Operator64> BestBasisSearch_Final(vector<pair<uint64_t, unsigned int>> Nvect, unsigned int n, unsigned int N, unsigned int k_max, bool bool_print = false);

/******************************************************************************/
/************************** MAIN **********************************************/
/******************************************************************************/

int main(int argc, char *argv[])
{
	/**********************     READ ARGUMENTS    *********************************/

    // argv[0] contains the name of the datafile, from the current folder (i.e. from the folder containing "data.h");
    // argv[1] contains the number of variables to read;
    if (argc == 3)
    {
      	string input_datafile_buffer = argv[1];
        string n_string_buffer = argv[2];

        input_datafile = "INPUT/" + input_datafile_buffer;
        //n_string_buffer = argv[2];
        n = stoul(n_string_buffer);
    }
    else if (argc != 1)
    {
        cout << "The number of arguments must be either 0 or 2" << endl;
        return 0;
    }

    cout << "--->> Create the \"OUTPUT\" Folder: (if needed) ";
    system(("mkdir -p " + OUTPUT_directory).c_str());
    cout << endl;

    // chrono variables:
	auto start = chrono::system_clock::now(); 
	auto end = chrono::system_clock::now(); 
    chrono::duration<double> elapsed = end - start; 

    cout << endl << "*******************************************************************************************";
    cout << endl << "***********************************  READ THE DATA:  **************************************";
    cout << endl << "*******************************************************************************************" << endl;

	unsigned int N=0;  // will contain the number of datapoints in the dataset

    vector<pair<uint64_t, unsigned int>> Nvect = read_datafile64_vect(input_datafile, &N, n); 

	if (N == 0) { return 0; } // Terminate program if the file can't be found or is empty


    cout << endl << "*******************************************************************************************";
    cout << endl << "*******************************************************************************************";
    cout << endl << "***************************  SEARCH IN THE ORIGINAL BASIS:  *******************************";
    cout << endl << "*******************************************************************************************";
    cout << endl << "*******************************************************************************************" << endl;

    bool bool_print = false;
    unsigned int k_max = 2;  // largest order of operators to take into account in each representation
    unsigned int B_it = 0;   // initial basis

    vector<Operator64> BestBasis = BestBasisSearch_FixedRepresentation(Nvect, n, N, k_max, B_it, bool_print);

    cout << endl << "*******************************************************************************************";
    cout << endl << "***********************  PRINT TERMINAL/FILE FINAL OPERATOR SET:  *************************";
    cout << endl << "*******************************************************************************************" << endl;

    PrintTerm_OpBasis(BestBasis, n, N);  
    //Is_Basis(BestBasis_k2, n);   // this function can check if a set of Operators is in independent set

    Histo_BasisOpOrder(BestBasis);

    cout << endl << "*******************************************************************************************";
    cout << endl << "*******************************************************************************************";
    cout << endl << "**********************************  EXHAUSTIVE SEARCH:  ***********************************";
    cout << endl << "*******************************************************************************************";
    cout << endl << "*******************************************************************************************" << endl;

    cout << "The Exhaustive Search is not recommended for dataset with more than n~20 variables." << endl << endl;

    BestBasis = BestBasis_ExhaustiveSearch(Nvect, n, N, bool_print);

    PrintTerm_OpBasis(BestBasis, n, N);  


    cout << endl << "*******************************************************************************************";
    cout << endl << "*******************************************************************************************";
    cout << endl << "**********************  SEARCH IN DIFFERENT REPRESENTATIONS:  *****************************";
    cout << endl << "********************  ! Stops when the found basis is Identity !  **************************";
    cout << endl << "*******************************************************************************************";
    cout << endl << "*******************************************************************************************" << endl;
/*
    k_max = 4;
    vector<Operator64> BestBasis = BestBasisSearch_Final(Nvect, n, N, k_max, bool_print);
*/

    return 0;
}




