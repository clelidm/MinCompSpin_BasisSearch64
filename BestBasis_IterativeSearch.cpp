#include <iostream>
#include <sstream>
#include <fstream>

#include <set>
#include <vector>
#include <list>

#include "data.h"

using namespace std;

/******************************************************************************/
/****************     Initial Choice of Operators for Basis    ****************/
/******************    All operators of order k or smaller    *****************/
/******************************************************************************/
set<Operator64> All_Op_k1(vector<pair<uint64_t, unsigned int>> Nvect, unsigned int n, unsigned int N, double *lowest_bias, bool print = false);

void Add_AllOp_kbits_MostBiased(set<Operator64>& OpSet, vector<pair<uint64_t, unsigned int>> Nvect, unsigned int n, unsigned int N, unsigned int k = 2, double Bias_LowerBound=0, bool print = false);

//void PrintTerm_OpSet(set<Operator64> OpSet);
void PrintTerm_OpSet(set<Operator64> OpSet, unsigned int n);
void PrintFile_OpSet(set<Operator64> OpSet, unsigned int n, string filename);

// Remove the Operator with too small Bias:
void CutSmallBias(set<Operator64>& OpSet, Struct_LowerBound LB);

unsigned int bitset_count(__int128_t bool_nb);
void int_to_digits_file(__int128_t bool_nb, unsigned int r, std::fstream &file);
std::string int_to_bstring(__int128_t bool_nb, unsigned int r);

/******************************************************************************/
/**************************     Select Best Basis    **************************/
/******************************************************************************/
vector<Operator64> BestBasis_inOpSet(set<Operator64> OpSet, unsigned int n, Struct_LowerBound* LowerBound, unsigned int m=1000);

/******************************************************************************/
/**************************     Basis  Tools  *********************************/
/******************************************************************************/
bool Is_Basis(vector<Operator64> Basis, unsigned int n);

void PrintTerm_OpBasis(vector<Operator64> OpVect_Basis, unsigned int n, unsigned int N);
void PrintFile_OpBasis(vector<Operator64> OpVect_Basis, unsigned int n, unsigned int N, string filename);

bool Check_Basis_Identity(vector<Operator64> Basis)
{
  bool check = true;

  for (auto& Op : Basis) {
    if(bitset_count(Op.bin) != 1)   { check = false; }
  }
  return check;
}

void SaveFile_Basis(vector<Operator64> Basis, unsigned int n, fstream &file)
{
  file << "########### New Basis: \t Total number of operators = " << Basis.size() << endl;  

  for (auto& Op : Basis)
  {
    file <<  int_to_bstring(Op.bin, n) << "\t Bias = " << Op.bias << "\t"; 
    int_to_digits_file(Op.bin, n, file);
  }

  file << endl;
}

/******************************************************************************/
/***************     Search in a Given Representation  Tools  *****************/
/******************************************************************************/

vector<Operator64> BestBasisSearch_FixedRepresentation(vector<pair<uint64_t, unsigned int>> Nvect, unsigned int n, unsigned int N, unsigned int k_max, unsigned int B_it, bool bool_print = false)
{
  k_max = (k_max<2)?2:k_max;  // k_max must be at least 2;

  cout << endl << "*******************  FIND THE SMALLEST BIAS OF THE CURRENT BASIS:  ************************";
  cout << endl << "*******************************************************************************************" << endl;

  // Compute the bias of the basis elements, and find the least biased one:
  // This value will serve as a lower bound for operators that will be kept later on.

  double Bias_LowerBound = 0.;  // Current lower bound (current lowest bias) is 0. --> we accept all possible bias

  set<Operator64> OpSet = All_Op_k1(Nvect, n, N, &Bias_LowerBound, bool_print); 

  //PrintTerm_OpSet(OpSet_B0, n);
  PrintFile_OpSet(OpSet, n, "B" + to_string(B_it) +"_k1");

  Struct_LowerBound LB; // Lower Bound Info
  LB.Bias = Bias_LowerBound;
  vector<Operator64> BestBasis;

/*  for (auto& Op: OpSet){
    BestBasis.push_back(Op);
  }
  cout << "Check Basis Identity: k=1: " << Check_Basis_Identity(BestBasis) << endl << endl; 
*/

  string filename_k = "";

  for (unsigned int k = 2; k <= k_max; k++)
  {
      filename_k = "B" + to_string(B_it) +"_k"+to_string(k);

      cout << endl << "****************************  ADD ALL OPERATORS for k = " << k << "  ********************************";
      cout << endl << "*******************************************************************************************" << endl;

      Add_AllOp_kbits_MostBiased(OpSet, Nvect, n, N, k, LB.Bias, bool_print);

      //PrintTerm_OpSet(OpSet, n);
      PrintFile_OpSet(OpSet, n, filename_k);

      cout << endl << "******************************  SEARCH FOR BEST BASIS:  ***********************************" << endl;
//    cout << endl << "*******************************************************************************************" << endl;

      BestBasis.clear();
      BestBasis = BestBasis_inOpSet(OpSet, n, &LB, 1000); // LB will be over-written with the updated values

      PrintFile_OpBasis(BestBasis, n, N, filename_k + "_BestBasis");

      CutSmallBias(OpSet, LB);

      PrintFile_OpSet(OpSet, n, filename_k + "_CutSmallBias");
  }

  cout << endl << "*****************************  DONE: SEARCH IN GIVEN REPRESENTATION  ****************************" << endl;

  return BestBasis;
}


/******************************************************************************/
/****************     Search in DIFFERENT REPRESENTATIONS   *******************/
/******************************************************************************/
vector<pair<uint64_t, unsigned int>> build_Kvect(vector<pair<uint64_t, unsigned int>> Nvect, list<uint64_t> Basis);

vector<Operator64> BestBasisSearch_Final(vector<pair<uint64_t, unsigned int>> Nvect, unsigned int n, unsigned int N, unsigned int k_max, bool bool_print = false)
{
    k_max = (k_max<2)?2:k_max;  // k_max must be at least 2;

    cout << endl << "*******************************************************************************************";
    cout << endl << "***************************  SEARCH IN THE ORIGINAL BASIS:  *******************************";
    cout << endl << "*******************************************************************************************" << endl;

    unsigned int B_it = 0;   // initial basis

    vector<Operator64> BestBasis = BestBasisSearch_FixedRepresentation(Nvect, n, N, k_max, B_it, bool_print);

//Save Basis:
    string Basis_filename = OUTPUT_directory + "All_Bases.dat";
    fstream Basis_file(Basis_filename, ios::out); 

    Basis_file << "### File containing all the successive Basis" << endl;
    Basis_file << "### Note that Bases are given in the successive representation, and not in the original representation" << endl << endl;

    PrintTerm_OpBasis(BestBasis, n, N);   
    SaveFile_Basis(BestBasis, n, Basis_file);

    cout << endl << "*******************************************************************************************";
    cout << endl << "***************************  SEARCH IN THE NEW BASIS:  ************************************";
    cout << endl << "*******************  Stops when the found basis is Identity  ******************************";
    cout << endl << "*******************************************************************************************" << endl;

    bool check = false;

    while( check == false ) // if the best basis is not the identity: then continue changing representation
    {
      cout << "-->> Change the representation of the data in the current Best Basis:" << endl;
      list<uint64_t> Basis_li;
      for(auto& Op:BestBasis)  { Basis_li.push_back(Op.bin);  }  // extract the integer representation of the basis operators:
    
      vector<pair<uint64_t, unsigned int>> Kvect = build_Kvect(Nvect, Basis_li);

      B_it = 1;   // New basis

      BestBasis.clear();
      BestBasis = BestBasisSearch_FixedRepresentation(Kvect, n, N, k_max, B_it, bool_print);

      PrintTerm_OpBasis(BestBasis, n, N);  
      SaveFile_Basis(BestBasis, n, Basis_file);

      check = Check_Basis_Identity(BestBasis);

      cout << "Check Basis Identity: k=1: " << Check_Basis_Identity(BestBasis) << endl << endl; 
    }

    Basis_file.close();

    cout << "-->> All successive Bases are saved in the file: \'" <<  Basis_filename << "\'" << endl;
    cout << "Note that Bases are given in the successive representation, and not in the original representation" << endl << endl;

    return BestBasis;
}






