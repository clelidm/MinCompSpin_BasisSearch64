#include <iostream>
#include <map>
#include <list>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;


/******************************************************************************/
/**********************     CONSTANTS  and TOOLS  *****************************/
/******************************************************************************/
#include "data.h"
const uint64_t one64 = 1;

std::string int_to_bstring(__int128_t bool_nb, unsigned int n);
void int_to_digits(__int128_t bool_nb, unsigned int n);
void int_to_digits_file(__int128_t bool_nb, unsigned int r, fstream &file);

unsigned int bitset_count(__int128_t bool_nb);

/******************************************************************************/
/*******************   Convert  Basis  to  F2 Matrix   ************************/
/******************************************************************************/
list<unsigned int> RREF_F2(bool** M, int n, int m);

// This function returns a binary matrix that has the operators of the basis 'Basis' as columns:
// Each Operator is a column of the matrix
// n = Number of spins = number of rows --> 1rst index          //     !! We placed the lowest bit (most to the right) in the top row !!
// m = Number of basis operators = number of columns --> 2nd index

bool** Basis_to_MatrixF2(vector<Operator64> Basis, unsigned int n)
{
  unsigned int m = Basis.size();
  uint64_t Op = 0;

  // Create a Boolean Matrix:
  bool** M = (bool**) malloc(n*sizeof(bool*));  // n rows --> 1rst index
  for (int i=0; i<n; i++)
    {   M[i] = (bool*) malloc(m * sizeof(bool));  }  // m columns --> 2nd index

  unsigned int i = 0; //iteration over the row;
  unsigned int j = 0; //iteration over the columns;

  // Copy the basis operators:
  for (auto& it_Op : Basis)
  {
    Op = it_Op.bin;

    // filling in each column:
    for (i=0; i<n; i++)
    { 
      M[i][j] = Op & one64;
      Op >>= 1;
    }

    j++; // next column
  } 

  return M;
}

/******************************************************************************/
/************   Check if a set of operators are all independent   *************/
/******************************************************************************/
// return 'True' if all the element in the matrix are independent

bool Is_Basis(vector<Operator64> Basis, unsigned int n)
{
  cout << "-->> Check if the set of operators are independent:" << endl;
  cout << "\t Number of operators analysed: " << Basis.size() << endl;
  // Convert the basis to a boolean matrix:
  bool** M_Basis = Basis_to_MatrixF2(Basis, n);

  // Row Reduction procedure:
  list<unsigned int> list_lead = RREF_F2(M_Basis, n, Basis.size());

  // Free memory:
  for (int i=0; i<n; i++)
    {   free(M_Basis[i]);   } 
  free(M_Basis);

  // Return if the operators are independent or not:
  if (list_lead.size() == Basis.size())
    { 
      cout << "All the operators are independents." << endl << endl;
      return true;  
    }

  else
    { 
      cout << "Not all the operators are independents:";
      cout << "   (number of independent operators = " << list_lead.size() << ") < (number of operators = " << Basis.size() << " ) " << endl << endl;
      return false; 
    }
}

/******************************************************************************/
/****************   Print Terminal Vector Best Operators  *********************/
/******************************************************************************/
void PrintTerm_OpBasis(vector<Operator64> Basis, unsigned int n, unsigned int N)
{
  int i = 1;
  double p1 = 1, LogLi = 0, LogL = 0, Nd = (double) N;

  cout << "-->> Print Basis Operators: \t Total number of operators = " << Basis.size() << endl << endl;  

  for (auto& Op : Basis)
  {
    p1 = ((double) Op.k1) / Nd;
    LogLi = p1*log(p1)+(1-p1)*log(1-p1);
    LogL += LogLi;

    cout << i << "\t" << int_to_bstring(Op.bin, n) << "\t" << Op.bias << "\t" << Op.k1  << "\t" << p1 << "\t Indices = "; //<< endl;
    int_to_digits(Op.bin, n);
    i++;
  }
  cout << endl;
  cout << "##  LogL / N = " << LogL << endl;  //<< setprecision (6) 
  cout << "## -LogL / N / log(2) = " << -LogL/log(2.) << " bits per datapoints "<< endl;  //<< setprecision (6) 
  cout << endl;
}

void PrintFile_OpBasis(vector<Operator64> Basis, unsigned int n, unsigned int N, string filename)
{
  string OpSet_filename = OUTPUT_directory + filename + "_OpSet.dat";

  cout << "-->> Print Basis Operators in the file: \'" <<  OpSet_filename << "\'" << endl;
  fstream file_OpBasis(OpSet_filename, ios::out);
  
  int i = 1;
  double p1 = 1, LogLi = 0, LogL = 0, Nd = (double) N;

  file_OpBasis << "Basis: Total number of operators = " << Basis.size() << endl << endl;  

  for (auto& Op : Basis)
  {
    p1 = ((double) Op.k1) / Nd;
    LogLi = p1*log(p1)+(1-p1)*log(1-p1);
    LogL += LogLi;

    file_OpBasis << i << "\t" << int_to_bstring(Op.bin, n) << "\t" << Op.bias << "\t" << Op.k1  << "\t" << p1 << "\t Indices = "; //<< endl;
    int_to_digits_file(Op.bin, n, file_OpBasis);
    i++;
  }
  file_OpBasis << endl;
  file_OpBasis << "##  LogL / N = " << LogL << endl;  //<< setprecision (6) 
  file_OpBasis << "## -LogL / N / log(2) = " << -LogL/log(2.) << " bits per datapoints "<< endl;  //<< setprecision (6) 
  file_OpBasis.close();

  cout << endl;
}

/******************************************************************************/
/*********************   Histo Order of basis Operators  **********************/
/******************************************************************************/
map<unsigned int, unsigned int> Histo_BasisOpOrder(vector<Operator64> Basis)
{
  map<unsigned int, unsigned int> histo_order;

  for (auto& Op : Basis) {
    histo_order[bitset_count(Op.bin)]+=1;
  }

  cout << "Number of basis operators of a given order k:" << endl;
  for (auto& histo: histo_order){
    cout << "\t k = " << histo.first << " :\t" << histo.second << " operators" << endl;
  }
  cout << endl;

  return histo_order;
}






