#include <iostream>
//#include <sstream>
//#include <fstream>

#include <set>
#include <vector>

#include <ctime> // for chrono
#include <ratio> // for chrono
#include <chrono> // for chrono

using namespace std;

#include "data.h"

/******************************************************************************/
/***************************   Constant variables   ***************************/
/******************************************************************************/
const __int128_t one128 = 1;
const uint64_t one64 = 1;

/******************************************************************************/
/***********************   All Operators with 1 bit only  *********************/
/************************   Find the lowest bias value  ***********************/
/******************************************************************************/

set<Operator64> All_Op_k1(vector<pair<uint64_t, unsigned int>> Nvect, unsigned int n, unsigned int N, double *lowest_bias, bool print = false);

Operator64 Value_Op(uint64_t Op_bin, vector<pair<uint64_t, unsigned int>> Nvect, double Nd);

vector<Operator64> BestBasis_inOpSet(set<Operator64> OpSet, unsigned int n, Struct_LowerBound* LowerBound, unsigned int m=1000);

/******************************************************************************/
/*******************************   All Operators  *****************************/
/****************   Keep only the one with bias larger than LB  ***************/
/******************************************************************************/

set<Operator64> All_Op_LBk1 (vector<pair<uint64_t, unsigned int>> Nvect, unsigned int n, unsigned int N, bool print = false)
{
  double lowest_bias = 0;

  set<Operator64> OpSet = All_Op_k1(Nvect, n, N, &lowest_bias, print);
  Operator64 Op;
  double Nd = (double) N;

  cout << "-->> Compute and rank all Operators:" << endl;

  unsigned int Op_bin_max =  (one64 << n) - 1;

  for (unsigned int Op_bin = 1; Op_bin <= Op_bin_max; Op_bin++)
  {
    Op = Value_Op(Op_bin, Nvect, Nd);
    if (Op.bias > lowest_bias) { OpSet.insert(Op); } 
  }

  return OpSet;
}

/******************************************************************************/
/***************************   Exhaustive Search  *****************************/
/******************************************************************************/

vector<Operator64> BestBasis_ExhaustiveSearch(vector<pair<uint64_t, unsigned int>> Nvect, unsigned int n, unsigned int N, bool bool_print = false)
{
  auto start = chrono::system_clock::now();

  cout << endl << "*******************************************************************************************";
  cout << endl << "************************  EXHAUSTIVE SEARCH FOR THE BEST BASIS:  **************************";
  cout << endl << "*******************************************************************************************" << endl;

  cout << "-->> Compute all Operators, with a smallest accepted biased fixed by the least informative first order operator:" << endl;
  set<Operator64> OpSet = All_Op_LBk1 (Nvect, n, N, bool_print);

// Time:
  auto end = chrono::system_clock::now();  chrono::duration<double> elapsed = end - start;
  cout << "\t Elapsed time (in s): " << elapsed.count() << endl << endl; 

  cout << "-->> Search for the Best Basis:" << endl;

  Struct_LowerBound LB; // Lower Bound Info
  LB.Bias = 0;
  vector<Operator64> BestBasis = BestBasis_inOpSet(OpSet, n, &LB, 1000); // LB will be over-written with the updated values

// Time:
  end = chrono::system_clock::now();  elapsed = end - start;
  cout << "\t Total elapsed time (in s): " << elapsed.count() << endl << endl; 

  return BestBasis;
}