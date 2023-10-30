#include <iostream>
#include <map>
#include <cmath>
#include <set>
#include <vector>
#include <fstream>

#include <ctime> // for chrono
#include <ratio> // for chrono
#include <chrono> // for chrono

using namespace std;


/******************************************************************************/
/**********************     CONSTANTS  and TOOLS  *****************************/
/******************************************************************************/
#include "data.h"
const uint64_t one64 = 1;

unsigned int bitset_count(__int128_t bool_nb);
string int_to_bstring(__int128_t bool_nb, unsigned int r);

void int_to_digits(__int128_t bool_nb, unsigned int r);
void int_to_digits_file(__int128_t bool_nb, unsigned int r, fstream &file);

//double min_double(double a, double b)
//{  return !(b<a)?a:b;   }

unsigned int Choose(unsigned int n, unsigned int k);

/******************************************************************************/
/************************   Print Terminal Operators  *************************/
/******************************************************************************/
void PrintTerm_OpSet(set<Operator64> OpSet, unsigned int n)
{
  cout << "--> Print Set of Operators: \t Total number of operators = " << OpSet.size() << endl << endl;  

  for (auto& Op : OpSet)
  {
    cout <<  int_to_bstring(Op.bin, n) << "\t Bias = " << Op.bias << "\t"; // << endl;
    int_to_digits(Op.bin, n);
  }
}

void PrintFile_OpSet(set<Operator64> OpSet, unsigned int n, string filename)
{
  string OpSet_filename = OUTPUT_directory + filename + "_OpSet.dat";

  cout << "-->> Print Operators in the file: \'" <<  OpSet_filename << "\'" << endl;

  fstream file_OpSet(OpSet_filename, ios::out);

  file_OpSet << "--> Print Set of Operators: \t Total number of operators = " << OpSet.size() << endl << endl;  

  for (auto& Op : OpSet)
  {
    file_OpSet <<  int_to_bstring(Op.bin, n) << "\t Bias = " << Op.bias << "\t"; //<< endl;
    int_to_digits_file(Op.bin, n, file_OpSet);
  }
  file_OpSet.close();
}

/******************************************************************************/
/********************     AVERAGES and OBSERVABLES   **************************/
/******************************************************************************/
// Number of times an operator is equal to 1 ( = <phi> in the {0,1} representation ) in the dataset
unsigned int K1_Op(vector<pair<uint64_t, unsigned int>> Nvect, uint64_t Op)  // Complexity = O(|Nset|)
{
  unsigned int K1=0;

  for (auto& it : Nvect)
    {    K1 += (bitset_count( ((it).first) & Op ) % 2)*((it).second);   } 

  return K1;
}

// ******* Data averages are taken using ISING convention: ******************** / 
Operator64 Value_Op(uint64_t Op_bin, vector<pair<uint64_t, unsigned int>> Nvect, double Nd)
{
  Operator64 Op;

  Op.bin = Op_bin;
  Op.k1 = K1_Op(Nvect, Op.bin);

  //Op.S = (Op.p1_D==0 || Op.p1_D==1)? 0: -( (Op.p1_D*log(Op.p1_D) + (1.-Op.p1_D)*log(1.-Op.p1_D)) ); // Entropy of Op in the data
  //Op.DKL = log(2.)-Op.S;   // DK = max-log-likelihood for each Op, translated by (log2) and divided by N

  //double p1 = ((double) Op.k1) / Nd;
  Op.bias = fabs((((double) Op.k1) / Nd) -0.5);  //fabs(p1-0.5);

  return Op;
}

/******************************************************************************/
/***********************   All Operators with 1 bit only  *********************/
/************************   Find the lowest bias value  ***********************/
/******************************************************************************/
// this value will serve as a lower bound for operators that we will keep later on.
set<Operator64> All_Op_k1(vector<pair<uint64_t, unsigned int>> Nvect, unsigned int n, unsigned int N, double *lowest_bias, bool print = false)
{
  auto start = chrono::system_clock::now();

  set<Operator64> OpSet;
  Operator64 Op;

  double Nd = (double) N;
  (*lowest_bias) = 0.5;

  cout << "--> Compute and rank all the observables of order 1 (fields).. " << endl; 

  uint64_t un_i = 1;
  for (int i=0; i<n; i++) // All Fields:
  { 
    Op = Value_Op(un_i, Nvect, Nd);
    OpSet.insert(Op);
    //if (Op.bias < (*lowest_bias)) { (*lowest_bias) = Op.bias; }
    if(print)
    { 
      cout <<  int_to_bstring(Op.bin, n) << "\t Bias = " << Op.bias << "\t";
      int_to_digits(Op.bin, n); 
    }
    un_i = un_i << 1;
  }

  cout << endl;

  cout << "Smallest bias = " << (*OpSet.rbegin()).bias << ", \t for the operator = " << int_to_bstring((*OpSet.rbegin()).bin, n) << "\t";
  int_to_digits((*OpSet.rbegin()).bin, n);
  cout << " Largest bias = " << (*OpSet.begin()).bias << ", \t for the operator = " << int_to_bstring((*OpSet.begin()).bin, n) << "\t";
  int_to_digits((*OpSet.begin()).bin, n);
  cout << endl;

  (*lowest_bias) = (*OpSet.rbegin()).bias;

// Compare with 3*sigma bound:
  double sigma3_bound = 3 * 0.5 / sqrt(Nd);
  if (sigma3_bound < (*lowest_bias)) {  cout << "Taking the \'small bias\' as a bound is MORE restrictive than the 3*sigma bound = " << sigma3_bound << endl; }
  else {  cout << "Taking the \'small bias\' as a bound is LESS restrictive than the 3*sigma bound = " << sigma3_bound << endl; }

// Time:
  auto end = chrono::system_clock::now();  chrono::duration<double> elapsed = end - start;
  cout << endl << "Elapsed time (in s): " << elapsed.count() << endl << endl; 

// Estimated time for next values of 'k':
  cout << "Estimated times (in s) for larger values of \'k\':" << endl;

  double time_estimate = 0;
  for (unsigned int k=2; k<6; k++) {
    time_estimate = elapsed.count()/n*Choose(n, k);
    cout << "\t k = " << k << ": \t" << time_estimate << " s \t" << time_estimate/60. << " min" << endl;
  }
  cout << endl;

  return OpSet;
}

/******************************************************************************/
/*************************   All Integers with k-bits   ***********************/
/******************************************************************************/
bool Incr_k_bits(unsigned int k, uint64_t *a, unsigned int n)  //only for (k <= n)
{
    uint64_t c = 1, c_stop = 1;

// Find position of the lowest bit:
    c = (((*a) - 1) ^ (*a)) & (*a);

// Increament "a" by 1 starting from the lowest bit:
    *a = *a + c;

// Fill with 1's, starting from bit 0, until "a" has k-bits:
    while(bitset_count(*a) < k)
        {   (*a)++;  c = 1;     }    // c = 1 --> lowest bit is now in bit 0
  
// Print:
    //std::string abit = int_to_bstring(*a, n);
    //std::cout << "final a = " << abit << "\t c = " << c << std::endl;

// Stopping criteria:  position of the lowest bit before taking a = a + c,  if it is "(n-k-1)", then we're done! --> stop
//if i_min = (n-k-1) the procedure is finished  
    c_stop = one64 << (n-k-1);

    return (c == c_stop); //((c & c_stop)); // return TRUE when the procedure is done --> STOP
}


void all_int_k_bits(unsigned int k, uint32_t *compt, unsigned int n)  //only for (k <= n)
{
    uint64_t a = (one64 << k) - 1;  // intialise "a" with the k first bits set at 1
    bool stop = false;

    (*compt)++;
    //std::cout << *compt << ": \t" << int_to_bstring(a, n) << std::endl;

    while (!stop)
    {
        stop = Incr_k_bits(k, &a, n);    
        (*compt)++;   
        //std::cout << *compt << "\t" << stop << ": \t" << int_to_bstring(a, n) << std::endl; 
    }

    //std::cout << std::endl << "total number of combinations = " << compt << std::endl;
}

/******************************************************************************/
/*********************   Add all Operators with 2 bit(s)  *********************/
/***************   All fields and all pairwise interactions   *****************/
/******************************************************************************/

void Add_AllOp_kbits_MostBiased(set<Operator64>& OpSet, vector<pair<uint64_t, unsigned int>> Nvect, unsigned int n, unsigned int N, unsigned int k = 2, double Bias_LowerBound=0, bool print = false)  
{
  auto start = chrono::system_clock::now(); 

  cout << "Current smallest 'Bias' = " << Bias_LowerBound ;
  cout << "\t --> all operator with smaller 'Bias' will be rejected" << endl;

  double Nd = (double) N;
  unsigned int OpSet_Size0 = OpSet.size();

  cout << "Start iteration: k = " << k << endl;

// intialise "Op_bin" with the "k" first bits set at 1 --> first operator:
  uint64_t Op_bin = (one64 << k) - 1; 
  Operator64 Op = Value_Op(Op_bin, Nvect, Nd);
  if (Op.bias > Bias_LowerBound)  { OpSet.insert(Op); }
  uint32_t compt = 1;

  if(print) {   
    cout << int_to_bstring(Op.bin, n) << "\t Bias = " << Op.bias << "\t";
    int_to_digits(Op.bin, n); 
    } //", \t Elapsed time (in s): " << elapsed.count() << endl; 

  bool stop = false;
  while (!stop)
  {
    stop = Incr_k_bits(k, &Op_bin, n); 
    Op = Value_Op(Op_bin, Nvect, Nd);
    if (Op.bias > Bias_LowerBound) { OpSet.insert(Op); }  
    if(print) {   
      cout << int_to_bstring(Op.bin, n) << "\t Bias = " << Op.bias << "\t";
      int_to_digits(Op.bin, n); 
    } //", \t Elapsed time (in s): " << elapsed.count() << endl; 
    compt++;
  } 

  cout << "End iteration: k = " << k << "\t total number of combinations = " << compt << "\t total number of accepted operators = " << OpSet.size() - OpSet_Size0 << endl;

  auto end = chrono::system_clock::now();  
  chrono::duration<double> elapsed = end - start;
  cout << endl << "Elapsed time (in s): " << elapsed.count() << endl << endl;  

// Estimated time for next values of 'k':
  cout << "Estimated times (in s) for larger values of \'k\':" << endl;

  double time_estimate = 0;
  double ch_k = Choose(n, k);

  for (unsigned int kk=k; kk<k+4; kk++) {
    time_estimate = elapsed.count() / ch_k * Choose(n, kk);
    cout << "\t k = " << kk << ": \t" << time_estimate << " s \t" << time_estimate/60. << " min" << endl;
  }
  cout << endl;
}


/******************************************************************************/
/********************   All Operators with k bits or less  ********************/
/*******************   All interactions of order k or less   ******************/
/******************************************************************************/
/*
set<Operator64> all_Op_UpTo_k_bits_MostBiased_rank64(vector<pair<uint64_t, unsigned int>> Nvect, unsigned int n, unsigned int N, unsigned int k, double Bias_LowerBound=0, bool print = false) 
{
  set<Operator64> OpSet; // Set of all the fields and pairwise operators ordered by bias 
  Operator64 Op;

  uint64_t Op_bin = 0;
  bool stop = false;
  double Nd = (double) N;
  //double sig_lim = alpha * 0.5 / sqrt(Nd);

  for (unsigned int ki=1; ki<=k; ki++)
  {
    Op_bin = (one64 << ki) - 1;  // intialise "a" with the ki first bits set at 1
    Op = Value_Op(Op_bin, Nvect, Nd);
    if (Op.bias > Bias_LowerBound) { OpSet.insert(Op); }

    //std::cout << int_to_bstring(Op_bin, n) << std::endl;
    if(print)
      {   std::cout << "Starts ki = " << ki << endl;  }

    stop = false;
    while (!stop)
      {
        stop = Incr_k_bits(k, &Op_bin, n); 
        Op = Value_Op(Op_bin, Nvect, Nd);
        if (Op.bias > Bias_LowerBound) { OpSet.insert(Op); }
        //std::cout << stop << ": \t" << int_to_bstring(Op_bin, n) << std::endl; 
      } 
  }
  return OpSet;
}
*/

/******************************************************************************/
/********************   REMOVE OPERATORS with SMALL BIAS  *********************/
/******************************************************************************/
void CutSmallBias(set<Operator64>& OpSet, Struct_LowerBound LB)
{
  cout << "-->> Remove operators with small bias:" << endl;
  cout << "\t Smallest Bias accepted = " << LB.Bias << endl; 
  cout << "\t Number of Operators left = " << (LB.Index + 1) << endl;  // indexing of the operators starts from '0', hence the '+1'

  set<Operator64>::iterator it = OpSet.begin();
  advance(it, (LB.Index+1));
  OpSet.erase (it, OpSet.end());

  cout << endl;
}

