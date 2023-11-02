#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <list>

#include <ctime> // for chrono
#include <ratio> // for chrono
#include <chrono> // for chrono

using namespace std;

/******************************************************************************/
/***************************   Constant variables   ***************************/
/******************************************************************************/
const __int128_t one128 = 1;
const uint64_t one64 = 1;

unsigned int bitset_count(__int128_t bool_nb);

/******************************************************************************/
/**************************     READ FILE    **********************************/
/******************************************************************************/
/**************    READ DATA and STORE them in Nset    ************************/

vector<pair<uint64_t, unsigned int>> read_datafile64_vect(string datafilename, unsigned int *N, unsigned int r)    // O(N)  where N = data set size
{
  auto start = chrono::system_clock::now();

  cout << endl << "--->> Read the datafile: \"" << datafilename << "\", \t Build Nset..." << endl;
  cout << "\t Number of variables to read: n = " << r << endl;

  string line, line2;     char c = '1';
  uint64_t state = 0, Op;
  (*N) = 0;            // N = dataset sizes

// ***** The data is stored in Nset as an histogram:  ********************************
  map<uint64_t, unsigned int> Nset; // Nset[mu] = #of time state mu appears in the data set

  ifstream myfile (datafilename.c_str());
  if (myfile.is_open())
  {
    while ( getline (myfile,line))
    {
      line2 = line.substr (0,r);          //take the r first characters of line
      Op = one64 << (r - 1);
      state = 0;
      for (auto &elem: line2)     //convert string line2 into a binary integer
      {
        if (elem == c) { state += Op; }
        Op = Op >> 1;
      }
      Nset[state] += 1;
      //cout << line << endl;   //cout << state << " :  " << int_to_bstring(state, r) << endl;
      (*N)++;
    }
  }
  else cout << endl << "--->> Unable to open file: Check datafilename and location." << endl << endl;
  
  myfile.close();

  if ((*N) == 0) 
    { 
    cout << endl << "--->> Failure to read the file, or file is empty:  Terminate." << endl << endl;
    }
  else
    {
    cout << endl << "--->> File has been read successfully:" << endl;
    cout << "\t Data size, N = " << (*N) << endl;
    cout << "\t Number of different states, Nset.size() = " << Nset.size() << endl << endl;
    }

  vector<pair<uint64_t, unsigned int>> Nvect(Nset.size());
  int i=0;
  for (auto& my_pair : Nset)
  {
    Nvect[i]=my_pair;
    i++;
  }

  auto end = chrono::system_clock::now();  
  chrono::duration<double>  elapsed = end - start;
  cout << endl << "Elapsed time (in s): " << elapsed.count() << endl << endl; 

  return Nvect;
}


vector<pair<__int128_t, unsigned int>> read_datafile128_vect(string datafilename, unsigned int *N, unsigned int r)    // O(N)  where N = data set size
{
  auto start = chrono::system_clock::now();

  cout << endl << "--->> Read the datafile: \"" << datafilename << "\", \t Build Nset..." << endl;
  cout << "\t Number of variables to read: n = " << r << endl;

  string line, line2;     char c = '1';
  __int128_t state = 0, Op;
  (*N) = 0;            // N = dataset sizes

// ***** The data is stored in Nset as an histogram:  ********************************
  map<__int128_t, unsigned int> Nset; // Nset[mu] = #of time state mu appears in the data set

  ifstream myfile (datafilename.c_str());
  if (myfile.is_open())
  {
    while ( getline (myfile,line))
    {
      line2 = line.substr (0,r);          //take the r first characters of line
      Op = one128 << (r - 1);
      state = 0;
      for (auto &elem: line2)     //convert string line2 into a binary integer
      {
        if (elem == c) { state += Op; }
        Op = Op >> 1;
      }
      Nset[state] += 1;
      (*N)++;
    }
    myfile.close();
  }
  else cout << endl << "--->> Unable to open file: Check datafilename and location." << endl << endl;

  if ((*N) == 0) 
    { 
    cout << endl << "--->> Failure to read the file, or file is empty:  Terminate." << endl << endl;
    }
  else
    {
    cout << endl << "--->> File has been read successfully:" << endl;
    cout << "\t Data size, N = " << (*N) << endl;
    cout << "\t Number of different states, Nset.size() = " << Nset.size() << endl << endl;
    }

  vector<pair<__int128_t, unsigned int>> Nvect(Nset.size());
  int i=0;
  for (auto& my_pair : Nset)
  {
    Nvect[i]=my_pair;
    i++;
  }

  auto end = chrono::system_clock::now();  
  chrono::duration<double>  elapsed = end - start;
  cout << endl << "Elapsed time (in s): " << elapsed.count() << endl << endl;  

  return Nvect;
}

/****************    PRINT Nset in file:    ************************/
/*void read_Nset (map<uint32_t, unsigned int> Nset, unsigned int N, string OUTPUTfilename)
// map.second = nb of time that the state map.first appears in the data set
{
  map<uint32_t, unsigned int>::iterator it;
  int Ncontrol = 0;

  fstream file(OUTPUTfilename.c_str(), ios::out);
  file << "#N = " << N << endl;
  file << "#Total number of accessible states = " << NOp_tot << endl;
  file << "#Number of visited states, Nset.size() = " << Nset.size() << endl;
  file << "#" << endl;
  file << "#1: state \t #2: nb of pts in state \t #3: Pba state" << endl;

  for (it = Nset.begin(); it!=Nset.end(); ++it)
  {
    file << it->first << ":\t" << bitset<n>(it->first) << " => " << it->second; // << endl;
    file << "  \t  P = " << it->second / (float) N << endl;
    Ncontrol += it->second;
  }

  if (Ncontrol != N) { cout << "Error function \'read_Nset\': Ncontrol != N" << endl;  }

  file.close();
}
*/

/******************************************************************************/
/*********************     CHANGE of BASIS: one datapoint  ********************/
/******************************************************************************/
// Given a choice of a basis (defined by the m-basis list) --> returns the new m-state (i.e. state in the new m-basis)
// Rem: must have m <= n 

// mu = old state
// final_mu = new state

uint64_t transform_mu_basis(uint64_t mu, list<uint64_t> basis)
{
  uint64_t un_i = 1, proj;
  uint64_t final_mu = 0;

  list<uint64_t>::iterator phi_i;

  for(phi_i = basis.begin(); phi_i != basis.end(); ++phi_i)
  {
    proj = (*phi_i) & mu;
    if ( (bitset_count(proj) % 2) == 1) // odd number of 1, i.e. sig_i = 1
    {
      final_mu += un_i;
    }
    un_i = (un_i << 1);
  }

  return final_mu;
}

/******************************************************************************/
/******************** CHANGE of BASIS: build K_SET ****************************/
/******************************************************************************/
// Build Kvect for the states written in the basis of the m-chosen independent 
// operator on which the SC model is based:

vector<pair<uint64_t, unsigned int>> build_Kvect(vector<pair<uint64_t, unsigned int>> Nvect, list<uint64_t> Basis)
// sig_m = sig in the new basis and cut on the m first spins 
// Kvect[sig_m] = #of time state mu_m appears in the data set
{
    map<uint64_t, unsigned int > Kvect_map;
    uint64_t sig_m;    // transformed state and to the m first spins

// ***** Build Kvect: *************************************************************************************
    cout << endl << "--->> Build Kvect..." << endl;

    for (auto const& it : Nvect)
    {
        sig_m = transform_mu_basis((it).first, Basis); // transform the initial state s=(it).first into the new basis
        Kvect_map[sig_m] += ((it).second); // ks = (it).second = number of time state s appear in the dataset
    }
    cout << endl;

// ***** Convert map to a vector:  for faster reading later on ********************************************
    vector<pair<uint64_t, unsigned int>> Kvect(Kvect_map.size());

    int i=0;
    for (auto& my_pair : Kvect_map)
    {
        Kvect[i]=my_pair;
        i++;
    }

    return Kvect;
}



