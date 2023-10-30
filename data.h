
/********************************************************************/
/**************************    PARAMETERS    ************************/
/********************************************************************/
// number of binary (spin) variables:
//unsigned int n = 50;

// INPUT DATA FILES (optional):  
// the input datafile can also be specified directly in the main() function, as an argument of the function "read_datafile()":
//std::string input_datafile = "INPUT/Big5PT.sorted_Ne5";  //"INPUT/Big5PT.sorted";  //"INPUT/SCOTUS_n9_N895_Data.dat"; //

// Cutoff criteria at alpha*sigma  (ex.  1*sigma or 3*sigma) for relevant basis operators
const unsigned int alpha = 3;

// OUTPUT FOLDER:
const std::string OUTPUT_directory = "OUTPUT/";

/********************************************************************/
/**************************    CONSTANTS    *************************/
/********************************************************************/
//const uint32_t NOp_tot = (uint32_t)((un << n) - 1);                     // number of operators = 2^n - 1

unsigned int bitset_count(__int128_t bool_nb);

struct Operator64
{
  uint64_t bin;     // binary representation of the operator
  unsigned int k1;  // nb of datapoints for which Op = 1 --> it's a R.V.:  k1 = sum(op[s^i])
  unsigned int r;   // nb of basis 

  double bias;     // bias = fabs(p1-0.5)
  bool operator < (const Operator64 &other) const   // for ranking Operators from the most to the less likely
    { return (bias > other.bias || (bias == other.bias && bitset_count(bin) < bitset_count(other.bin)) || (bias == other.bias && bitset_count(bin) == bitset_count(other.bin) && bin < other.bin)); }
};

struct Struct_LowerBound
{
  double Bias = 0.;
  unsigned int Index = 0.;
};


class MatrixF2 {
  public:
  unsigned int n; // = Number of spins = number of rows --> 1rst index  
  unsigned int m; // = Number of selected operators = number of columns --> 2nd index
  bool** M;       // = boolean Matrix

  unsigned int OpSet_offset=0; // number of Operators already analysed

  // constructor: 
  MatrixF2(unsigned int n_, unsigned int m_) {
    n = n_;
    m = m_;

    // Create a Boolean Matrix:
    M = (bool**) malloc(n*sizeof(bool*));  // n rows --> 1rst index
    for (int i=0; i<n; i++)
      {   M[i] = (bool*) malloc(m * sizeof(bool));  }  // m columns --> 2nd index
    }

  // destructor:
  ~MatrixF2() {
    //std::cout << "Object MatrixF2 is being deleted" << std::endl;
    for (int i=0; i<n; i++)
      {   free(M[i]);  }
    free(M);
  }
};

/*
class Data {   // ***** The data is stored in Nvect as an histogram:  
  unsigned int n;  // number of binary variables
  unsigned int N;  // total number of datapoints

  vector<pair<uint64_t, unsigned int>> Nvect; // Nvect[mu] = #of time state mu appears in the data set
};
*/





