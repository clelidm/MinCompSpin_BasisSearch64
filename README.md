# Search for the Best Basis Representation of a Binary Dataset
# or: Best Independent Minimally Complex Models (MCM)

This program searches for the **Best Basis representation for a chosen binary dataset**, while taking into account possible **high-order patterns of the data**. It was develop for the paper [*Statistical Inference of Minimally Complex Models*](https://arxiv.org/abs/2008.00520) [1]. More details on the general algorithm can be found in the paper.

[1]  C. de Mulatier, P. P. Mazza, M. Marsili, *Statistical Inference of Minimally Complex Models*, [arXiv:2008.00520](https://arxiv.org/abs/2008.00520)

## General information

The Best basis for a binary data with `n` variables is the one for which the independent model formed by `n` field operators has the largest log-likelihood (and therefore the largest log-evidence, as all independent model with the same number of operators are equivalent -- see Ref[1]).

There are three main functions that you can use to search for the best basis from the `main.cpp`:

 1) **Exhaustive Search:** This function will compute all $2^n-1$ operators and will search for the best basis among them with a Greedy approach (i.e. rank them from the most to least biased and extract the set of `n` most bias independent operators starting from the most biased one):
```c++
vector<Operator64> BestBasis_ExhaustiveSearch(vector<pair<uint64_t, unsigned int>> Nvect, unsigned int n, unsigned int N, bool bool_print = false)
```

 2) **In a fixed representation:** This function searches for the best Basis among all operators up to order `kmax` in the a given representation (the one in which the dataset stored in Nvect is written in):
```c++
BestBasisSearch_FixedRepresentation(vector<pair<uint64_t, unsigned int>> Nvect, unsigned int n, unsigned int N, unsigned int k_max, unsigned int B_it, bool bool_print = false)
```
If you take the largest order to be equal to the number of variables (`kmax = n`), then this function will perform an exhaustive search for the best basis among all possible operators. Note: we advise doing such search only for small systems (up to ~15 variables).

 3) **In varying representations:**
```c++
vector<Operator64> BestBasisSearch_Final(vector<pair<uint64_t, unsigned int>> Nvect, unsigned int n, unsigned int N, unsigned int k_max, bool bool_print = false)
```
This function performs the search procedure described in Ref.[1]. The codes searches for the best basis up to order `k_max`; the data is then successively transformed in the representation given by the previously found best basis and a new best basis searched for in this representation. The algorithm stops when the new basis found is identitity (i.e. the basis has not changed).

## Requirements

The code uses the C++11 version of C++.

## Usage with Makefile:

Run the following commands in your terminal, from the main folder (folder containing the `makefile` document):

 - **To compile:** `make`

 - **To Execute:** `make run` . This will use the datafile and variables that are specified in the makefile.

To change datafile: open the makefile and replace the values of the two following variables at the very top of the file (an example is provided):
>  - `datafile`: name of your datafile; this file must be placed in the folder `INPUT`
>  - `n`: number of variables in your file; maximum possible value `n = 128`.

You can also execute the code by running in your terminal the command (from the main folder):
```bash
./BestBasis.out  datafilename  n
```

where you must replace `datafilename` by the name of your datafile and `n` by your number of variables.

 - **To clean:** `make clean` (to use only once you're done using the code)

## Examples

All the functions that can be called from `int main()` are declared at the beginning of the `main.cpp` file. The most useful functions are described in the section "General information" above. 

In the `INPUT` folder, we provided two examples:
  - the binary dataset `SCOTUS_n9_N895_Data.dat`, which is the dataset of votes of the US Supreme court analysed in Ref.[3] and used as an example in Ref.[1]. 
  - the binary dataset `Big5-IPC1_VS3_Ne5.dat`: this is a binarized version of the first `100 000` samples of the Big 5 dataset [4], which has `50` variables. See paper [1] for comments on the Best Basis obtained for this dataset.

Each of these two datasets can be the one run as an example by commenting/uncommenting the correct datafile choice at the beginning of the `makefile`.

For hands-on and simple tests of the program, please check the examples in the function `int main()` of the `main.cpp` file. 

[3] E.D. Lee, C.P. Broedersz, W. Bialek, Statistical Mechanics of the US Supreme Court. [J Stat Phys 160, 275–301 (2015)](https://link.springer.com/article/10.1007/s10955-015-1253-6).

[4] Raw data from [Open-Source Psychometrics Project](https://openpsychometrics.org/_rawdata/) in the line indicated as "Answers to the IPIP Big Five Factor Markers"; [here](https://openpsychometrics.org/_rawdata/IPIP-FFM-data-8Nov2018.zip) is a direct link to the same zipfile.


## License

This code is an open source project under the GNU GPLv3.
