# MinCompSpin_BasisSearch
This codes looks for the Best Basis representation for binary datasets

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
