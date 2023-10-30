########################################################################################################################
####################################      CONSTANT TO SPECIFY     ######################################################
########################################################################################################################
#### If the following two variables are left empty:           ##########################################################
####     then the program will automatically use the filename and number of variables specifed in 'main.cpp'   #########
########################################################################################################################

#### EXAMPLE: Big 5 dataset:
#datafilename := Big5-IPC1_VS3_Ne5.dat  # datafile name ### IMPORTANT: this file must be in the 'INPUT' folder
#n := 50		# number of binary variables in the datafile 

#### EXAMPLE 2: US Supreme Court dataset:
datafilename := SCOTUS_n9_N895_Data.dat  
n := 9		

########################################################################################################################
######## ENTER THE FOLLOWING IN YOUR TERMINAL:
#### TO COMPILE:  	make
#### TO RUN: 		make run
#### TO CLEAN:  	make clean    --> to use only when you are completely done
########################################################################################################################


########################################################################################################################
##########################################      DO NOT MODIFY     ######################################################
########################################################################################################################
CC = g++ 	# Flag for implicit rules: used for linker
CXX = g++ 	# Flag for implicit rules: compilation of c++ files
CXXFLAGS = -std=c++11 #-O2  #Extra flags to give to the C++ compiler

objects = tools.o ReadDataFile.o #Basis_Choice.o 

objectsWdataH = Init_OpSet.o ExtractBasis_inOpSet.o BasisTools.o BestBasis_IterativeSearch.o BestBasis_ExhaustiveSearch.o
#directory =  Folder

# Compilation -- Implicite rule:
#%.o : %.c   
#		$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $^ -o $@

# Link -- Implicite rule:
#$(CC) $(LDFLAGS) $^ $(LOADLIBES) $(LDLIBS) -o $@

#main: $(objects) $(objectsWdataH) main.o #main.o tools.o  # implicit link #g++ main.o tools.o -o main

BestBasis.out: $(objects) $(objectsWdataH) main.o 
	g++ $(CXXFLAGS) main.o $(objects) $(objectsWdataH) -o BestBasis.out

main.o: main.cpp data.h
	g++ $(CXXFLAGS) -c main.cpp -o main.o   # Compile main.cpp

Init_OpSet.o: Init_OpSet.cpp data.h
	g++ $(CXXFLAGS) -c Init_OpSet.cpp -o Init_OpSet.o   # Compile main.cpp

ExtractBasis_inOpSet.o: ExtractBasis_inOpSet.cpp data.h
	g++ $(CXXFLAGS) -c ExtractBasis_inOpSet.cpp -o ExtractBasis_inOpSet.o   # Compile main.cpp

BasisTools.o: BasisTools.cpp data.h
	g++ $(CXXFLAGS) -c BasisTools.cpp -o BasisTools.o   # Compile main.cpp

BestBasis_IterativeSearch.o: BestBasis_IterativeSearch.cpp data.h
	g++ $(CXXFLAGS) -c BestBasis_IterativeSearch.cpp -o BestBasis_IterativeSearch.o   # Compile main.cpp

BestBasis_ExhaustiveSearch.o: BestBasis_ExhaustiveSearch.cpp data.h
	g++ $(CXXFLAGS) -c BestBasis_ExhaustiveSearch.cpp -o BestBasis_ExhaustiveSearch.o   # Compile main.cpp

clean:
	rm main.o $(objects) $(objectsWdataH) BestBasis.out #main

run:
	time ./BestBasis.out $(datafilename) $n
	#time ./main $(datafilename) $n


