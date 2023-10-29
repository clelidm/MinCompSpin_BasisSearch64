########################################################################################################################
####################################      CONSTANT TO SPECIFY     ######################################################
########################################################################################################################
#### If the following two variables are left empty:           ##########################################################
####     then the program will automatically use the filename and number of variables specifed in 'main.cpp'   #########
########################################################################################################################

#### EXAMPLE: Big 5 dataset:
datafilename := Big5-IPC1_VS3_Ne5.dat  # datafile name ### IMPORTANT: this file must be in the 'INPUT' folder
n := 50		# number of binary variables in the datafile 

#### EXAMPLE 2: US Supreme Court dataset:
#datafilename := SCOTUS_n9_N895_Data.dat  
#n := 9		

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

objectsWdataH = BestBasis_Init_SetOp.o BestBasis.o BasisTools.o BestBasisSearch.o
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

BestBasis_Init_SetOp.o: BestBasis_Init_SetOp.cpp data.h
	g++ $(CXXFLAGS) -c BestBasis_Init_SetOp.cpp -o BestBasis_Init_SetOp.o   # Compile main.cpp

BestBasis.o: BestBasis.cpp data.h
	g++ $(CXXFLAGS) -c BestBasis.cpp -o BestBasis.o   # Compile main.cpp

BasisTools.o: BasisTools.cpp data.h
	g++ $(CXXFLAGS) -c BasisTools.cpp -o BasisTools.o   # Compile main.cpp

BestBasisSearch.o: BestBasisSearch.cpp data.h
	g++ $(CXXFLAGS) -c BestBasisSearch.cpp -o BestBasisSearch.o   # Compile main.cpp

clean:
	rm main.o $(objects) $(objectsWdataH) BestBasis.out #main

run:
	time ./BestBasis.out $(datafilename) $n
	#time ./main $(datafilename) $n


