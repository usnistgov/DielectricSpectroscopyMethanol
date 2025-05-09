Identification of chains and tracking chains through time

This process has two parts separated into two pieces of code: (1) identification of chains via id_chains.c and (2) tracking of chains through time via track_chain.c. Both can be compiled through the gcc compiler like so:
gcc id_chains.c -lm
gcc track_chains.c -lm
Either of these will yield an executable named a.out in the same directory which can then be run either with or without command line inputs. If command line inputs are not provided, the usage is:
./a.out
Then, the code will run and prompt the user for necessary inputs. Otherwise, command line inputs can be provided. This looks like:
./a.out input1 input2 input3
The inputs for the identification of chains are:
(1) number of molecules in the simulation (integer)
The inputs for the tracking of chains in time are:
(1) number of molecules in the simulation (integer)
The output of the identification of chains is a series of files:
chainmol# where the number is 1 through the final time step of the simulation
The output of the tracking of chains is a series of files:
chain#_start -> identifies the starting time step of chain #
chain#_end -> identiies the ending time step of chain #
chain#_#mol -> identifies the molecule IDs of molecules in chain # and time step #