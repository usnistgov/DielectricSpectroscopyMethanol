Identification of hydrogen bonds

The identification of hydrogen bonds in the simulation involves an algorithm specified in the accompanying publication. Notably, a threshold length and angle are both required beyond which a particular configuration will not be considered. The length is the maximal allowed distance between the two atoms being considered beyond which their interactions should not be considered a hydrogen bond. The angle ensures that the orientation of the two atoms is actually toward each other (i.e., two hydroxyl groups may be near enough to look like they are hydrogen bonding but have the hydrogens of the two groups pointed at each other, which is not a hydrogen bond). This code uses both of these criteria to identify hydrogen bonds in the simulation from the positions of specified atom types.

The code to undertake this identification is find_hbonds.c. It can be compiled like so:
gcc find_hbonds.c -lm
This creates an executable a.out in the folder in which the code was compiled. Then, the code can be run with or without command line inputs. Without command line inputs, the code is run like so:
./a.out
With command line inputs, the code is run like so:
./a.out input1 input2 input3
If the code is run without command line inputs, the code will instead prompt the user to input these values. The expected inputs of this code are:
(1) simbox.txyz -> name of the tinker input configuration file
(2) maxdistance -> maximal distance in Angstrom allowed for a hydrogen bond
The code will use the name of the tinker input configuration file to construct the expected name of the tinker output configuration file, the tinker output static polarization file, and the tinker output induced polarization file:
simbox.arc -> tinker output configuration file
simbox.ustc -> tinker output atomic static polarization file
simbox.uind -> tinker output atomic induced polarization file
If these files are not present in the folder in which the code is being executed, the code will not proceed.

In addition to these inputs, the code expects a file named hb_types in the same directory as the execution of the code. This file should contain two lines:
(1) type1
(2) type2
The types specified here are the tinker types of the atoms through which hydrogen bonding occurs. For example, methanol hydrogen bonds via the hydroxyl hydrogen of one molecule to the hydroxyl oxygen of another molecule. If the hydroxyl hydrogen has tinker type 1 and the hydroxyl oxygen has tinker type 2, the contents of this file should be
1
2
This will ensure that only the atoms which are capable of hydrogen bonding are considered for this identification algorithm.