#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

/*

    This code identifies hydrogen bonds from the configurations
    of a molecular dynamics simulation by considering two types 
    of atoms within the simulation in terms of their distance 
    apart and orientation toward each other. If both criteria 
    are met, the two are considered to be part of a hydrogen
    bond.

    Input file:    hb_types
                        Line 1: type1
                        Line 2: type2
    Command line or requested inputs:
        simbox.txyz -> tinker input configuration file
        maxdist     -> maximal allowed distance for consideration


    Output files: (for each time step # between 1 and total time steps)
        hbmol# -> molecule IDs of molecules that are HB
        hbatom# -> atom IDs through which molecules are HB
        hb_indpol# -> induced polarization vector projected onto HB
        hb_stapol# -> static polarization vector projected onto HB
        hb_length# -> length of hydrogen bonds

*/


// Subroutines -------------------------------------------------------
// -------------------------------------------------------------------

// File existence check
int chk_file(char name[]) {
    // Function takes string array containing file name and locations
    // Function outputs flag signalling if indicated file exists
    int pacc = 0;
    if ( access( name,F_OK ) == -1) {
        pacc = 0; // No File found
    } else {
        pacc = 1; // File found
    }
    return pacc;
}

// Check for exit request
int check_exit(char name[]) {
    int flag = 0; // Flag for exit
    if (strcmp(name,"exit") == 0) {
        flag = -1; // Exit requested
    } else if (strcmp(name,"Exit") == 0) {
	flag = -1; // Exit requested
    } else if (strcmp(name,"EXIT") == 0) {
	flag = -1; // Exit requeste
    }
    return flag;
}

// Main function -----------------------------------------------------
// -------------------------------------------------------------------
void main(int argc,char* argv[]) {

    // Print header ---------------------------------------------------

    // Initial declarations  ------------------------------------------
 
    int x,t,i,j,k,comp,molid,at,ind_close; // Integer variables for loops and flags
    int natom,nmol; // Integer variables for size of arrays
    int type1,type2; // Integer variables of tinker types of HBing atoms
    int nhb; // Integer counter of hydrogen bonds found
    double dist1, dist2, dist3, dist, maxdist,angle; // Doubles for distance calculation
    double proj_ind1,proj_ind2,proj_sta1,proj_sta2; // Final projection onto HB axis
    double veclen; // Vector length
    double xmin,xmax,ymin,ymax,zmin,zmax,dx,dy,dz; // Doubles for boundaries
    char fname[100]; // Character array for file names
    char fline[1200]; // Character array for file reading
    char delim[] = " "; // Expected file delimiter
    char fxyz[100]; // Character array for name of tinker xyz file
    char farc[100]; // Character array for name of tinker arc file
    char fustc[100]; // Character array for name of tinker ustc file
    char fuind[100]; // Character array for name of tinker uind file
    double PI = 3.14159265358979323; // constant

    // Process input ---------------------------------------------------

    // Read or request input parameters
    if (argc < 2) { // No command line inputs
	printf("No input parameters detected.\n"); // Report to terminal

	// Tinker initial configuration (xyz) file
	printf("Tinker xyz file: "); // Request initial configuration file
	scanf("%s",&fxyz); // User input from command line
	x = chk_file(fxyz); // Validate that input file exists
	while (x == 0) { // If the file doesn't exist, keep asking
	    x = check_exit(fxyz); // Check for exit request
	    if (x == -1) { // Hard error
		exit(EXIT_FAILURE);
	    }
	    printf("Tinker xyz file: "); // Keep asking
	    scanf("%s",&fxyz); // User input from command line
	    x = chk_file(fxyz); // Validte user input
	}

	// Allowed length of hydrogen bond
	printf("Maximum length of hydrogen bond (Angstrom): "); //Request max length
	scanf("%lf",&maxdist); // User input from command line
        while (maxdist <= 0) {
	    printf("Invalid maximum length of hydrogen bond.\n"); // Warning
	    printf("Maximum length of hydrogen bond (Angstrom): "); // Keep asking
	    scanf("%lf",&maxdist);
	}
    } else if (argc < 3) { // Too few command line inputs
	// Assume nohup usage and exit
	printf("ERROR:\n Too few command line arguments detected.\n"); // Hard error
	exit(EXIT_FAILURE);
    } else if (argc == 3) { // Correct number of command line inputs
	
	// Tinker xyz file
	sprintf(fxyz,"%s",argv[1]); // Assign first input
        x = chk_file(fxyz); // Validate input
        if (x == 0) { // Input file does not exist
	    // Assume nohup usage and exit
	    printf("ERROR:\n Specified tinker xyz file does not exist.\n"); // Hard error
	    exit(EXIT_FAILURE);
	}

        // Maximum length of hydrogen bond
	maxdist = strtod(argv[2],NULL); // Assign third input
	if (maxdist <= 0) { // Validate input
	    // Assume nohup usage and exit
	    printf("ERROR:\n Specified maximal hydrogen bond length is invalid.\n"); // Hard error
	    exit(EXIT_FAILURE);
	}
    } else { // Too many command line inputs
	// Assume nohup usage and exit
	printf("ERROR:\n Too many command line inputs detected.\n"); // Hard error
	exit(EXIT_FAILURE);
    }

    //printf("Parsing file names.\n"); // Report to terminal -> DEBUG

    // Parse xyz file name for arc, uind, ustc  file names and validate
    sprintf(fname,"%s",fxyz); // Print file name to character string for parsing
    char *ptr = strtok(fname,"."); // Parsing name up to first period
    sprintf(farc,"%s.arc",ptr); // Print arc file name to character string
    x = chk_file(farc); // Check for file existence
    if (x == 0) { // Hard error
	printf("ERROR:\n No arc file found.\n");
	printf("All files must be in the working directory.\n");
	exit(EXIT_FAILURE);
    }
    sprintf(fuind,"%s.uind",ptr); // Print uind file name to character string
    x = chk_file(fuind);
    if (x == 0) { // Hard error
	printf("ERROR:\n No uind file found.\n");
	printf("All files must be in the working directory.\n");
	exit(EXIT_FAILURE);
    }
    sprintf(fustc,"%s.ustc",ptr); // Print ustc file name to character string
    x = chk_file(fustc);
    if (x == 0) { // Hard error
	printf("ERROR:\n No usc file found.\n");
        printf("All files must be in the working directory.\n");
	exit(EXIT_FAILURE);
    }

    // Read in tinker types of Hydrogen Bonding atoms (e.g., O and H)
    sprintf(fname,"%s","hb_types"); // Print file name to character string
    x = chk_file(fname); // Check for file existence
    if (x == 0) { // Hard error
	printf("ERROR:\n No hydrogen bonding tinker types file detected.\n");
	printf("All files must be in the working directory.\n");
        exit(EXIT_FAILURE);
    }
    FILE * hbi; // Pointer to file
    hbi = fopen(fname,"r"); // Open file for reading
    fscanf(hbi,"%d\n %d\n",&type1,&type2); // Read indeces
    fclose(hbi); // Close file

    // Initialization -------------------------------------------------

    printf("Beginning initialization.\n"); // Report to terminal

    // Open tinker xyz file
    FILE * xyz; // Pointer to file
    xyz = fopen(fxyz,"r"); // Open file for reading
    fgets(fline,1200,xyz); // Retrieve lin
    ptr = strtok(fline,delim); // Parse first line entry
    natom = atoi(ptr); // Number of atoms

    printf("%d atoms in simulation\n",natom); // Report to terminal

    // Allocate array for tinker types
    int* ttype;
    ttype = (int*) malloc(natom*sizeof(int));
    for (i = 0; i < natom; i++) {
	ttype[i] = 0;
    }
    //printf("%d\n",sizeof(ttype)/sizeof(ttype[0]));

    // Allocate array for Molecule ID
    int *mol;
    mol = (int*) malloc(natom*sizeof(int));
    for (i = 0; i < natom; i++) {
	mol[i] = 0;
    }
    //printf("%d\n",mol[0]);
    //printf("%d\n",sizeof(mol)/sizeof(mol[0]));

    printf("Determining number of molecules.\n");

    // Read tinker types and determine molecule id from bonds
    nmol = 0; // Zero out number of molecules
    for (i = 0; i < natom; i++) {
	//printf("i=%d\n",i);
	//printf("nmol=%d\n",nmol);
	if (mol[i] == 0) { // if molecule ID is not yet assigned
	    nmol++; // Increment number of molecules
	    mol[i] = nmol; // Assign new molecule number
	}
	fgets(fline,1200,xyz); // Retrieve line
	ptr = strtok(fline,delim); // Atomic number
	ptr = strtok(NULL,delim); // Atomic symbol
	ptr = strtok(NULL,delim); // x(i)
	ptr = strtok(NULL,delim); // y(i)
	ptr = strtok(NULL,delim); // z(i)
	ptr = strtok(NULL,delim); // tinker id
        ttype[i] = atoi(ptr);
	ptr = strtok(NULL,delim); // Possible bound atom
	while (ptr != NULL) { // If end of line has net yet been reached
	    j = atoi(ptr); // Parse atom number
	    if (j != 0) {
	        //if (mol[j-1] != 0) {
                //    nmol = mol[j-1];
		//    mol[i] = mol[j-1]; // Assign molecule id
	        //} else {
	        mol[j-1] = nmol; // Assign molecule ID
	        //}
	        ptr = strtok(NULL,delim); // Possible bound atom
            } else {
                ptr = NULL;
	    }
        }
    }
    printf("%d molecules detected.\n",nmol);
    fclose(xyz); // Close xyz file

    printf("Maximum length of HB: %lf Angstrom\n",maxdist);

    // Allocate array for position (natom x 3 array)
    double **pos; // Pointer to array
    pos = (double**) malloc(natom*sizeof(double*)); // Allocate first level
    for (i = 0; i < natom; i++) { // Loop through first level
	pos[i] = (double*) malloc(3*sizeof(double)); // Allocate second level
    }
    // Zero out allocated array
    for (i = 0; i < natom; i++) {
	for (j = 0; j < 3; j++) {
	    pos[i][j] = 0;
	}
    }

    // Allocate array for static polarization (natom x 3 array)
    double **stapol; // Pointer to array
    stapol = (double**) malloc(natom*sizeof(double*)); // Allocate first leel
    for (i = 0; i < natom; i++) { // Loop through first level
	stapol[i] = (double*) malloc(3*sizeof(double)); // Allocate second level
    }
    // Zero out allocated array
    for (i = 0; i < natom; i++) {
	for (j = 0; j < 3; j++) {
	    stapol[i][j] = 0;
	}
    }

    // Allocate array for induced polarization (natom x 3 array)
    double **indpol; // Pointer to array
    indpol = (double**) malloc(natom*sizeof(double*)); // Allocate first level
    for (i = 0; i < natom; i++) { // Loop through first level
	indpol[i] = (double*) malloc(3*sizeof(double)); // Allocate second level
    }
    // Zero out allocated array
    for (i = 0; i < natom; i++) {
	for (j = 0; j < 3; j++) {
	    indpol[i][j] = 0;
	}
    }

    // Open arc file for reading
    FILE * arc; // Pointer to file
    arc = fopen(farc,"r"); // Open file

    // Open ustc file for reading
    FILE * ustc; // Pointer to file
    ustc = fopen(fustc,"r"); // open file

    // Open uind file for reading
    FILE * uind; // Pointer to file
    uind = fopen(fuind,"r"); // Open file

    // Open file for number of hydrogen bonds
    sprintf(fname,"count_hb"); // name of file to character string
    FILE * ct_hb; // Pointer to file
    ct_hb = fopen(fname,"w"); // Open file for writing

    // Declare vector for polarization projections
    double axis[3];
    axis[0] = 0;
    axis[1] = 0;
    axis[2] = 0;

    // Declare vector for relative location of closest HB-type atom
    double relpos[3];
    relpos[0] = 0;
    relpos[1] = 0;
    relpos[2] = 0;

    // Declare vector 1 for relative angle of HB
    double vec1[3];
    vec1[0] = 0;
    vec1[1] = 0;
    vec1[2] = 0;

    // Declare vector 2 for relative angle of HB
    double vec2[3];
    vec2[0] = 0;
    vec2[1] = 0;
    vec2[2] = 0;

    printf("Initialization complete.\n"); // Report to terminal

    // Analysis -------------------------------------------------------
    
    printf("Beginning analysis.\n"); // Report to terminal

    t = 0; // Time frame number
    comp = 0; // Flag for end of file
    while ((comp == 0) && (t < 500000000)) {
	
	// Update time
	t++; // Increment time frame number
        //printf("t=%d\n",t); // Report to terminal -> DEBUG

        // Positions ----------------------------------

	// Find first line of current frame of atomic positions
	if (fgets(fline,1200,arc) == NULL) { // Check for end of file
	    comp = 1; // End of file flag
	    break; // Break from analysis loop
	}
	ptr = strtok(fline,delim); // Parse first line entry
	sprintf(fname,"%s",ptr); // First entry to character string
	while (strcmp(fname,"1") != 0) { // Until first atom index found
	    if (fgets(fline,1200,arc) == NULL) { // Check for end of file
	        comp = 1; // End of file flag
		break; // Break rom current loop
	    }
	    ptr = strtok(fline,delim); // Parse first entry
	    sprintf(fname,"%s",ptr); // First entry to character string
	}
	if (comp == 1) { // If end of file flag active
            break; // Break from analysis loop
	}

        // Read atom positions
	for (i = 0; i < natom; i++) {
	    if (i != 0) { // On subsequent atoms, retrieve next line
		if (fgets(fline,1200,arc)  == NULL) { // Check for end of file
		    comp = 1; // End of file flag
		    break; // Break from current loop
		}
		ptr = strtok(fline,delim); // Atom number
	    }
	    ptr = strtok(NULL,delim); // Atomic symbol
	    ptr = strtok(NULL,delim); // x(i)
            pos[i][0] = strtod(ptr,NULL); // Assign to array position
	    ptr = strtok(NULL,delim); // y(i)
	    pos[i][1] = strtod(ptr,NULL); // Assign to array position
	    ptr = strtok(NULL,delim); // z(i)
	    pos[i][2] = strtod(ptr,NULL); // Assign to array position
	}

        //printf("Atomic positions read.\n"); // Report to terminal -> DEBUG

	// Determine axis extents
	xmin = 10000000; // Set arbitrary start
	xmax = -10000000; // Set arbitrary start
	ymin = 10000000; // Set arbitrary start
	ymax = -10000000; // Set arbitrary start
	zmin = 10000000; // Set arbitrary start
	zmax = -10000000; // Set arbitrary start
	for (i = 0; i < natom; i++) {
	    if (pos[i][0] < xmin) { // Test xmin
		xmin = pos[i][0]; // Replace value
	    }
	    if (pos[i][0] > xmax) { // Test max
		xmax = pos[i][0]; // Replace value
	    }
	    if (pos[i][1] < ymin) { // Test ymin
		ymin = pos[i][1]; // Replace value
	    }
	    if (pos[i][1] > ymax) { // Test ymax
		ymax = pos[i][1]; // Replace value
	    }
	    if (pos[i][2] < zmin) { // Test zmin
		zmin = pos[i][2]; // Replace value
	    }
	    if (pos[i][2] > zmax) { // Test zmax
		zmax = pos[i][2]; // Replace value
	    }
	}
	dx = xmax - xmin; // x axis extent
	dy = ymax - ymin; // y axis extent
	dz = zmax - zmin; // z axis extent

        // Move atoms in molecules such that they are all on one side of all pbc
	molid = 0; // Current molecule id (-1)
	at = 0; // current starting atom for current molecule
	for (i = 0; i < natom; i++) { // Loop through each atom in simulation
	    if ((mol[i] == molid+1) && (i != at)) { // Check for if current atom is first atom of new molecule
		// x axis adjustment
		dist1 = fabs(pos[at][0] - pos[i][0]); // No x adjustment
		dist2 = fabs(pos[at][0] - pos[i][0] - dx); // + dx
		dist3 = fabs(pos[at][0] - pos[i][0] + dx); // - dx
		if ((dist1 < dist2) && (dist1 < dist3)) {
		    // Current x position is not across boundary
		} else if ((dist2 < dist1) && (dist2 < dist3)) {
		    // Add dx
		    pos[i][0]+=dx;
		} else {
		    // Subtract dx
		    pos[i][0]-=dx;
		}
		// y axis adjustment
		dist1 = fabs(pos[at][1] - pos[i][1]); // No y adjustment
		dist2 = fabs(pos[at][1] - pos[i][1] - dy); // + dy
		dist3 = fabs(pos[at][1] - pos[i][1] + dy); // - dy
		if ((dist1 < dist2) && (dist1 < dist3)) {
		    // Current y position is not across boundary
		} else if ((dist2 < dist1) && (dist2 < dist3)) {
		    // Add dy
		    pos[i][1]+=dy;
		} else {
		    // Subtract dy
		    pos[i][1]-=dy;
		}
		//z axis adjustment
		dist1 = fabs(pos[at][2] - pos[i][2]); // No z adjustment
		dist2 = fabs(pos[at][2] - pos[i][2] - dz); // + dz
		dist3 = fabs(pos[at][2] - pos[i][2] + dz); // - dz
		if ((dist1 < dist2) && (dist1 < dist3)) {
		    // Current z position is not across boundary
		} else if ((dist2 < dist1) && (dist2 < dist3)) {
		    // Add dz
		    pos[i][2]+=dz;
		} else {
		    // Subtract dz
		    pos[i][2]-=dz;
		}
	    }
	}

        // Induced polarization -----------------------------------

        // Find beginning of next time frame
	if (fgets(fline,1200,uind)  == NULL) { // Check for end of file
	    comp = 1; // End of file flag
	    break; // Break from analysis loop
	}
	if (fgets(fline,1200,uind)  == NULL) { // Check for end of file
	    comp = 1; // End of file flag
	    break; // Break from analysis loop
	}
	ptr = strtok(fline,delim); // Parse first line entry
	sprintf(fname,"%s",ptr); // Print to character string
	x = strcmp(fname,"1"); // Check if first entry is "1"
	while (x != 0) {
	    if (fgets(fline,1200,uind)  == NULL) { // Check for end of file
                comp = 1; // End of file flag
		x = 0; // Flag to break current loop
		break; // Break from current loop
	    }
	    ptr = strtok(fline,delim); // Parse first line entry
	    sprintf(fname,"%s",ptr);  // Print to character string
	    x = strcmp(fname,"1"); // Check if first entry is "1"
	}

	// Read per atom values from file
	for (i = 0; i < natom; i++) { // Loop through atoms in simulation
	    if (i != 0) { // If current atom is not the first
		if (fgets(fline,1200,uind) == NULL) { // Check for end of file
		    comp = 1; // End of file flag
		    break; // Break from current loop
		}
		ptr = strtok(fline,delim); // Parse atom number
	    }
	    ptr = strtok(NULL,delim); // Atomic syhmbol
	    ptr = strtok(NULL,delim); // uind(x)
	    indpol[i][0] = strtod(ptr,NULL); // Record
	    ptr = strtok(NULL,delim); // uind(y)
	    indpol[i][1] = strtod(ptr,NULL); // Record
	    ptr = strtok(NULL,delim); // uind(z)
	    indpol[i][2] = strtod(ptr,NULL); // Record
	}
	if (comp == 1) { // If end of file flag is active
	    break; // Break from analysis loop
	}

        // Static polarization -----------------------------------
	
	// Find beginning of next time frame
	if (fgets(fline,1200,ustc) == NULL) { // Check for end of file
	    comp = 1; // End of file flag
	    break; // Break from analysis loop
	}
	if (fgets(fline,1200,ustc) == NULL) { // Check for end of file
	    comp = 1; // End of file flag
	    break; // Break from analysis loop
	}
	ptr = strtok(fline,delim); // Parse first entry
	sprintf(fname,"%s",ptr); // Print to character string
	x = strcmp(fname,"1"); // Check if first entry is "1"
	while (x != 0) {
	    if (fgets(fline,1200,ustc) == NULL) { // Check for end of file
	        comp = 1; // End of file flag
		break; // Break from current loop
	    }
	    ptr = strtok(fline,delim); // Parse first entry
	    sprintf(fname,"%s",ptr); // Print to character string
	    x = strcmp(fname,"1"); // Check if first entry is "1"
	}
	if (comp == 1) { // Check for active end of file flag
	    break; // Break from analysis loop
	}

	// Read per atom values from file
	for (i = 0; i < natom; i++) { // Loop through atoms in simulation
            if (i != 0) { // On subsequent, read next line
		if (fgets(fline,1200,ustc) == NULL) {
                    comp = 1; // End of file flag
		    break; // Break from current loop
		}
		ptr = strtok(fline,delim); // Parse atom nunmber
	    }
	    ptr = strtok(NULL,delim); // Atomic symbol
	    ptr = strtok(NULL,delim); // ustc(x)
	    stapol[i][0] = strtod(ptr,NULL); // Record
	    ptr = strtok(NULL,delim); // ustc(y)
	    stapol[i][1] = strtod(ptr,NULL); // Record
	    ptr = strtok(NULL,delim); // ustc(z)
	    stapol[i][2] = strtod(ptr,NULL); // Record
	}
	if (comp == 1) { // If end of file flag is active
	    break; // Break from analysis loop
	}

        // Hydrogen bond identification --------------------------

        // Open file for HB molecules at current step
	sprintf(fname,"hbmol%d",t); // File name to character string
	FILE * currhb_mol; // Pointer to file
	currhb_mol = fopen(fname,"w"); // Open file

	// Open file for HB atoms at current step
	sprintf(fname,"hbatom%d",t); // File name to character string
	FILE * currhb_atom; // Pointer to file
	currhb_atom = fopen(fname,"w"); // Open file

        // Open file for HB static polarization at current step
	sprintf(fname,"hb_stapol%d",t); // File name to character string
	FILE * currsta; // Pointer to file
	currsta = fopen(fname,"w"); // Open file
	
	// Open file for HB induced polarization at current step
        sprintf(fname,"hb_indpol%d",t); // File name to character string
        FILE * currind; // Pointer to file
	currind = fopen(fname,"w"); // Open file

        // Open file for length of HB at current step
	sprintf(fname,"hb_length%d",t); // File name to character string
	FILE * currlen; // Pointer to file
        currlen = fopen(fname,"w"); // Open file

        // Search for hydrogen bonds
        for (i = 0; i < natom-1; i++) { // Loop 1 to N-1
            if (ttype[i] == type1) { // Compare current atom type to HB type -> Oxygen
                for (j = i+1; j < natom; j++) { // Loop i+1 to N
		    if (mol[i] != mol[j]) { // If not the same molecule
	                if (ttype[j] == type2) { // Compare second atom to other HB type -> Hydrogen
		            // Check x position vs boundary
			    dist1 = fabs(pos[i][0] - pos[j][0]); // No x adjust
			    dist2 = fabs(pos[i][0] - pos[j][0] - dx); // + dx
			    dist3 = fabs(pos[i][0] - pos[j][0] + dx); // - dx
			    if ((dist1 < dist2) && (dist1 < dist3)) {
			        dist = pow(dist1,2); // No adjust
			    } else if ((dist2 < dist1) && (dist2 < dist3)) {
				dist = pow(dist2,2); // + dx
			    } else {
				dist = pow(dist3,2); // - dx
			    }
			    // Check y position vs boundary
			    dist1 = fabs(pos[i][1] - pos[j][1]); // No y adjust
			    dist2 = fabs(pos[i][1] - pos[j][1] - dy); // + dy
			    dist3 = fabs(pos[i][1] - pos[j][1] + dy); // - dy
			    if ((dist1 < dist2) && (dist1 < dist3)) {
				dist+=pow(dist1,2); // No adjust
		            } else if ((dist2 < dist1) && (dist2 < dist3)) {
				dist+=pow(dist2,2); // + dy
			    } else {
				dist+=pow(dist3,2); // - dy
			    }
			    // Check z position vs boundary
			    dist1 = fabs(pos[i][2] - pos[j][2]); // No z adjust
			    dist2 = fabs(pos[i][2] - pos[j][2] - dz); // + dz
			    dist3 = fabs(pos[i][2] - pos[j][2] + dz); // - dz
			    if ((dist1 < dist2) && (dist1 < dist3)) {
				dist+=pow(dist1,2); // No adjust
			    } else if ((dist2 < dist1) && (dist2 < dist3)) {
				dist+=pow(dist2,2); // + dz
			    } else {
				dist+=pow(dist3,2); // - dz
			    }
			    // Determine distance
			    dist = sqrt(dist);
			    
                            // Compare to maximal allowed distance
                            if (dist < maxdist) {
                                
                                // Relative position of second atom-> x coordinate
		                dist1 = fabs(pos[i][0] - pos[j][0]); // No x adjust
		                dist2 = fabs(pos[i][0] - pos[j][0] - dx); // + dx
		                dist3 = fabs(pos[i][0] - pos[j][0] + dx); // - dx
		                if ((dist1 < dist2) && (dist1 < dist3)) {
			            relpos[0] = pos[j][0]; // No x adjust
		                } else if ((dist2 < dist1) && (dist2 < dist3)) {
			            relpos[0] = pos[j][0] + dx; // + dx
		                } else {
			            relpos[0] = pos[j][0] - dx; // - dx
		                }
		                // Relative position of second atom-> y coordinate
		                dist1 = fabs(pos[i][1] - pos[j][1]); // No y adjust
		                dist2 = fabs(pos[i][1] - pos[j][1] - dy); // + dy
		                dist3 = fabs(pos[i][1] - pos[j][1] + dy); // - dy
		                if ((dist1 < dist2) && (dist1 < dist3)) {
			            relpos[1] = pos[j][1]; // No y adjust
		                } else if ((dist2 < dist1) && (dist2 < dist3)) {
			            relpos[1] = pos[j][1] +dy; // + dy
		                } else {
			            relpos[1] = pos[j][1]-dy; // - dy
		                }
		                // Relative position of second atom-> z coordinate
		                dist1 = fabs(pos[i][2] - pos[j][2]); // No z adjust
		                dist2 = fabs(pos[i][2] - pos[j][2] - dz); // + dz
		                dist3 = fabs(pos[i][2] - pos[j][2] + dz); // - dz
		                if ((dist1 < dist2) && (dist1 < dist3)) {
			            relpos[2] = pos[j][2]; // No z adjust
		                } else if ((dist2 < dist1) && (dist2 < dist3)) {
			            relpos[2] = pos[j][2] + dz; // + dz
		                } else {
			            relpos[2] = pos[j][2] - dz; // - dz
		                }

                                // Check angle between vectors

                                // FInd oxygen of molecule 2 for angle comparison
                                angle = 0;
                                for (k = 0; k < natom; k++) {
                                    if (mol[j] == mol[k]) {
                                        if (ttype[k] == type1) {
                                            // Determine O(1)...H(2)-O(2) angle
                                            // Vector 1: pos(i) - pos(j)
                                            vec1[0] = pos[i][0]-relpos[0]; // x component
                                            vec1[1] = pos[i][1]-relpos[1]; // y coponent
                                            vec1[2] = pos[i][2]-relpos[2]; // z component

                                            // Vector 2: pos(k) - pos(j)
                                            // x component
                                            dist1 = fabs(pos[k][0] - pos[j][0]); // No adjustment to x
                                            dist2 = fabs(pos[k][0] - pos[j][0] - dx); // + dx
                                            dist3 = fabs(pos[k][0] - pos[j][0] + dx); // - dx
                                            if ((dist1 < dist2) && (dist1 < dist3)) {
                                                vec2[0] = pos[k][0] - pos[j][0]; // no adjustment to x
                                            } else if ((dist2 < dist1) && (dist2 < dist3)) {
                                                vec2[0] = pos[k][0] - pos[j][0] - dx; // + dx
                                            } else { 
                                                vec2[0] = pos[k][0] - pos[j][0] + dx; // - dx
                                            }
                                            // y component
                                            dist1 = fabs(pos[k][1] - pos[j][1]); // No adjustment to y
                                            dist2 = fabs(pos[k][1] - pos[j][1] - dy); // + dy
                                            dist3 = fabs(pos[k][1] - pos[j][1] + dy); // - dy
                                            if ((dist1 < dist2) && (dist1 < dist3)) {
                                                vec2[1] = pos[k][1] - pos[j][1]; // No adjustment to y
                                            } else if ((dist2 < dist1) && (dist2 < dist3)) {
                                                vec2[1] = pos[k][1] - pos[j][1] - dy; // + dy
                                            } else {
                                                vec2[1] = pos[k][1] - pos[j][1] + dy; // - dy
                                            }
                                            // z component
                                            dist1 = fabs(pos[k][2] - pos[j][2]); // No adjustment to z
                                            dist2 = fabs(pos[k][2] - pos[j][2] - dz); // + dz
                                            dist3 = fabs(pos[k][2] - pos[j][2] + dz); // - dz
                                            if ((dist1 < dist2) && (dist1 < dist3)) {
                                                vec2[2] = pos[k][2] - pos[j][2]; // No adjustment to z
                                            } else if ((dist2 < dist1) && (dist2 < dist3)) {
                                                vec2[2] = pos[k][2] - pos[j][2] - dz; // + dz
                                            } else {
                                                vec2[2] = pos[k][2] - pos[j][2] + dz; // - dz
                                            }
                                            angle = 1;
                                            break; // Break from loop
                                        }
                                    }
                                }
                                if (angle == 1) { // Check for flag that atom was found
                                    // Determine angle O(1,i)-H(2,j)-O(2,k)
                                    angle = vec1[0] * vec2[0]; // v1(1)*v2(1)
                                    angle+=(vec1[1]*vec2[1]); // +v1(2)*v2(2)
                                    angle+=(vec1[2]*vec2[2]); // +v1(3)*v2(3)
                                    // magnitude of vector 1
                                    veclen = pow(vec1[0],2); // x component
                                    veclen+=pow(vec1[1],2); // y component
                                    veclen+=pow(vec1[2],2); // z component
                                    angle/=sqrt(veclen);
                                    // magnitude of vector 2
                                    veclen = pow(vec2[0],2); // x component
                                    veclen+=pow(vec2[1],2); // y component
                                    veclen+=pow(vec2[2],2); // z component
                                    angle/=sqrt(veclen);
                                    // acos(angle)
                                    angle = acos(angle);

                                    // Check O(1,i)-H(2,j)-O(2,k) angle
                                    if (angle > PI/2) { // Must be at least 90 deg. = PI/2 rad

                                        // Check H(1)-O(1)-H(2) angle
                                        angle = 0;
                                        for (k = 0; k < natom; k++) {
                                            if (mol[i] == mol[k]) {
                                                if (ttype[k] == type2) {
                                                    // Determine vectors for H(1,k)-O(1,i)-H(2,j) angle
                                                    // Vector 1 is the same but opposite in sign
                                                    vec1[0]*=-1;
                                                    vec1[1]*=-1;
                                                    vec1[2]*=-1;
                                                    // Vector 2: pos(k)-pos(i)
                                                    // x component
                                                    dist1 = fabs(pos[k][0] - pos[i][0]);
                                                    dist2 = fabs(pos[k][0] - pos[i][0] - dx);
                                                    dist3 = fabs(pos[k][0] - pos[i][0] + dx);
                                                    if ((dist1 < dist2) && (dist1 < dist3)) {
                                                        vec2[0] = pos[k][0] - pos[i][0]; // No x adjustment
                                                    } else if ((dist2 < dist1) && (dist2 < dist3)) {
                                                        vec2[0] = pos[k][0] - pos[i][0] - dx; // + dx
                                                    } else {
                                                        vec2[0] = pos[k][0] - pos[i][0] + dx; // - dx
                                                    }
                                                    // y component
                                                    dist1 = fabs(pos[k][1] - pos[i][1]);
                                                    dist2 = fabs(pos[k][1] - pos[i][1] - dy);
                                                    dist3 = fabs(pos[k][1] - pos[i][1] + dy);
                                                    if ((dist1 < dist2) && (dist1 < dist3)) {
                                                        vec2[1] = pos[k][1] - pos[i][1]; // No y adjustment
                                                    } else if ((dist2 < dist1) && (dist2 < dist3)) {
                                                        vec2[1] = pos[k][1] - pos[i][1] - dy; // + dy
                                                    } else {
                                                        vec2[1] = pos[k][1] - pos[i][1] + dy; // - dy
                                                    }
                                                    // z component
                                                    dist1 = fabs(pos[k][2] - pos[i][2]);
                                                    dist2 = fabs(pos[k][2] - pos[i][2] - dz);
                                                    dist3 = fabs(pos[k][2] - pos[i][2] + dz);
                                                    if ((dist1 < dist2) && (dist1 < dist3)) {
                                                         vec2[2] = pos[k][2] - pos[i][2]; // No z adjustment
                                                    } else if ((dist2 < dist1) && (dist2 < dist3)) {
                                                         vec2[2] = pos[k][2] - pos[i][2] - dz; // + dz
                                                    } else {
                                                         vec2[2] = pos[k][2] - pos[i][2] + dz; // - dz
                                                    }
                                                    angle = 1; // Flag for found atom
                                                    break; // Break from loop
                                                }
                                            }
                                        }
                                        if (angle == 1) { // Check for flag that atom was found
                                            // Determine angle H(1)-O(1)-H(2)
                                            // dot product of vectors
                                            angle = vec1[0]*vec2[0];
                                            angle+=(vec1[1]*vec2[1]);
                                            angle+=(vec1[2]*vec2[2]);
                                            // magnitude of vector 1
                                            veclen = pow(vec1[0],2);
                                            veclen+=pow(vec1[1],2);
                                            veclen+=pow(vec1[2],2);
                                            angle/=sqrt(veclen);
                                            // magnitude of vector 2
                                            veclen = pow(vec2[0],2);
                                            veclen+=pow(vec2[1],2);
                                            veclen+=pow(vec2[2],2);
                                            angle/=sqrt(veclen);
                                            // acos(angle)
                                            angle = acos(angle);
               
                                            // Angle must be greater than 90 deg. = pi/2 radians
                                            if (angle > PI/2) {

                                                nhb++; // Increment count of HB

                                                // Axis of HB for atom 1
                                                axis[0] = relpos[0] - pos[i][0];
		                                axis[1] = relpos[1] - pos[i][1];
		                                axis[2] = relpos[2] - pos[i][2];

		                                // Project static polarization onto HB of atom 1
                                                // length of static polarization vector
		                                veclen = pow(stapol[i][0],2);
		                                veclen+=pow(stapol[i][1],2);
		                                veclen+=pow(stapol[i][2],2);
                                                veclen = sqrt(veclen);
                                                // dot product of static polariation vector and axis of HB
		                                proj_sta1 = stapol[i][0] * axis[0];
		                                proj_sta1+=(stapol[i][1]*axis[1]);
		                                proj_sta1+=(stapol[i][2]*axis[2]);
		                                proj_sta1/=veclen;

		                                // Project induced polarization onto HB of atom 1
                                                // length of induced polarization vector
                                                veclen = pow(indpol[i][0],2);
		                                veclen+=pow(indpol[i][1],2);
		                                veclen+=pow(indpol[i][2],2);
		                                veclen = sqrt(veclen);
                                                // dot product of induced polarization vector and axis of HB
		                                proj_ind1 = indpol[i][0] * axis[0];
		                                proj_ind1+=(indpol[i][1]*axis[1]);
		                                proj_ind1+=(indpol[i][2]*axis[2]);
		                                proj_ind1/=veclen;

		                                // Axis of HB for atom 2
		                                axis[0] = pos[i][0] - relpos[0];
		                                axis[1] = pos[i][1] - relpos[1];
		                                axis[2] = pos[i][2] - relpos[2];
		    
                                                // Project static polarization onto HB of atom 2
                                                // length of static polarization vector
                                                veclen = pow(stapol[j][0],2);
		                                veclen+=pow(stapol[j][1],2);
		                                veclen+=pow(stapol[j][2],2);
		                                veclen = sqrt(veclen);
                                                // dot product of static polarization vector and axis of HB
		                                proj_sta2 = stapol[j][0]*axis[0];
		                                proj_sta2+=(stapol[j][1]*axis[1]);
		                                proj_sta2+=(stapol[j][2]*axis[2]);
		                                proj_sta2/=veclen;

		                                // Project induced polarization onto HB of atom 2
                                                // length of induced polarization vector
                                                veclen = pow(indpol[j][0],2);
		                                veclen+=pow(indpol[j][1],2);
		                                veclen+=pow(indpol[j][2],2);
		                                veclen = sqrt(veclen);
                                                // dot product of induced polarization vector and axis of HB
		                                proj_ind2 = indpol[j][0]*axis[0];
		                                proj_ind2+=(indpol[j][1]*axis[1]);
		                                proj_ind2+=(indpol[j][2]*axis[2]);
		                                proj_ind2/=veclen;
                               
		                                // Report
                                                fprintf(currhb_mol,"%d %d\n",mol[i],mol[j]); // Report molecules
	                                        fprintf(currhb_atom,"%d %d\n",i+1,j+1); // Report atoms
                                                fprintf(currsta,"%lf %lf\n",proj_sta1,proj_sta2); // Report static polarization projection
		                                fprintf(currind,"%lf %lf\n",proj_ind1,proj_ind2); // Report induced polarization projection
                                                fprintf(currlen,"%lf\n",dist); // Report current length of HB
                                            }
                                        }
                                    }
		                }
                            }
                        }
                    }
                }
            } else if (ttype[i] == type2) { // Compare current atom type i to tinker type of hydrogen in HB
	        for (j = i+1; j < natom; j++) { // Loop i+1 to N
		    if (mol[i] != mol[j]) { // If not the same molecule
		        if (ttype[j] == type1) { // Compare current atom type j to tinker type of oxygen in HB
                            // Check x distance vs boundary
		            dist1 = fabs(pos[i][0] - pos[j][0]);
                            dist2 = fabs(pos[i][0] - pos[j][0] - dx);
                            dist3 = fabs(pos[i][0] - pos[j][0] + dx);
                            if ((dist1 < dist2) && (dist1 < dist3)) {
                                dist = pow(dist1,2); // No x adjustment
                            } else if ((dist2 < dist1) && (dist2 < dist3)) {
                                dist = pow(dist2,2); // + dx
                            } else {
                                dist = pow(dist3,2); // - dx
                            }
			    // Check y distance vs boundary
                            dist1 = fabs(pos[i][1] - pos[j][1]);
                            dist2 = fabs(pos[i][1] - pos[j][1] - dy);
                            dist3 = fabs(pos[i][1] - pos[j][1] + dy);
                            if ((dist1 < dist2) && (dist1 < dist3)) {
                                dist+=pow(dist1,2); // No y adjustment
                            } else if ((dist2 < dist1) && (dist2 < dist3)) {
                                dist+=pow(dist2,2); // + dy
                            } else {
                                dist+=pow(dist3,2); // - dy
                            }
			    // Check z distance vs boundary
                            dist1 = fabs(pos[i][2] - pos[j][2]);
                            dist2 = fabs(pos[i][2] - pos[j][2] - dz);
                            dist3 = fabs(pos[i][2] - pos[j][2] + dz);
                            if ((dist1 < dist2) && (dist1 < dist3)) {
                                dist+=pow(dist1,2); // No z adjustment
                            } else if ((dist2 < dist1) && (dist2 < dist3)) {
                                dist+=pow(dist2,2); // + dz
                            } else {
                                dist+=pow(dist3,2); // - dz
                            }
			    // Determine actual distance
                            dist = sqrt(dist);

                            // Compare distance to maximal allowed distance
                            if (dist < maxdist) {

                                // Determine relative position of atom 2
                                // x position
		                dist1 = fabs(pos[i][0] - pos[j][0]);
		                dist2 = fabs(pos[i][0] - pos[j][0] - dx);
		                dist3 = fabs(pos[i][0] - pos[j][0] + dx);
		                if ((dist1 < dist2) && (dist1 < dist3)) {
		                    relpos[0] = pos[j][0];
		                } else if ((dist2 < dist1) && (dist2 < dist3)) {
			            relpos[0] = pos[j][0] + dx;
		                } else {
			            relpos[0] = pos[j][0] - dx;
		                }
                                // y position
		                dist1 = fabs(pos[i][1] - pos[j][1]);
		                dist2 = fabs(pos[i][1] - pos[j][1] - dy);
		                dist3 = fabs(pos[i][1] - pos[j][1] + dy);
		                if ((dist1 < dist2) && (dist1 < dist3)) {
			            relpos[1] = pos[j][1];
                                } else if ((dist2 < dist1) && (dist2 < dist3)) {
			            relpos[1] = pos[j][1] + dy;
		                } else {
			            relpos[1] = pos[j][1] - dy;
		                }
                                // z position
		                dist1 = fabs(pos[i][2] - pos[j][2]);
		                dist2 = fabs(pos[i][2] - pos[j][2] - dz);
		                dist3 = fabs(pos[i][2] - pos[j][2] + dz);
		                if ((dist1 < dist2) && (dist1 < dist3)) {
                                    relpos[2] = pos[j][2];
		                } else if ((dist2 < dist1) && (dist2 < dist3)) {
			            relpos[2] = pos[j][2] + dz;
		                } else {
			            relpos[2] = pos[j][2] - dz;
		                }
                                    
                                // Check angle between vectors-> H(1,i)...O(2,j)-H(2,k)
                                angle = 0;
                                for (k = 0; k < natom; k++) {
                                    if (mol[j] == mol[k]) { // If j and k are in the same molecule
                                        if (ttype[k] == type2) { // If k is hydrogen
                                            // Determine H(i)...O(j)-H(k) angle
                                            // Vetor 1: H(i)...O(j)
                                            vec1[0] = pos[i][0] - relpos[0];
                                            vec1[1] = pos[i][1] - relpos[1];
                                            vec1[2] = pos[i][2] - relpos[2];
                                            // Vector 2: O(j)-H(k)
                                            // x component
                                            dist1 = fabs(pos[k][0] - pos[j][0]);
                                            dist2 = fabs(pos[k][0] - pos[j][0] - dx);
                                            dist3 = fabs(pos[k][0] - pos[j][0] + dx);
                                            if ((dist1 < dist2) && (dist1 < dist3)) {
                                                vec2[0] = pos[k][0] - pos[j][0]; // No x adjustment
                                            } else if ((dist2 < dist1) && (dist2 < dist3)) {
                                                vec2[0] = pos[k][0] - pos[j][0] - dx; // + dx
                                            } else {
                                                vec2[0] = pos[k][0] - pos[j][0] + dx; // - dx
                                            }
                                            // y component
                                            dist1 = fabs(pos[k][1] - pos[j][1]);
                                            dist2 = fabs(pos[k][1] - pos[j][1] - dy);
                                            dist3 = fabs(pos[k][1] - pos[j][1] + dy);
                                            if ((dist1 < dist2) && (dist1 < dist3)) {
                                                vec2[1] = pos[k][1] - pos[j][1]; // No y adjustment
                                            } else if ((dist2 < dist1) && (dist2 < dist3)) {
                                                vec2[1] = pos[k][1] - pos[j][1] - dy; // + dy
                                            } else {
                                                vec2[1] = pos[k][1] - pos[j][1] + dy; // - dy
                                            }
                                            // z component
                                            dist1 = fabs(pos[k][2] - pos[k][2]);
                                            dist2 = fabs(pos[k][2] - pos[k][2] - dz);
                                            dist3 = fabs(pos[k][2] - pos[k][2] + dz);
                                            if ((dist1 < dist2) && (dist1 < dist3)) {
                                                vec2[2] = pos[k][2] - pos[j][2]; // No z adjustment
                                            } else if ((dist2 < dist1) && (dist2 < dist3)) {
                                                vec2[2] = pos[k][2] - pos[j][2] - dz; // + dz
                                            } else {
                                                vec2[2] = pos[k][2] - pos[j][2] + dz; // - dz
                                            }
                                            angle = 1; // Flag for angle found
                                            break; // Break from loop
                                        }
                                    }
                                }

                                if (angle == 1) { // Check for flag
                                    // Determine angle between vectors
                                    // dot product of vectors
                                    angle = vec1[0]*vec2[0];
                                    angle+=(vec1[1]*vec2[1]);
                                    angle+=(vec1[2]*vec2[2]);
                                    // magnitude of vector 1
                                    veclen = pow(vec1[0],2);
                                    veclen+=pow(vec1[1],2);
                                    veclen+=pow(vec1[2],2);
                                    angle/=sqrt(veclen);
                                    // magnitude of vector 2
                                    veclen = pow(vec2[0],2);
                                    veclen = pow(vec2[1],2);
                                    veclen = pow(vec2[2],2);
                                    angle/=sqrt(veclen);
                                    // acos(angle)
                                    angle = acos(angle);
                                    if (angle > PI/2) {

                                        // Check second angle between vectors: O(1,k)-H(1,i)...O(2,j)
                                        angle = 0;
                                        for (k = 0; k < natom; k++) {
                                            if (mol[k] == mol[i]) {
                                                if (ttype[k] == type1) {
                                                    // Determine vectors
                                                    // Vector 1: pos(j)-pos(i) -> same as previous but opposite sign
                                                    vec1[0]*=-1;
                                                    vec1[1]*=-1;
                                                    vec1[2]*=-1;
                                                    // Vector 2: pos(k) - pos(i)
                                                    // x component
                                                    dist1 = fabs(pos[k][0] - pos[i][0]);
                                                    dist2 = fabs(pos[k][0] - pos[i][0] - dx);
                                                    dist3 = fabs(pos[k][0] - pos[i][0] + dx);
                                                    if ((dist1 < dist2) && (dist1 < dist3)) {
                                                        vec2[0] = pos[k][0] - pos[i][0]; // No x adjustment
                                                    } else if ((dist2 < dist1) && (dist2 < dist3)) {
                                                        vec2[0] = pos[k][0] - pos[i][0] - dx; // + dx
                                                    } else {
                                                        vec2[0] = pos[k][0] - pos[i][0] + dx; // - dx
                                                    }
                                                    // y component
                                                    dist1 = fabs(pos[k][1] - pos[i][1]);
                                                    dist2 = fabs(pos[k][1] - pos[i][1] - dy);
                                                    dist3 = fabs(pos[k][1] - pos[i][1] + dy);
                                                    if ((dist1 < dist2) && (dist1 < dist3)) {
                                                        vec2[1] = pos[k][1] - pos[i][1]; // No y adjustment
                                                    } else if ((dist2 < dist1) && (dist2 < dist3)) {
                                                        vec2[1] = pos[k][1] - pos[i][1] - dy; // + dy
                                                    } else {
                                                        vec2[1] = pos[k][1] - pos[i][1] + dy; // - dy
                                                    }
                                                    // z component
                                                    dist1 = fabs(pos[k][2] - pos[i][2]);
                                                    dist2 = fabs(pos[k][2] - pos[i][2] - dz);
                                                    dist3 = fabs(pos[k][2] - pos[i][2] + dz);
                                                    if ((dist1 < dist2) && (dist1 < dist3)) {
                                                        vec2[2] = pos[k][2] - pos[i][2]; // No z adjustment
                                                    } else if ((dist2 < dist1) && (dist2 < dist3)) {
                                                        vec2[2] = pos[k][2] - pos[i][2] - dz; // + dz
                                                    } else {
                                                        vec2[2] = pos[k][2] - pos[i][2] + dz; // - dz
                                                    }
                                                    angle = 1; // Flag for atom found
                                                    break; // Break from loop
                                                }
                                            }
                                        }

                                        if (angle == 1) { // Check for atom found flag
                                            // Determine angle between vectors
                                            // dot product
                                            angle = vec1[0]*vec2[0];
                                            angle*=(vec1[1]*vec2[1]);
                                            angle*=(vec1[2]*vec2[2]);
                                            // magnitude of vector 1
                                            veclen = pow(vec1[0],2);
                                            veclen+=pow(vec1[1],2);
                                            veclen+=pow(vec1[2],2);
                                            angle/=sqrt(veclen);
                                            // magnitude of vector 2
                                            veclen = pow(vec2[0],2);
                                            veclen+=pow(vec2[1],2);
                                            veclen+=pow(vec2[2],2);
                                            angle/=sqrt(veclen);
                                            // acos(angle)
                                            angle = acos(angle);

                                            if (angle > PI/2) { // Angle should be geater than 90 deg = pi/2 radians
                                                nhb++; // Increment count of HB

		                                // Determine HB axis for atom 1
                                                axis[0] = relpos[0] - pos[i][0];
		                                axis[1] = relpos[1] - pos[i][1];
		                                axis[2] = relpos[2] - pos[i][2];

                                                // Determine static polarization projection onto HB for atom 1
                                                // length of static polarization vector
		                                veclen = pow(stapol[i][0],2);
		                                veclen+=pow(stapol[i][1],2);
		                                veclen+=pow(stapol[i][2],2);
		                                veclen = sqrt(veclen);
                                                // dot product of static polarization vector and axis of HB
		                                proj_sta1 = stapol[i][0]*axis[0];
		                                proj_sta1+=(stapol[i][1]*axis[1]);
                                                proj_sta1+=(stapol[i][2]*axis[2]);
		                                proj_sta1/=veclen;

		                                // Determine induced polarization projection onto HB for atom 1
                                                // length of induced polarization vector
                                                veclen = pow(indpol[i][0],2);
		                                veclen+=pow(indpol[i][1],2);
		                                veclen+=pow(indpol[i][2],2);
		                                veclen = sqrt(veclen);
                                                // dot product of induced polarization vector and axis of HB
		                                proj_ind1 = indpol[i][0]*axis[0];
		                                proj_ind1+=(indpol[i][1]*axis[1]);
                                                proj_ind1+=(indpol[i][2]*axis[2]);
		                                proj_ind1/=veclen;

                                                // Determine HB axis for atom 2
		                                axis[0] = pos[i][0] - relpos[0];
		                                axis[1] = pos[i][1] - relpos[1];
		                                axis[2] = pos[i][2] - relpos[2];

		                                // Determine static polarization projection onto HB for atom 2
                                                // length of static polarization vector
		                                veclen = pow(stapol[j][0],2);
		                                veclen+=pow(stapol[j][1],2);
		                                veclen+=pow(stapol[j][2],2);
		                                veclen = sqrt(veclen);
                                                // dot product of static polarization vector and axis of HB
		                                proj_sta2 = stapol[j][0]*axis[0];
		                                proj_sta2+=(stapol[j][1]*axis[1]);
                                                proj_sta2+=(stapol[j][2]*axis[2]);
		                                proj_sta2/=veclen;

		                                // Determine induced polarization projection onto HB for atom 2
                                                // length of induced polarization vector
                                                veclen = pow(indpol[j][0],2);
		                                veclen+=pow(indpol[j][1],2);
		                                veclen+=pow(indpol[j][2],2);
		                                veclen = sqrt(veclen);
                                                // dot product of induced polarization vector and axis of HB
		                                proj_ind2 = indpol[j][0]*axis[0];
		                                proj_ind2+=(indpol[j][1]*axis[1]);
                                                proj_ind2+=(indpol[j][2]*axis[2]);
		                                proj_ind2/=veclen;

		                                // Report
                                                fprintf(currhb_mol,"%d %d\n",mol[i],mol[j]); // Report molecules
                                                fprintf(currhb_atom,"%d %d\n",i+1,j+1); // Report atoms
		                                fprintf(currsta,"%lf %lf\n",proj_sta1,proj_sta2); // Report static polarization projection
		                                fprintf(currind,"%lf %lf\n",proj_ind1,proj_ind2); // Report induced polarization projection
		                                fprintf(currlen,"%lf\n",dist); // Report current length of HB
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
	 	}
	    }
	}

        // Close current time step report files
	fclose(currhb_atom);
	fclose(currhb_mol);
	fclose(currsta);
	fclose(currind);
	fclose(currlen);

        // Record number of hydrogen bonds
	fprintf(ct_hb,"%d\n",nhb);

	// Zero out for next time step
        nhb = 0;

    }

    printf("Analysis complete.\n"); // Report to terminal

    // Clean up -------------------------------------------------------

    // Close files
    fclose(arc);
    fclose(ct_hb);
    
    // De-allocate arrays
    free(ttype);
    free(mol);
    for (i = 0; i < natom; i++) {
	free(pos[i]);
    }
    free(pos);
    for (i = 0; i < natom; i++) {
	free(stapol[i]);
    }
    free(stapol);
    for (i = 0; i < natom; i++) {
	free(indpol[i]);
    }
    free(indpol);

    // Final exit -----------------------------------------------------

    printf("Complete.\n");
    exit(EXIT_SUCCESS);

}
