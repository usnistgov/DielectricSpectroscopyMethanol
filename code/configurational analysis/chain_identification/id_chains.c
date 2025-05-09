#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

/*

    This code takes as input the series of files in which the hydrogen 
    bonds at each time step of the simulation are identified by the 
    molecule ID of the molecules involved in each bond. From this input, 
    this code determines the chains of hydrogen bonded molecules that 
    exist at each time step of the simulation. Note: this code does 
    not track these chains through time.



 *       Expected input:
 *
 *       nmol      ->    number of molecules in simulation
 *
         Expected output:

              Files named

                   chainmol# where # is 1 through total time steps in simulation


 -Rebecca A. Bone

 */

// Subroutines -----------------------------------------------------------
// -----------------------------------------------------------------------

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

// Main function --------------------------------------------------------
// ----------------------------------------------------------------------

void main(int argc,char* argv[]) {

    // Print header -----------------------------------------------------

    printf("Tracking chains of hydrogen bonded molecules through time.\n");

    // User input handling -----------------------------------------------

    // Variable declarations
    int nmol,ct; // Integers
    char fline[1200]; // Character array to facilitate file reading

    // Handle provided inputs or request from command line
    if (argc == 2) { // Correct number of command line inputs
	
	// nmol    ->    number of molecules in simulation
	nmol = atoi(argv[1]); // Parse to integer
	ct = 0; // Zero out counter o cycle iterations
	while (nmol < 1) { // Number of molecules must be at least 1
            // Handle cycle counting
	    ct++; // Increment cycle counter
	    if (ct >= 10) { // Check for hitting threshold of cycle iterations
		printf("ERROR:\n Requested user input too many times.\n"); // Report to terminal
		exit(EXIT_FAILURE); // Exit
	    }
	    // Request and validate user input
	    printf("Number of molecules: "); // Request input from user
	    scanf("%s",&fline); // Read command line input
	    nmol = atoi(fline); // Parse to integer
	}
    } else { // Incorrect number of command line arguments
	
	// nmol    ->    number of molecules in simulation
	printf("Number of molecules: "); // Request input from user
	scanf("%s",&fline); // Read command line input
	nmol = atoi(fline); // Parse to integer
	ct = 0; // Zero out counter of cycle iterations
	while (nmol < 1) { // Number of molecules must be at least 1
	    // Handle cycle counting
	    ct++; // Increment cycle counter
	    if (ct >= 10) { // Check for hitting threshold of cycle iterations
		printf("ERROR:\n Requested user input too many times.\n"); // Report to terminal
		exit(EXIT_FAILURE); // Exit
	    }
	    // Request and validate user input
	    printf("Number of molecules: "); // Request input from user
	    scanf("%s",&fline); // Read command line input
	    nmol = atoi(fline); // Parse to integer
	}
    }

    // Initialization ----------------------------------------------------

    // Variable declarations
    int i,j,k; // Integers

    // Allocate array for number index of extant chains (nmol x 1)
    int *index; // Pointer to array
    index = (int*) malloc(nmol*sizeof(int)); // Allocate array
    // Zero out allocated array
    for (i = 0; i < nmol; i++) {
	index[i] = 0;
    }

    // Allocate flag array for marking extant chains as identified at current time step
    int *found; // Pointer to array
    found = (int*) malloc(nmol*sizeof(int)); // Allocate array
    // Zero out allocated array
    for (i = 0; i < nmol; i++) {
	found[i] = 0;
    }

    // Allocate array for molecules in extant chains (nmol x nmol)
    int **extant_molec; // Pointer to array
    extant_molec = (int**) malloc(nmol*sizeof(int*)); // Allocate first level of array
    for (i = 0; i < nmol; i++) { // Loop through first level of array
	extant_molec[i] = (int*) malloc(nmol*sizeof(int)); // Allocate second level of array 
    }
    // Zero out allocated array
    for (i = 0; i < nmol; i++) {
	for (j = 0; j < nmol; j++) {
	    extant_molec[i][j] = 0;
	}
    }

    // Allocate array for hydrogen bonds in extant chains (nmol x 2*nmol x 2)
    int ***extant_hb; // Pointer to array
    extant_hb = (int***) malloc(nmol*sizeof(int**)); // Allocate first level of array
    for (i = 0; i < nmol; i++) { // Loop through first level of array
	extant_hb[i] = (int**) malloc(5*nmol*sizeof(int*)); // Allocate second level of array
    }
    for (i = 0; i < nmol; i++) { // Loop through first level of array
	for (j = 0; j < 5*nmol; j++) { // Loop through second level of array
	    extant_hb[i][j] = (int*) malloc(2*sizeof(int)); // Allocate third level of array
	}
    }
    // Zero out allocated array
    for (i = 0; i < nmol; i++) {
	for (j = 0; j < 5*nmol; j++) {
	    for (k = 0; k < 2; k++) {
		extant_hb[i][j][k] = 0;
	    }
	}
    }

    //printf("Test 1.\n"); // DEBUG : report to terminal

    // Allocate array for reading in chains at current time step (nmol x 1)
    int *chain; // Pointer to array
    chain = (int*) malloc(nmol*sizeof(int)); // Allocate array
    // Zero out allocated array
    for (i = 0; i < nmol; i++) {
	chain[i] = 0;
    }
    
    // Allocate array for reading in hbs of chain at current time step (2*nmol x 2)
    int **hblist; // Pointer to array
    hblist = (int**) malloc(5*nmol*sizeof(int*)); // Allocate first level of array
    for (i = 0; i < 5*nmol; i++) { // Loop through first level
	hblist[i] = (int*) malloc(5*sizeof(int)); // Allocate second level of array
    }
    // Zero out allocated array
    for (i = 0; i < 5*nmol; i++) {
	for (j = 0; j < 2; j++) {
	    hblist[i][j] = 0;
	}
    }

    //printf("Test 2.\n"); // DEBUG  report to terminal

    // Allocate array for recording fully new chains to be recorded after full comparison to extant chains (nmol x nmol)
    int **new_molec; // Pointer to array
    new_molec = (int**) malloc(nmol*sizeof(int*)); // Allocate first level of array
    for (i = 0; i < nmol; i++) { // Loop through first level of array
	new_molec[i] = (int*) malloc(nmol*sizeof(int)); // Allocate second level of array
    }
    // Zero out allocated array
    for (i = 0; i < nmol; i++) {
	for (j = 0; j < nmol; j++) {
	    new_molec[i][j] = 0;
	}
    }

    // Allocate array for recording fully new hydrogen bonds in chains to be recorded after fully comparison to extant chains (nmol x nmol*2 x 2)
    int ***new_hb; // Pointer to array
    new_hb = (int***) malloc(nmol*sizeof(int**)); // Allocate first level of array
    for (i = 0; i < nmol; i++) { // Loop through first level of array
	new_hb[i] = (int**) malloc(5*nmol*sizeof(int*)); // Allocate second level of array
    }
    for (i = 0; i < nmol; i++) { // Loop through first level o array
	for (j = 0; j < 5*nmol; j++) { // Loop through second level of array
	    new_hb[i][j] = (int*) malloc(5*sizeof(int)); // Allocate third level of array
	}
    }
    // Zero out allocated array
    for (i = 0; i < nmol; i++) {
	for (j = 0; j < 5*nmol; j++) {
	   for (k = 0; k < 2; k++) {
	       new_hb[i][j][k] = 0;
	   }
        }
    }

    //printf("Test 3.\n"); // DEBUG : report to terminal

    // Allocate array to record index of extant chain whose history is mimicked in new chain (nmol x 1)
    int *mimic; // Pointer to array
    mimic = (int*) malloc(nmol*sizeof(int)); // Allocate array
    // Zero out allocated array
    for (i = 0; i < nmol; i++) {
	mimic[i] = 0;
    }

    printf("Initialiation complete.\n"); // Report to terminal

    // Establish extant lists at first time step -------------------------

    // Variable declarations
    int x,extant_count,new_count,max_index; // Integers
    char hbline[1200]; // Character array for reading in hydrogen bond lists
    char fname[200]; // Character array for constructing file names
    char delim[] = " "; // Expected file delimiter
    char hb_delim[] = ", "; // Expected hydrogen bond file sub-delimiter
    char *ptr; // Character pointer to assist parsing items of file lines

    // Open first time step chain molecules specification file
    sprintf(fname,"chainmol1"); // Construct expected file name
    x = chk_file(fname); // Check for file existence
    if (x == 0) { // No file found
	printf("ERROR:\n No initial chain molecules file found.\n"); // Report to terminal
	printf("Expected file name: %s\n",fname); // Report to terminal
	exit(EXIT_FAILURE); // Exit
    }
    FILE * molin1; // Pointer to file
    molin1 = fopen(fname,"r"); // Open file for reading
    if (molin1 == NULL) { // Check or correct file opening
	printf("ERROR:\n Initial file could not be opened.\n"); // Report to terminal
	printf("File name: %s\n",fname); // Report to terminal
	perror("Error"); // Report reason for error
	exit(EXIT_FAILURE);
    }

    // Open first time step chain hbs specification file
    sprintf(fname,"chainhb1"); // Construct expected file name
    x = chk_file(fname); // Check for file existence
    if (x == 0) { // No file found
	printf("ERROR:\n No initial chain hydrogen bonds file found.\n"); // Report to terminal
	printf("Expected file name: %s\n",fname); // Report to terminal
	exit(EXIT_FAILURE); // Exit
    }
    FILE * hbin1; // Pointer to file
    hbin1 = fopen(fname,"r"); // Open file for reading
    if (hbin1 == NULL) { // Check for correct file opening
	printf("ERROR:\n Initial file could not be opened.\n"); // Report to terminal
	printf("File name: %s\n",fname); // Report to terminal
	perror("ERROR"); // Report reasonf ro error
	exit(EXIT_FAILURE); // Exit
    }

    //printf("Initial files found.\n"); // DEBUG : report to terminal

    extant_count = 0; // Set count of extant chains to zero
    new_count = 0; // Set count of new chains to zero

    // Read through first time step data and assign as extant chains
    while (fgets(fline,1200,molin1) != NULL) { // While additional lines in file
	//printf("Extant chain count: %d\n",extant_count); // DEBUG : report to terminal
	// Read in hb data....and check for existence
	if (fgets(hbline,1200,hbin1) == NULL) { // Check that hb data also exists
            printf("ERROR:\n Incomplete initial chain hydrogen bonds file.\n"); // Report to terminal
	    printf("File name: %s\n",fname); // Report to terminal
	    exit(EXIT_FAILURE);
	}
	//printf("%s\n",fline); // DEBUG
	//printf("%s\n",hbline); // DEBUG
        // Record current chain molecules
	ptr = strtok(fline,delim); // Parse first entry of molecule line
	if (ptr == NULL) { // Empty line found
	    break; // Leave while loop
	}
	extant_molec[extant_count][0] = atoi(ptr); // Parse to integer and assign
	ct = 0; // Count of molecules in current chain
	ptr = strtok(NULL,delim); // Parse next line item
	while (ptr != NULL) { // Continue until no more items in current line
	    ct++; // Increment count of molecules in current chain
	    extant_molec[extant_count][ct] = atoi(ptr); // Parse to integer and assign
	    ptr = strtok(NULL,delim); // Parse next line item
	}
	x = ct; // Record number of molecules in current chain for reading hydrogen bonds
	// Fill remainder of current row with non-molecule id values
	while (ct < nmol) { // Loop through remaining row positions
	    ct++; // Increment count
	    extant_molec[extant_count][ct] = -1; // Assign non-molecule flag value
	}
	if (x != 0) { // Check for number of molecules in chain...vs single molecule
            // Actual chain
	    //printf("Chain located.\n"); // DEBUG : report to terminal
	    //printf("Length: %d\n",x+1); // DEBUG : report to terminal
            // Print to file molecules in current chain at first time step
            sprintf(fname,"chain%d_1mol",extant_count+1); // Construct file name
            FILE * molecout1; // Pointer to file
            molecout1 = fopen(fname,"w"); // Open file for writing
	    if (molecout1 == NULL) { // Check for correct file opening
		printf("ERROR:\n File could not be opened.\n"); // Report to terminal
		printf("File name: %s\n",fname); // Report to terminal
		perror("ERROR"); // Report reason for error
		exit(EXIT_FAILURE); // Exit
	    }
            for (i = 0; i < x; i++) { // Loop through all but last molecule of current chain
		//printf("i=%d\n",i); // DEBUG : report to terminal
                fprintf(molecout1,"%d ",extant_molec[extant_count][i]); // Print to file
		//printf("%d ",extant_molec[extant_count][i]); // DEBUG
            }
	    //printf("Last term print.\n"); // DEBUG : report to terminal
            fprintf(molecout1,"%d\n",extant_molec[extant_count][x]); // Print last molecule to file
	    //printf("%d\n",extant_molec[extant_count][x]); // DEBUG
            fclose(molecout1); // Close file
	    //printf("Molecules reported.\n"); // DEBUG : report to terminal
	    // Record current chain hydrogen bonds
	    ptr = strtok(hbline,hb_delim); // Parse first entry of hydrogen bond line
	    if (ptr == NULL) { // Empty line found
	        printf("ERROR:\n Incomplete hydrogen bond information.\n"); // Report to terminal
		printf("File name: chainhb1\n"); // Report to terminal
	        exit(EXIT_FAILURE); // Exit
	    }
            extant_hb[extant_count][0][0] = atoi(ptr); // Parse to integer and assign
	    ptr = strtok(NULL,hb_delim); // Parse second molecule in current hydrogen bond
	    extant_hb[extant_count][0][1] = atoi(ptr); // Parse to integer and assign
	    ct = 0; // Count of hydrogen bonds in current chain
	    ptr = strtok(NULL,hb_delim); // Parse next line item
            while (ptr != NULL) { // Continue until no more items in current line
                ct++; // Increment count of hydrogen bonds in current chain
	        extant_hb[extant_count][ct][0] = atoi(ptr); // Parse to integer and assign
	        ptr = strtok(NULL,hb_delim); // Parse next line item
	        extant_hb[extant_count][ct][1] = atoi(ptr); // Parse to integer and assign
	        ptr = strtok(NULL,hb_delim); // Parse next line item
	    }
            x = ct;
	    while (ct < nmol-1) {
	        ct++; // Increment position
	        extant_hb[extant_count][ct][0] = -1; // Mark with non-molecule value
	        extant_hb[extant_count][ct][1] = -1; // Mark with non-molecule value
	    }
            // Record hydrogen bonds in current chain at first time step
            sprintf(fname,"chain%d_1hb",extant_count+1); // Construct file name
            FILE * hbout1; // Pointer to file
            hbout1 = fopen(fname,"w"); // Open file for writing
	    if (hbout1 == NULL) { // Check for correct file opening
                printf("ERROR:\n File could not be opened.\n"); // Report to terminal
		printf("File name: %s\n",fname); // Report to terminal
		perror("ERROR"); // Report reason for error
		exit(EXIT_FAILURE); // Exit
	    }
            for (i = 0; i < x+1; i++) {
                fprintf(hbout1,"%d,%d\n",extant_hb[extant_count][i][0],extant_hb[extant_count][i][1]); // Print
		//printf("%d,%d ",extant_hb[extant_count][i][0],extant_hb[extant_count][i][1]); // DEBUG
            }
	    //printf("\n"); // DEBUG
            fclose(hbout1); // Close file
	    // Record index of current chain
	    index[extant_count] = extant_count+1;
	    //printf("Index %d of chain %d\n",index[extant_count],extant_count); // DEBUG
            extant_count++; // Increment count of extant chains
	    // Record starting point as first time step
	    sprintf(fname,"chain%d_start",extant_count); // Construct file name
	    FILE * startout1; // Pointer to file
	    startout1 = fopen(fname,"w"); // Open file for writing
	    if (startout1 == NULL) { // Check or correct file opening
		printf("ERROR:\n File could not be opened.\n"); // Report to terminal
		printf("File name: %s\n",fname); // Report to terminal
		perror("ERROR"); // Report reason for error
		exit(EXIT_FAILURE); // Exit
	    }
	    fprintf(startout1,"1"); // Report initial time frame as start
	    fclose(startout1); // Close file
	}
    }

    max_index = extant_count; // Record current extant count as max index of found

    // Close first time step files
    fclose(hbin1);
    fclose(molin1);

    printf("First time frame read in.\n");
    printf("%d extant chains.\n",extant_count);

    // DEBUG : print starting extant moleucles and hydrogen bonds
    /*for (i = 0; i < extant_count; i++) {
	printf("Chain %d: ",i+1);
	for (j = 0; j < nmol; j++) {
	    if (extant_molec[i][j] != -1) {
		printf("%d ",extant_molec[i][j]);
	    } else {
		printf("\n");
		break;
	    }
	}
	printf("HB: ");
	for (j = 0; j < 5*nmol; j++) {
	    if (extant_hb[i][j][0] == -1) {
		printf("\n");
		break;
	    } else {
		printf("%d,%d ",extant_hb[i][j][0],extant_hb[i][j][1]);
	    }
	}
    }*/

    // Analysis ----------------------------------------------------------

    // Variable deciarations
    int comp,t; // Integers

    t = 1; // Set initial time frame
    comp = 0; // Set flag for end of files

    while (comp == 0) {
        
        t++; // Increment time rame
	printf("t=%d\n",t); // Report to terminal

        /*if (t == 151) { // DEBUG : end early
	    break;
	}*/

        if (extant_count > nmol) {
	    printf("WARNING:\n Invalid number of extant chains.\n");
	//} else {
	    //printf("%d extant chains.\n",extant_count);
	}

	// Zero out counters
	new_count = 0;

	// Open current time step molecules in chains
	sprintf(fname,"chainmol%d",t); // Construct expected file name
	x = chk_file(fname); // Check for file existence
	if (x == 0) { // No file found
	    comp = 1; // Flip flag for no more files
	    printf("End of files found at time t=%d\n",t); // Report to terminal
	    break; // Break from while loop
	}
	FILE * molin; // Pointer to file
	molin = fopen(fname,"r"); // Open file for reading
	if (molin == NULL) { // Check for correct file opening
	    printf("ERROR:\n File could not be opened.\n"); // Report to terminal
	    printf("File name: %s\n",fname); // Report to terminal
	    perror("Error"); // Report reason for error
	    exit(EXIT_FAILURE); // Exit
	}

	// Open current time step hydrogen bonds in chains
	sprintf(fname,"chainhb%d",t); // Construct expected file name
	x = chk_file(fname); // Check for file existence
	if (x == 0) { // No file found
	    comp = 1; // Flip flag for no more files
	    printf("End of files found at time t=%d\n",t); // Report to terminal
	    break; // Break from while loop
	}
	FILE * hbin; // Pointer to file
	hbin = fopen(fname,"r"); // Open file for reading
	if (hbin == NULL) { // Check or correct file opening
	    printf("ERROR:\n File could not be opened.\n"); // Report to terminal
	    printf("File name: %s\n",fname); // Report to terminal
	    perror("Error"); // Report reason for error
	    exit(EXIT_FAILURE); // Exit
	}

        // Make sure found is ze3roed and mimic is -1'ed
	for (i = 0; i < nmol; i++) {
	    found[i] = 0;
	    mimic[i] = -1;
	}

        // Read in molecules and hydrogen bonds and compare to extant chains
	while (fgets(fline,1200,molin) != NULL) { // While additional lines of file
            // Retrieve hydrogen bond line or flag and exit loop
	    if (fgets(hbline,1200,hbin) == NULL) { // Missing hydrogen bond line for molecule line
                printf("ERROR:\n Missing hydrogen bond data.\n"); // Report to terminal
		comp = -1; // Flip flag for end of files
		break; // Leave while loop
	    }
	    //printf("New chain raw: %s\n",fline); // DEBUG
	    //printf("New hb raw: %s\n",hbline); // DEBUG
	    // Parse molecules into chain array
            ptr = strtok(fline,delim); // Parse first line entry
	    chain[0] = atoi(ptr); // Parse and assign first value
	    ptr = strtok(NULL,delim); // Grab next line entry
	    ct = 0; // Zero out chain length
	    while (ptr != NULL) { // While additional members of current line
                ct++; // Increment chain length
		chain[ct] = atoi(ptr); // Parse and assign next member of chain
		ptr = strtok(NULL,delim); // Grab next line entry
	    }
	    x = ct; // Record number of molecules in current chain
	    if (x > 0) { // Only continue analysis if more than one molecule
		//printf("It's a chain.\n"); // DEBUG : report to terminal
		//printf("Length %d\n",x); // DEBUG : report to terminal
	        // Mark remaining parts of chain array as non-molecules
	        while (ct < nmol-1) { // Continue through remaining elements to mark as non-molecules
		    ct++;
		    chain[ct] = -1;
	        }
	        // Parse hydrogen bonds into array
	        ptr = strtok(hbline,hb_delim); // Parse first line entry
	        if (ptr != NULL) {
	            hblist[0][0] = atoi(ptr); // Parse and assign first value
	            ptr = strtok(NULL,hb_delim); // Parse second part of first entry
		    hblist[0][1] = atoi(ptr); // Parse and assign second value of first entry
		    ct = 0; // Set count of hydrogen bonds involved
		    ptr = strtok(NULL,hb_delim); // Parse next entry
		    while (ptr != NULL) {
	                ct++; // Increment number of hydrogen bonds involved
		        hblist[ct][0] = atoi(ptr); // Parse and assign first part of hydrogen bond
		        ptr = strtok(NULL,hb_delim); // Parse second half
		        hblist[ct][1] = atoi(ptr); // Parse and assign second part of hydrogen bond
		        ptr = strtok(NULL,hb_delim); // Parse next hydrogen bond part 1
		    }
	        } else {
		    ct = -1; // Mark spot for no hydrogen bonds
	        }
	        // Mark remaining spots in array as non-molecules
	        while (ct < 5*nmol-1) {
		    ct++; // Increment position
		    hblist[ct][0] = -1;
		    hblist[ct][1] = -1;
	        }

                // DEBUG : check read in and parse of molecules and hb
		/*printf("New molecules: ");
		for (i = 0; i < nmol; i++) {
		    if (chain[i] != -1) {
			printf("%d ",chain[i]);
		    } else {
			printf("\n");
			break;
		    }
		}
		printf("New hbs: ");
		for (i = 0; i < 5*nmol; i++) {
		    if (hblist[i][0] != -1) {
			printf("%d,%d ",hblist[i][0],hblist[i][1]);
		    } else {
			printf("\n");
			break;
	            }
		}*/

	        // Compare read in chain to extant chains
		x = -1; // Flag for if current chain has been found at all
		for (i = 0; i < extant_count; i++) { // Loop through extant chains
		    if (found[i] != -1) {
		        ct = 0; // Count of hb in common
		        for (j = 0; j < 5*nmol; j++) { // Loop through HB list of extant chain
			    if (extant_hb[i][j][0] != -1) {
			        for (k = 0; k < 5*nmol; k++) { // Loop through HB list of read in chain
				    if (hblist[k][0] != -1) {
				        if ((extant_hb[i][j][0] == hblist[k][0]) && (extant_hb[i][j][1] == hblist[k][1])) { // HB matches
					    ct++; // Increment count of matching HB
				        } else if ((extant_hb[i][j][0] == hblist[k][1]) && (extant_hb[i][j][1] == hblist[k][0])) { // HB matches
					    ct++; // Increment count of matching HB
				        }
				    } else {
				        break; // Break from loop of read in HB list
				   }
			        }
			    } else {
			        break; // Break from loop of extant HB list
			    }
		        }
		        if (ct > 0) { // If there are HB in common with current extant chain
			    if (x == -1) { // Current chain has not been found as part of other extant chain
				//printf("Linked to extant chain %d.\n",index[i]); // DEBUG : report to terminal
			        x = i; // Record spot first found
			        if (found[i] == 0) { // If identified extant chain has not been found elsewhere yet
				    found[i] = 1; // Mark extant chain as found
				    // Print current chain molecules as next part of extant chain history
				    sprintf(fname,"chain%d_%dmol",index[i],t); // Construct file name
				    FILE * molout; // Pointer to file
				    molout = fopen(fname,"w"); // Open file for writing
				    if (molout == NULL) { // Check for correct file opening
					printf("ERROR:\n File could not be opened.\n"); // Report to terminal
					printf("File name: %s\n",fname); // Report to terminal
					perror("Error"); // Report reason for error
					exit(EXIT_FAILURE); // Exit
				    }
				    for (j = 0; j < nmol; j++) {
				        if (chain[j+1] != -1) {
					    fprintf(molout,"%d ",chain[j]); // Intermediate print
				        } else {
					    fprintf(molout,"%d\n",chain[j]); // Final print
					    break; // Break from printing for loop
				        }
				    }
				    fclose(molout); // Close file
				    // Print current chain hydrogen bonds as next part of extant chain history
				    sprintf(fname,"chain%d_%dhb",index[i],t); // Construct file name
				    FILE * hbout; // Pointer to file
				    hbout = fopen(fname,"w"); // Open ile or writing
				    if (hbout == NULL) { // Check for correct file opening
					printf("ERROR:\n File could not be opened.\n"); // Report to terminal
					printf("File name: %s\n",fname); // Report to terminal
					perror("Error"); // Report reason for error
					exit(EXIT_FAILURE); // Exit
				    }
			            for (j = 0; j < 5*nmol; j++) {
				        if (hblist[j][0] != -1) {
					    fprintf(hbout,"%d %d\n",hblist[j][0],hblist[j][1]); // Print
				        } else {
					    break; // Break from printing for loop
				        }
				    }
				    fclose(hbout); // Close file
			        } else if (found[i] == 1) { // Identified extant chain has been found previously
			            //printf("Found fragment of extant chain %d.\n",index[i]); // DEBUG : report to terminal
				    // Will record new entry of extant that is split history of this extant chain
				    // Print current chain to new chain list
				    for (j = 0; j < nmol; j++) {
				        new_molec[new_count][j] = chain[j];
				    }
				    // Print current chain HB to new chain HB list
				    for (j = 0; j < 5*nmol; j++) {
				        new_hb[new_count][j][0] = hblist[j][0];
				        new_hb[new_count][j][1] = hblist[j][1];
				    }
				    // Record index of extant chain for mimicking
				    mimic[new_count] = i;
				    // Increment number of new chains
				    new_count++;
			        } // Or do nothing if extant chain marked for extinction
			    } else { // Current chain is result of joinder of two chains
				// Mark extant chain for extinction
				//printf("Extinction.\n"); // DEBUG : report to terminal
				found[i] = -1;
			    }
		        } // Do nothing if no HB in common
		    } // Do nothing if extant chain marked for extinction
		} // Finish looking through extant chains
		if (x == -1) { // If current chain was not found in any extant chain
	            //printf("Brand new chain.\n"); // DEBUG : report to terminal
	            // Record chain molecules to new chain list
		    for (i = 0; i < nmol; i++) {
			new_molec[new_count][i] = chain[i];
		    }
		    // Record chain hbs to new hb list
		    for (i = 0; i < 5*nmol; i++) {
			new_hb[new_count][i][0] = hblist[i][0];
			new_hb[new_count][i][1] = hblist[i][1];
		    }
                    // Increment count of new chains
		    new_count++;
                    //printf("New count: %d\n",new_count); // DEBUG : report to terminal
		}
	    } else {
		//printf("Single molecule.\n"); // DEBUG : report to terminal
	    }
        }

	// Close input files
	fclose(hbin);
	fclose(molin);

        //printf("Chain search done.\n"); // DEBUG : report to terminal

	// Check for end of files flag
	if (comp != 0) { // Check for end of file
	    printf("End of files found at time t=%d\n",t); // Report to terminal
	    break; // Break from main while loop
	}

        // Remove extant chains that were not found in current time frame
        for (i = extant_count-1; i >= 0; i--) {
	    if (found[i] != 1) {
		// Record ending time of chain
		sprintf(fname,"chain%d_end",index[i]); // Construct file name
		FILE * endout; // Pointer to file
		endout = fopen(fname,"w"); // Open file for writing
		if (endout == NULL) { // Check for correct file opening
		    printf("ERROR:\n File could not be opened.\n"); // Report to terminal
		    printf("File name: %s\n",fname); // Report to terminal
		    perror("Error"); // Report reason for error
		    exit(EXIT_FAILURE); // Exit
		}
		fprintf(endout,"%d\n",t); // Print last time of chain as previous time step
		fclose(endout); // Close file
		// Move extant chain molecules up to overwrite ended chain
		if (i != extant_count-1) {
		    for (j = i+1; j < extant_count; j++) {
			for (k = 0; k < nmol; k++) {
			    extant_molec[j-1][k] = extant_molec[j][k];
			}
		    }
		}
		// Move extant chain hydrogen bonds up to overwrite ended chain
		if (i != extant_count-1) {
		    for (j = i+1; j < extant_count; j++) {
			for (k = 0; k < 5*nmol; k++) {
			    extant_hb[j-1][k][0] = extant_hb[j][k][0];
			    extant_hb[j-1][k][1] = extant_hb[j][k][1];
			}
		    }
		}
		// Move index list up to overwrite ended chain
		if (i != extant_count-1) {
		    for (j = i+1; j < extant_count; j++) {
			index[j-1] = index[j];
		    }
		}
		// Check if any mimicking above current point
		if (new_count > 0) {
		    for (j = 0; j < new_count; j++) {
			if (mimic[j] > i) {
			    mimic[j]--;
			}
		    }
		}

		// Decrement number of extant chains
		extant_count--;
	    }
	}

        //printf("Extinct chains removed.\n"); // DEBUG : report to terminal
        //printf("Extant chains: %d\n",extant_count); // DEBUG

	// Add new chains to extant chain list
        if (new_count > 0) {
	    //printf("New chains: %d\n",new_count); // DEBUG
	    if (new_count > nmol) {
		printf("WARNING:\n More new chains reported than molecules.\n"); // Report to terminal
	    }
	    for (i = 0; i < new_count; i++) {
	        if (mimic[i] == -1) { // If chain is completely new
		    //printf("New chain %d is completely new and recorded as chain %d.\n",i+1,max_index+1); // DEBUG
		    // Record molecule list to extant molecule list
		    for (j = 0; j < nmol; j++) {
			extant_molec[extant_count][j] = new_molec[i][j];
		    }
		    // Record hydrogen bond list to extant hydrogen bond list
		    for (j = 0; j < 5*nmol; j++) {
			extant_hb[extant_count][j][0] = new_hb[i][j][0];
			extant_hb[extant_count][j][1] = new_hb[i][j][1];
		    }
		    // Record indexing number
		    index[extant_count] = max_index+1;
		    max_index++; // Increment max index count
		    // Increment number of extant chains
		    extant_count++;
		    // Record start time of new extant chain
		    sprintf(fname,"chain%d_start",index[extant_count-1]); // Construct file name
		    FILE * startout; // Pointer to file
		    startout = fopen(fname,"w"); // Open file for writing
		    if (startout == NULL) { // Check for correct file opening
			printf("ERROR:\n File could not be opened.\n"); // Report to terminal
			printf("File name: %s\n",fname); // Report to terminal
			perror("Error"); // Report reason for error
			exit(EXIT_FAILURE); // Exit
	            }
		    fprintf(startout,"%d\n",t); // Print current time frame as start point
		    fclose(startout); // Close file
		    // Print molecule list of extant chain
		    sprintf(fname,"chain%d_%dmol",index[extant_count-1],t); // Construct file name
		    FILE * molout_new; // Pointer to file
		    molout_new = fopen(fname,"w"); // Open file for writing
		    if (molout_new == NULL) { // Check for correct file opening
			printf("ERROR:\n File could not be opened.\n"); // Report to terminal
			printf("File name: %s\n",fname); // Report to terminal
			perror("Error"); // Report reason for error
			exit(EXIT_FAILURE); // Exit
		    }
		    for (j = 0; j < nmol; j++) {
			if (new_molec[i][j+1] == -1) {
			    fprintf(molout_new,"%d\n",new_molec[i][j]); // Final print
			    break; // Break from print loop
			} else {
			    fprintf(molout_new,"%d ",new_molec[i][j]); // Normal print
			}
		    }
		    fclose(molout_new); // Close file
		    // Print hydrogen bond list of extant chain
		    sprintf(fname,"chain%d_%dhb",index[extant_count-1],t); // Construct file name
		    FILE * hbout_new; // Pointer to file
		    hbout_new = fopen(fname,"w"); // Open file for writing
		    if (hbout_new == NULL) { // Check or correct file opening
			printf("ERROR:\n File could not be opened.\n"); // Report to terminal
			printf("File name: %s\n",fname); // Report to terminal
			perror("Error"); // Report reason or error
			exit(EXIT_FAILURE); // Exit
		    }
		    for (j = 0; j < 5*nmol; j++) {
			if (new_hb[i][j][0] != -1) {
			    fprintf(hbout_new,"%d %d\n",new_hb[i][j][0],new_hb[i][j][1]); // Print
			} else {
			    break; // Break from print loop
			}
		    }
		    fclose(hbout_new); // Close file
		} else { // If chain is splinter of existing chain
		   // printf("New chain %d is splinter of extant chain %d and recorded as chain %d.\n",i+1,index[mimic[i]],max_index+1); // DEBUG
		    // Record molecule list to extant molecule list
		    for (j = 0; j < nmol; j++) {
		        extant_molec[extant_count][j] = new_molec[i][j];
		    }
		    //printf("Molecule list recorded to extant list.\n"); // DEBUG
		    // Record hydrogen bond list to extant hydrogen bond list
		    for (j = 0; j < 5*nmol; j++) {
			extant_hb[extant_count][j][0] = new_hb[i][j][0];
			extant_hb[extant_count][j][1] = new_hb[i][j][1];
		    }
		    //printf("Hydrogen bond list recorded to extant list.\n"); // DEBUG
		    // Record new index
		    index[extant_count] = max_index+1;
		    max_index++; // Increment max index
		    //printf("New highest chain index recorded to extant index list.\n"); // DEBUG
		    // Increment number of extant chains
		    extant_count++;
		    // Record start time of ``new'' extant chain from mimicked chain
		    // Grab start time of chain being mimicked
		    sprintf(fname,"chain%d_start",index[mimic[i]]); // Construct file name or mimicked chain
		    x = chk_file(fname); // Check for file existence
		    if (x == 0) {
			printf("ERROR:\n Missing starting time of chain file.\n"); // Report to terminal
			printf("File name: %s\n",fname); // Report to terminal
			perror("Error"); // Report reason for error
			exit(EXIT_FAILURE); // Exit
	            }
		    FILE * startin_mimic; // Pointer to file
		    startin_mimic = fopen(fname,"r"); // Open file for reading
		    if (startin_mimic == NULL) { // Check or proper file opening
			printf("ERROR:\n File not properly opened.\n"); // Report to terminal
			printf("File name: %s\n",fname); // Report to terminal
			perror("Error"); // Report reason for error
			exit(EXIT_FAILURE); // Exit
		    }
		    fgets(fline,1200,startin_mimic); // Grab line
		    ptr = strtok(fline,delim); // Parse first line entry
		    x = atoi(ptr); // Starting time of mimicked chain
		    fclose(startin_mimic); // Close file
		    //printf("Starting time frame of mimicked chain is %d.\n",x); // DEBUG
		    // Report mimicked chain start time to ``new'' chain
		    sprintf(fname,"chain%d_start",index[extant_count-1]); // Construct file name
		    FILE * startout_mimic; // Pointer to file
		    startout_mimic = fopen(fname,"w"); // Open file for writing
		    if (startout_mimic == NULL) { // Check for correct file opening
			printf("ERROR:\n File could not be opened.\n"); // Report to terminal
			printf("File name: %s\n",fname); // Report to terminal
			perror("Error"); // Report reason for error
			exit(EXIT_FAILURE); // Exit
		    }
		    fprintf(startout_mimic,"%d\n",x); // Print mimicked start time to file
		    fclose(startout_mimic); // Close file
		    //printf("Starting time of splinter chain reported.\n"); // DEBUG
		    // Print molecule lists at previous times mimicking extant chain
		    for (j = x; j < t; j++) { // Loop through time mimicked chain existed up to now
			// Open file to be copied
			sprintf(fname,"chain%d_%dmol",index[mimic[i]],j); // Construct jth time file name
			ct = chk_file(fname); // Check for file existence
			if (ct == 0) {
			    printf("ERROR:\n Missing molecules of chain file.\n"); // Report to terminal
			    printf("File name: %s\n",fname); // Report to terminal
			    exit(EXIT_FAILURE); // Exit
			}
			FILE * molin_mimic; // Pointer to file
			molin_mimic = fopen(fname,"r"); // Open file for reading
			if (molin_mimic == NULL) { // Check or correct file opening
			    printf("ERROR:\n File not opened as expeted.\n"); // Report to terminal
                            printf("File name: %s\n",fname); // Report to terminal
			    perror("Error"); // Report reason for error
			    exit(EXIT_FAILURE); // Exit
			}
			//printf("Molecule file of original chain at time %d opened.\n",j); // DEBUG
			// Open file being copied to
			sprintf(fname,"chain%d_%dmol",index[extant_count-1],j); // Construct jth time file name
			FILE * molout_mimic; // Pointer to file
			molout_mimic = fopen(fname,"w"); // Open file for writing
			if (molout_mimic == NULL) { // Check for correct file opening
			    printf("ERROR:\n File not opened as expected.\n"); // Report to terminal
			    printf("File name: %s\n",fname); // Report to terminal
			    perror("Error"); // Report reason for error
			    exit(EXIT_FAILURE); // Exit
			}
			//printf("Molecule file of new chain at time %d opened.\n",j); // DEBUG
			// Read data, parse, and print to new file
			fgets(fline,1200,molin_mimic); // Retrieve line from file
			//printf("Molecules from original chain: %s\n",fline);
			ptr = strtok(fline,delim); // Parse first line entry
			//printf("First molecule: %d\n",atoi(ptr)); // DEBUG
			if (ptr == NULL) {
			    printf("ERROR:\n Invalid molecules of chain file.\n"); // Report to terminal
			    printf("File name: %s\n",fname); // Report to terminal
			    perror("Error"); // Report reason for error
			    exit(EXIT_FAILURE);
			}
			fprintf(molout_mimic,"%d ",atoi(ptr));
			ptr = strtok(NULL,delim); // Parse next line entry
			while (ptr != NULL) { // While additional line entries
		            //printf("Next molecule: %d\n",atoi(ptr)); // DEBUG
			    fprintf(molout_mimic,"%d ",atoi(ptr)); // Print to file
			    ptr = strtok(NULL,delim); // Parse next line entry
			}
			// Close files
			fclose(molin_mimic);
			fclose(molout_mimic);
                    }
		    //printf("Historic extant molecule lists printed for splinter chain.\n"); // DEBUG
		    // Print current molecule list as new molecule list
		    sprintf(fname,"chain%d_%dmol",index[extant_count-1],t); // Construct file name
		    FILE * molout_mimiclast; // Pointer to file
		    molout_mimiclast = fopen(fname,"w"); // Open file for writing
		    if (molout_mimiclast == NULL) { // Check for correct file opening
			printf("ERROR:\n File could not be opened.\n"); // Report to terminal
			printf("File name: %s\n",fname); // Report to terminal
			perror("Error"); // Report reason for error
			exit(EXIT_FAILURE); // Exit
		    }
		    for (j = 0; j < nmol; j++) {
			if (new_molec[i][j+1] == -1) {
			    fprintf(molout_mimiclast,"%d\n",new_molec[i][j]); // Last print
			    break; // Break rom print loop
			} else {
			    fprintf(molout_mimiclast,"%d ",new_molec[i][j]); // Normal print
			}
		    }
		    fclose(molout_mimiclast); // Close file
		    //printf("Current time splinter chain molecule list reported.\n"); // DEBUG
		    // Print hydrogen bond lists at previous times mimicking extant chain
		    for (j = x; j < t; j++) { // Loop through times mimicked chain existed
			// Open mimicked chain hydrogen bond file
			sprintf(fname,"chain%d_%dhb",index[mimic[i]],j); // Construct file name
			ct = chk_file(fname); // Check for file existence
			if (ct == 0) {
			    printf("ERROR:\n Missing hydrogen bond information.\n"); // Report to terminal
			    printf("File name: %s\n",fname); // Report to terminal
			    exit(EXIT_FAILURE); // Exit
			}
			FILE * hbin_mimic; // Pointer to file
			hbin_mimic = fopen(fname,"r"); // Open file for reading
			if (hbin_mimic == NULL) { // Check or correct opening
			    printf("ERROR:\n File not opened as expected.\n"); // Report to terminal
			    printf("File name: %s\n",fname); // Report to terminal
			    perror("Error"); // Report reason for error
			    exit(EXIT_FAILURE); // Exit
			}
			//printf("Mimic chain hydrogen bond info at time %d opened.\n",j); // DEBUG
			// Open new hydrogen bond chain file
			sprintf(fname,"chain%d_%dhb",index[extant_count-1],j); // Construct file name
			FILE * hbout_mimic; // Pointer to file
			hbout_mimic = fopen(fname,"w"); // Open file for writing
			if (hbout_mimic == NULL) { // Check for ile opened correctly
		            printf("ERROR:\n File not opened as expected.\n"); // Report to terminal
			    printf("File name: %s\n",fname); // Report to terminal
			    perror("Error"); // Report reason for error
			    exit(EXIT_FAILURE); // Exit
			}
			//printf("Splinter chain hydrogen bond info at time %d opened.\n",j); // DEBUG
			// Read, parse, and print
		        while (fgets(fline,1200,hbin_mimic) != NULL) { // While we have additional lines
			    ptr = strtok(fline,hb_delim); // Parse first line entry
			    fprintf(hbout_mimic,"%d ",atoi(ptr)); // Parse to integer and print
                            ptr = strtok(NULL,hb_delim); // Parse second line entry
			    fprintf(hbout_mimic,"%d\n",atoi(ptr)); // Parse to integer and print
			}
			// Close files
			fclose(hbin_mimic);
			fclose(hbout_mimic);
		    }
		    //printf("Previous time hydrogen bond reporting complete.\n"); // DEBUG
		    // Print current hydrogen bond list as new hydrogen bond list
                    sprintf(fname,"chain%d_%dhb",index[extant_count-1],t); // Construct file name
		    FILE * hbout_mimiclast; // Pointer to file
		    hbout_mimiclast = fopen(fname,"w"); // Open file for writing
		    if (hbout_mimiclast == NULL) { // Check for correct file opening
			printf("ERROR:\n File could not be opened.\n"); // Report to terminal
			printf("FIle name: %s\n",fname); // Report to terminal
			perror("Error"); // Report reason for error
			exit(EXIT_FAILURE); // Exit
		    }
		    for (j = 0; j < 5*nmol; j++) {
			if (new_hb[i][j][0] != -1) {
			    fprintf(hbout_mimiclast,"%d %d\n",new_hb[i][j][0],new_hb[i][j][1]); // Print
			} else {
			    break; // Break from print loop
			}
		    }
		    fclose(hbout_mimiclast); // Close file
		    //printf("Current time hydrogen bond reporting complete.\n"); // DEBUG
		}
	    }
	}

        //printf("New chains added.\n"); // DEBUG : report to terminal
        //printf("Extant chains: %d\n",extant_count); // DEBUG

        // Zero arrays
	for (i = 0; i < nmol; i++) {
	    found[i] = 0;
	    mimic[i] = -1;
	}

        // Check for current time extant chain reporting files
	for (i = 0; i < extant_count; i++) {
	    sprintf(fname,"chain%d_%dmol",index[i],t); // Build file name in character array
	    x = chk_file(fname); // Check for file existence
	    if (x == 0) {
		printf("ERROR:\n Missing reporting file for chain %d at time %d.\n",index[i],t); // Report to terminal
		printf("File name: %s\n",fname); // Report to terminal
		exit(EXIT_FAILURE); // Exit
	    }
	    sprintf(fname,"chain%d_%dhb",index[i],t); // Build file name in character array
	    x = chk_file(fname); // Check for file existence
	    if (x == 0) {
		printf("ERROR:\n Missing reporting file for chain %d at time %d.\n",index[i],t); // Report to terminal
		printf("File name: %s\n",fname); // Report to terminal
		exit(EXIT_FAILURE); // Exit
	    }
	}

	// Check indeces vs max index
	for (i = 0; i < nmol; i++) {
	    if (index[i] > max_index) {
		printf("ERROR:\n Invalid index of %dth chain as %d when max index is %d.\n",i+1,index[i],max_index); // Report to terminal
		exit(EXIT_FAILURE); // Exit
	    }
	}

    }

    // Clear extant list
    if (extant_count > 0) {
	for (i = 0; i < extant_count; i++) {
	    sprintf(fname,"chain%d_end",index[i]);
	    FILE * last_end;
	    last_end = fopen(fname,"w");
	    fprintf(last_end,"%d",t);
	    fclose(last_end);
	}
    }

    // Clean up ----------------------------------------------------------

    // Close opened files

    // De-allocate allocated arrays
    free(index);
    free(found);
    for (i = 0; i < nmol; i++) {
	free(extant_molec[i]);
    }
    free(extant_molec);
    for (i = 0; i < nmol; i++) {
	for (j = 0; j < 5*nmol; j++) {
	    free(extant_hb[i][j]);
	}
    }
    for (i = 0; i < nmol; i++) {
	free(extant_hb[i]);
    }
    free(extant_hb);
    free(chain);
    for (i = 0; i < 5*nmol; i++) {
	free(hblist[i]);
    }
    free(hblist);
    for (i = 0; i < nmol; i++) {
	free(new_molec[i]);
    }
    free(new_molec);
    for (i = 0; i < nmol; i++) {
	for (j = 0; j < 5*nmol; j++) {
	    free(new_hb[i][j]);
	}
    }
    for (i = 0; i < nmol; i++) {
	free(new_hb[i]);
    }
    free(new_hb);
    free(mimic);

    // Final report and exit
    printf("Complete.\n"); // Report to terminal
    exit(EXIT_SUCCESS); // Exit

}
