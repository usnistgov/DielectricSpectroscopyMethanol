#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

/*

    This code takes as input the identified chains of hydrogen 
    bonded molecules at each time step of the simulation. This 
    code then tracks these chains through time. At each time step, 
    this code compares the chains that exist to the chains that 
    existed at the previous time step of the simulation. If there 
    are hydrogen bonds held in common, the chains are related in 
    history. However, an extant chain may split into multiple 
    chains or merge with other chains. Thus, a more nuanced algorithm 
    is applied to account for these overlapping histories. Addition 
    information is available in the accompanying publication.


    Expected input:

        nmol -> number of molecules in simulation

    Expected output:

        series of files for each chain identified:

            chain#_start -> contains starting index of chain #
            chain#_end   -> contains ending index of chain #
            chain#_#mol  -> contains molecule IDs of molecules in
                            chain # at time step #



-Dr. Rebecca A. Bone

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
    return pacc; // Return flag
}

// Exit request check
int chk_exit(char name[]) {
    // Function takes string array containing ile name and locations
    // Function outputs flag signalling if indicated file exists
    int request = 0;
    if (strcmp(name,"exit") == 0) {
	request = 1;
    } else if (strcmp(name,"Exit") == 0) {
	request = 1;
    } else if (strcmp(name,"EXIT") == 0) {
	request = 1;
    } else if (strcmp(name,"quit") == 0) {
	request = 1;
    } else if (strcmp(name,"Quit") == 0) {
	request = 1;
    } else if (strcmp(name,"QUIT") == 0) {
	request = 1;
    }
    return request;
}

// Main function-------------------------------------------------------
// --------------------------------------------------------------------

void main(int argc,char* argv[]) {

    printf("Analysis code for chains and loops of molecules in 1D HB networks.\n");

    // Declarations ---------------------------------------------------

    int i,j,k,q,t,comp,nmol,mol1,mol2,x,ct,nstack,current,nchain,lenchain;
    char fname[100];
    char fline[1200];
    char delim[] = " ";
    
    // User input -----------------------------------------------------

    // Need: nmol (number of molecules in simulation)

    if (argc == 1) { // 0/1 inputs 
        printf("Number of molecules in simulation: "); // Request user input
	scanf("%s",&fline);
	nmol = atoi(fline);
	if (nmol < 1) {
	    nmol = 0;
	}
	ct = 0; // Zero out cycle request counter
        while (nmol == 0) {
            ct++; // Increment cycle counter
            if (ct >= 10) {
                printf("ERROR:\n Too many input requests.\n"); // Report reason for error
                printf("Exiting.\n"); // Report exit
                exit(EXIT_FAILURE); // Exit
            }
            printf("Number of molecules in simulation: "); // Request user input
            scanf("%s",&fline); // User command line input
            if (chk_exit(fline) == 1) { // Check for user exit request
                printf("Exit requested.\n"); // Report to terminal
                exit(EXIT_FAILURE); // Exit
            }
            nmol = atoi(fline);
            if (nmol < 1) { // Number of atoms must be positive
                nmol = 0;
            }
        }
    } else if (argc == 2) { // 1/1 inputs
	sprintf(fline,"%s",argv[2]); // Print second input to character array
	nmol = atoi(fline); // Parse to integer
	if (nmol < 1) { // Number of atoms must be positive
	    nmol = 0;
	}
	ct = 0; // Zero out request cycle counter
	while (nmol == 0) {
	    ct++; // Increment cycle counter
	    if (ct >= 10) {
		printf("ERROR:\n Too many input requests.\n"); // Report reason for error
		printf("Exiting.\n"); // Report exit
		exit(EXIT_FAILURE); // Exit
            }
	    printf("Number of molecules in simulation: "); // Request user input
	    scanf("%s",&fline); // User command line input
	    if (chk_exit(fline) == 1) { // Check for user exit request
		printf("Exit requested.\n"); // Report to terminal
		exit(EXIT_FAILURE); // Exit
	    }
	    nmol = atoi(fline);
	    if (nmol < 1) { // Number of atoms must be positive
	        nmol = 0;
	    }
	}
    } else if (argc > 2) { // 2+/1 inputs
        printf("Number of molecules in simulation: "); // Request user input
        scanf("%s",&fline);
        nmol = atoi(fline);
        if (nmol < 1) {
            nmol = 0;
        }
        ct = 0; // Zero out cycle request counter
        while (nmol == 0) {
            ct++; // Increment cycle counter
            if (ct >= 10) {
                printf("ERROR:\n Too many input requests.\n"); // Report reason for error
                printf("Exiting.\n"); // Report exit
                exit(EXIT_FAILURE); // Exit
            }
            printf("Number of molecules in simulation: "); // Request user input
            scanf("%s",&fline); // User command line input
            if (chk_exit(fline) == 1) { // Check for user exit request
                printf("Exit requested.\n"); // Report to terminal
                exit(EXIT_FAILURE); // Exit
            }
            nmol = atoi(fline);
            if (nmol < 1) { // Number of atoms must be positive
                nmol = 0;
            }
        }
    }

    // Initialization -------------------------------------------------

    printf("Initializing.\n"); // Report to terminal

    // Allocate array for HB connections between molecules (2D: nmol x nmol)
    int** connect; // Pointer
    connect = (int**) malloc(nmol*sizeof(int*)); // First level pointers
    for (i = 0; i < nmol; i++) {
	connect[i] = (int*) malloc(nmol*sizeof(int)); // Second level pointers
    }
    // Zero out allocated array
    for (i = 0; i < nmol; i++) {
	for (j = 0; j < nmol; j++) {
	    connect[i][j] = 0;
	}
    }

    // Allocate array for indicating if in DFS the current node has been found (1D: nmol)
    int *found; // Declare pointer
    found = (int*) malloc(nmol*sizeof(int)); // Allocation
    // Zero out allocated array
    for (i = 0; i < nmol; i++) {
	found[i] = 0;
    }

    // Allocate array for indicating parent node in DFS of each node (1D: nmol)
    int *parent; // Declare pointer
    parent = (int*) malloc(nmol*sizeof(int)); // Allocation
    // Zero out allocated array
    for (i = 0; i < nmol; i++) {
	parent[i] = 0;
    }

    // Allocate array for stack for DFS (1D: 10*nmol)
    int *stack; // Declare pointer
    stack = (int*) malloc(10*nmol*sizeof(int)); // Allocation
    // Zero out allocated array
    for (i = 0; i < 10*nmol; i++) {
	stack[i] = 0;
    }

    // Allocate array for chain reconstrution (1D: nmol)
    int *chainlist;
    chainlist = (int*) malloc(nmol*sizeof(int)); // Allocation
    // Zero out allocated array
    for (i = 0; i < nmol; i++) {
	chainlist[i] = 0;
    }

    // Open file for number of chains reporting
    sprintf(fname,"number_chains_hb");
    FILE * out_chain;
    out_chain = fopen(fname,"w");
    if (out_chain == NULL) {
	printf("ERROR:\n File to report number of chains of HB did not open correctly.\n");
	perror("Failed- ");
	exit(EXIT_FAILURE);
    }

    // Analysis -------------------------------------------------------

    printf("Beginning analysis.\n"); // Report to terminal

    comp = 0; // Flag for no more files
    t = 0; // Time frame index

    while ((comp == 0) && (t < 200000001)) {

        t++; // Increment time frame number
	printf("t=%d\n",t);

	// Construct expected HB analysis file name
	sprintf(fname,"hbmol%d",t);
	x = chk_file(fname); // Ensure current file exists
	if (x == 0) {
	    comp = 1;
	    break;
	}

	// Open current HB analysis file
	FILE * hbout; // Pointer to file
	hbout = fopen(fname,"r"); // Open file for reading
        if (hbout == NULL) {
	    printf("ERROR:\n Hbmol file not opened successfully.\n");
	    perror("Failed- ");
	    exit(EXIT_FAILURE);
	}

        // Read in HB Molecules and record to connection array
	x = 0; // Flag for end of file
	while (x == 0) {
	    if (fgets(fline,1200,hbout) == NULL) {
		x = 1;
		break;
	    } else {
		char *ptr = strtok(fline,delim);
                mol1 = atoi(ptr);
		ptr = strtok(NULL,delim);
		mol2 = atoi(ptr);
		connect[mol1-1][mol2-1] = 1;
		connect[mol2-1][mol1-1] = 1;
	    }
	}

	fclose(hbout);

        //printf("Current HB read in.\n"); // DEBUG : report to terminal

        for (i = 0; i < nmol; i++) {
	    if (found[i] == 0) {
		for (j = 0; j < 10*nmol; j++) {
		    stack[j] = 0;
		}
		stack[0] = i;
		nstack = 1;

	        // Breadth First Search
	        while (nstack > 0) { // Continue while stack is not empty

                    // Push first index of stack
	            current = stack[0];

	            // Move stack to disclude first index
	            for (i = 1; i < nstack; i++) {
		        stack[i-1] = stack[i];
	            }
	            stack[nstack] = 0;
	            nstack--;

	            // Check for previously found
                    if (found[current] == 0) {
		        found[current] = 1;

		        // Check for parent already identified
		        if (parent[current] == 0) {
		            parent[current] = current;
		            nchain++;
		        }

                        // Search connections of current
                        for (i = 0; i < nmol; i++) {
	                    if (connect[current][i] == 1) {
			        if (found[i] == 0) {
                                    for (j = nstack; j > 0; j--) {
			                stack[j] = stack[j-1];
			            }
			            stack[0] = i;
			            nstack++;
			            if (parent[i] == 0) {
				        parent[i] = current;
				    }
			        }
			    }
		        }
		    }
	        }
	    }
        }
	
        //printf("BFS complete.\n"); // DEBUG : report to terminal

        // Open file to report chain lengths
        sprintf(fname,"chainlen%d",t);
	FILE * outlen;
	outlen = fopen(fname,"w");

	// Open file to report reconstructed chains
	sprintf(fname,"chains_t%d",t);
	FILE * reconst;
	reconst = fopen(fname,"w");

	// Analyze parent list to reconstruct chains
	nchain = 0;
	for (i = 0; i < nmol; i++) {
	    if (parent[i] == i) {
		// First element of new chain found
		nchain++;
		chainlist[0] = i;
		lenchain = 1;
		// Search for new elements along parent list
		x = 0;
		while (x == 0) {
		    x = 1;
		    for (j = i+1; j < nmol; j++) {
			if (parent[j] == chainlist[lenchain-1]) {
			    if (lenchain == 1) {
				chainlist[lenchain] = j;
				lenchain++;
				x = 0;
				break;
			    } else {
				if (chainlist[lenchain-2] != j) {
                                    if (chainlist[0] == j) {
					x = -1;
					break;
				    } else {
				        chainlist[lenchain] = j;
				        lenchain++;
				        x = 0;
				        break;
				    }
				}
			    }
			}
		    }
		}
		// Extend parent list direction from initial chain element again through possible other molecule
		if (x != -1) {
		    x = 0;
		    while (x == 0) {
		        x = 1;
		        for (j = i+1; j < nmol; j++) {
			    if (parent[j] == chainlist[0]) {
			        if (chainlist[1] != j) {
		                    if (chainlist[lenchain-1] != j) {
				        // Shift chain list for newly found element
				        for (k = lenchain-1; k > 0; k--) {
				            chainlist[k] = chainlist[k-1];
				        }
				        chainlist[0] = j;
					x = 0;
					break;
				    } else {
					x = -1;
					break;
				    }
				}
			    }
			}
		    }
		}

		if (lenchain > 1) {
		    // Report length of chain
                    fprintf(outlen,"%d\n",lenchain);

                    if (lenchain > 2) {
		        // Report reconstructed chain
		        for (j = 1; j < lenchain-1; j++) {
	                    fprintf(reconst,"%d ",chainlist[j]+1);
		        }
		        fprintf(reconst,"%d\n",chainlist[lenchain-1]+1);
		    } else {
			fprintf(reconst,"%d %d\n",chainlist[0]+1,chainlist[1]+1);
		    }
		} else {
		    nchain--;
		}
            }
	}

	fclose(reconst);
        fclose(outlen);

        // Report number of chains
	fprintf(out_chain,"%d\n",nchain);	

        //printf("Reporting complete.\n"); // DEBUG : report to terminal

	// Re-zero connection array
	for (i = 0; i < nmol; i++) {
	    for (j = 0; j < nmol; j++) {
                connect[i][j] = 0;
	    }
	}

	// Re-zero parent array
	for (i = 0; i < nmol; i++) {
	    parent[i] = 0;
	}

	// Re-zero found array
	for (i = 0; i < nmol; i++) {
	    found[i] = 0;
	}

	// Re-zero chain reconstruction array
	for (i = 0; i < nmol; i++) {
	    chainlist[i] = 0;
	}

	//printf("Preparation for next time step complete.\n"); // DEBUG : report to terminal
    }

    // Clean up -------------------------------------------------------

    printf("End of input file located.\n"); // Report to terminal

    // Close opened files
    fclose(out_chain);

    // De-allocate arrays
    for (i = 0; i < nmol; i++) {
	free(connect[i]);
    }
    free(connect);
    free(found);
    free(parent);
    free(chainlist);

    printf("Complete.\n");
}
