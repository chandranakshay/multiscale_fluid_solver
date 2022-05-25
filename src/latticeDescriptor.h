#ifndef _LATTICEDESCRIPTOR_
#define _LATTICEDESCRIPTOR_

/** @file */

/*! \brief Defines the lattice.
 *  It defines the lattice based on the input parameters 
 *  - number of grid points per processor 
 *  - number of processors in each direction 
 *  - number of ghost nodes (Dummy nodes on each side of the grid. Will be used for boundary conditions)  
 */

typedef struct latticeInfo {
    
    /*! \brief Default constructor 
     * Looks for the files simparam.data which should contain the following data.
     * It should have the number of real points, number of processors and ghost nodes for X,Y and Z directions in the same order.
     * Eg: Nx,Ny,Nz,Px,Py,Pz,Gx,Gy,Gz
     */    
    latticeInfo()  
    {
        std::ifstream fileSimParam; 
        fileSimParam.open("simParam.data");
        
        if(!fileSimParam)throw("Input param file missing"); 
        fileSimParam >> n1;     
        fileSimParam >> n2;
        fileSimParam >> n3;
        fileSimParam >> nP1;    
        fileSimParam >> nP2;
        fileSimParam >> nP3; 
        fileSimParam >> nGhost1; 
        fileSimParam >> nGhost2;
        fileSimParam >> nGhost3;
        
        fileSimParam.close(); 
    }
    
    /*! \brief Copy Constructor
     *  Takes the lattice parameters as input directly instead of reading them from the simparam.data file.
     */    
    latticeInfo(unsigned long long int N1,unsigned long long int N2,unsigned long long int N3,unsigned long long int NP1,unsigned long long int NP2,unsigned long long int NP3,unsigned long long int NG1,unsigned long long int NG2,unsigned long long int NG3)
    {
        n1 = N1;
        n2 = N2;
        n3 = N3;
        
        nP1 = NP1;
        nP2 = NP2;
        nP3 = NP3;
        
        nGhost1 = NG1;
        nGhost2 = NG2;
        nGhost3 = NG3;
    }
    // Number of Grid points on local processor in each Direction
    unsigned long long int n1, n2, n3;
    // Number of Processors in Each Direction
    int nP1, nP2, nP3;
    // Ghost Nodes
    unsigned long long int nGhost1,nGhost2,nGhost3; 
} latticeDescriptor;



#endif
