#ifndef  _GRID_BCC_3D_BASIC_H_
#define  _GRID_BCC_3D_BASIC_H_ 

/** @file */

#include<mm_malloc.h>
#include<iostream>
#include<paniniAliases.h>

/*! \brief Defines the Grid.
 *  It defines the grid for each processor based on the lattice parameters 
 *  - It has m1 number of real points 
 *  - 4 Ghost nodes on each side of lattice i.e., 8 Ghost nodes per direction
 */
template <int numfield,int numblock, typename dataType> 
class gridBCC3D 
{
public:
    gridBCC3D(int size1=1,int size2=1,int size3 =1,int stride1=0,int stride2=0,int stride3=0)
    {
        m1       = size1;                   
        m2       = size2;
        m3       = size3;
        actualGridSize = m1*m2*m3;
        ndB1     = stride1;               
        ndB2     = stride2;
        ndB3     = stride3;
        nB1      = ndB1 + stride1;        
        nB2      = ndB2 + stride2;
        nB3      = ndB3 + stride3;
        n1       = size1 + 4*stride1;     
        n2       = size2 + 4*stride2;
        n3       = size3 + 4*stride3;
        nE1      = size1 + nB1 -1;        
        nE2      = size2 + nB2 -1;
        nE3      = size3 + nB3 -1;
        ndE1     = nE1 + stride1;       
        ndE2     = nE2 + stride2;
        ndE3     = nE3 + stride3;
        sizeGrid = n1*n2*n3;            
        numField = numfield; 
        numBlock = numblock;
        numElem  = sizeGrid*numField*numBlock ; /**< Total number of elemnts in the grid per Node or Cell*/

/*! \brief data in bytes reuired for the grid per processor.
 * 
 *   Crystallographic boltzmann assumes that there are two grids at an offset of 0.5 in each direction.
 *   Multiplication by two signifies that there are two grids (Node and Cell)
 *   For vectorization purpose aligned version of memory allocation with vector size as 4096
 */
        data = (dataType*) _mm_malloc(2*numElem*sizeof(dataType), 4096);
        
        for(int i=0; i<numBlock;i++)
        {
            node[i]=nodeTYPE::NODE;
            cell[i]=nodeTYPE::CELL;
        } 
        
        
        
    }
    
    ~gridBCC3D()
    { 
        delete data; 
    }
    
    unsigned long long int size() const
    {
        return numElem;
    }

/*! \brief returns an index value which is unique for a population at given point on a grid.
 *  The 3D grid is made into a one-dimensional grid.
 *  The velocties are counted first, then the points in direction 1, direction 2, direction 3.
 *  Then the node values are colse to the cell values.
 *  After these two the populations from the next block are placed to this block.
 *  It can be assumed that each block is like an independent grid with Node and Cell values
*/    
    unsigned long long int getIndex (const int i1,const int i2,const int i3,  const int block,const int nodeType,const int dv) const
    { 
        return  ((((block*2+nodeType)*n3+i3)*n2 + i2)*n1 + i1)*numField + dv;
    }
    
    dataType  operator()(const int i1,const int i2,const int i3,  const int block,const int nodeType,const int dv) const  { return data[ getIndex(i1,i2,i3,block,nodeType,dv)];} 
    dataType& operator()(const int i1,const int i2,const int i3,  const int block,const int nodeType,const int dv)        { return data[ getIndex(i1,i2,i3,block,nodeType,dv)];} 
    
    dataType value(const unsigned long long int i) const { return data[i];} 
    dataType& value(const unsigned long long int i)      { return data[i];} 
    
    
///*! \brief initialises integers to the grid
// * Each population is given it's index value which is unique.
//*/        
//    void initializeIntegers()
//    {
//        for ( int block=0; block<numBlock; block++)
//            for ( int i3=ndB3; i3<=ndE3; i3++)
//                for ( int i2=ndB2; i2<=ndE2; i2++) 
//                    for ( int i1=ndB1; i1<=ndE1; i1++)   
//                        for ( int dv=0; dv<numfield; dv++)
//                        {
//                            data[ getIndex(i1,i2,i3,block,nodeTYPE::NODE,dv)]= getIndex(i1,i2,i3,block,nodeTYPE::NODE,dv);
//                            data[ getIndex(i1,i2,i3,block,nodeTYPE::CELL,dv)]= getIndex(i1,i2,i3,block,nodeTYPE::CELL,dv);
//                        }
//    }
    
    template < int N, int bL, typename dataType1>
    friend std::ostream& operator<< (std::ostream& os,const gridBCC3D<N, bL, dataType1> &thisVec);
    
    
    
    
public:
    int m1 /** Number of real points in direction 1*/, m2, m3;                
    int nB1/**The Beginning point of the real grid in direction 1*/, nB2, nB3;           
    int nE1/**The End point of the real grid in direction 1*/, nE2, nE3;             
    int ndB1/**The first dummy point or the beginning of the grid*/, ndB2, ndB3;             
    int ndE1/**The last  dummy point or the end of the grid*/, ndE2, ndE3;             
    int n1/** Total number of points in direction 1*/, n2, n3;                
    int numBlock;
    unsigned long long int sizeGrid;   /**<total number of computational grid points n1*n2*n3 */             
    unsigned long long int numField;
    unsigned long long int numElem;                   
    unsigned long long int actualGridSize;/**total number of real grid points m1*m2*m3 */  
    int node[numblock];
    int cell[numblock];
    dataType *data; 
};


#endif
