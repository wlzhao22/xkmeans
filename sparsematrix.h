#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

struct SparseMatrix
{
    float *data;//data
    unsigned int *col;//each col's first argument position
    unsigned int *index;//col index;
};


#endif
