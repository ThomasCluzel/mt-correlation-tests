//Added from http://wiki.cs.princeton.edu/index.php/TNT
#include <cassert>
#include "tnt_array1d.h"
#include "tnt_array2d.h"
#include "jama_lu.h"

template<class T>
TNT::Array2D<T> invert(const TNT::Array2D<T> &M)
{
	assert(M.dim1() == M.dim2()); // square matrices only please

	// solve for inverse with LU decomposition
	JAMA::LU<T> lu(M);

	// create identity matrix
	TNT::Array2D<T> id(M.dim1(), M.dim2(), (T)0);
	for (int i = 0; i < M.dim1(); i++) id[i][i] = 1;

        // solves A * A_inv = Identity
	return lu.solve(id);
}

template<class T>
TNT::Array2D<T> transpose(const TNT::Array2D<T> &M)
{
	TNT::Array2D<T> tran(M.dim2(), M.dim1() );
	for(int r=0; r<M.dim1(); ++r)
		for(int c=0; c<M.dim2(); ++c)
			tran[c][r] = M[r][c];
	return tran;
}
