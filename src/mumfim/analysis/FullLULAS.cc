#include "FullLULAS.h"
#include "MatrixUtil.h"
namespace amsi
{
  FullLULAS::FullLULAS(int num_rows_in) :
    vector(num_rows_in, 0.0),
    solution(num_rows_in, 0.0),
    num_rows(num_rows_in)
  {
    matrix = allocate_matrix(num_rows);
    memset(&matrix[0][0],0,num_rows*num_rows*sizeof(double));
  }
  FullLULAS::~FullLULAS ( )
  {
    free_matrix(matrix);
  }
  void FullLULAS::Solve()
  {
    direct_solver(matrix, &vector[0], &solution[0], solution.size());
  }
}
