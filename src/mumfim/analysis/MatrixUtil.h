#ifndef MATRIXUTIL_H_
#define MATRIXUTIL_H_
namespace amsi
{
  void direct_solver(double ** matrix,
                     double * vector,
                     double * solution,
                     int size);
  void free_matrix( double ** matrix);
  double ** allocate_matrix(int size);
}
#endif
