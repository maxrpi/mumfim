#ifndef AMSI_DOF_HOLDER_H_
#define AMSI_DOF_HOLDER_H_
#include <amsiMPI.h>
namespace amsi
{
  class DofHolder
  {
  protected:
    int dof;
    bool lcl;
  public:
    DofHolder()
      : dof()
      , lcl(false)
    {
      int sz = 0;
      int rnk = -1;
      MPI_Comm_rank(AMSI_COMM_SCALE,&rnk);
      MPI_Comm_size(AMSI_COMM_SCALE,&sz);
      if(rnk == sz - 1)
        lcl = true;
    }
    bool isLocal() {return lcl;}
    void setLocal(bool l) {lcl = l;}
    int getDof() {return dof;}
    void setDof(int d) {dof = d;}
  };
}
#endif
