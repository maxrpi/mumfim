#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
namespace amsi {
  using std::abs;
  using std::vector;
  using std::fill;
  using std::copy;
  using std::fprintf;
  const double TINY = 1.0e-20;
  template <typename IntArray>
  void ludcmp(double **a, int n, IntArray& indx, double& d);
  template <typename IntArray, typename DoubleArray>
  void lubksb(double **a, int n, IntArray& indx, DoubleArray& b);
  void invmat(double **a, double **y, int n);
  void nrerror(const char* error);
  double ** allocate_matrix(int n)
  {
    int i;
    double ** v = new double *[n];
    double *help = new double[n*n];
    for (i=0; i!= n; i++){
      v[i] = &help[i*n];
    }
    for(i=0;i<n*n;i++)help[i] = 0.0;
    return v;
  }
  void free_matrix(double **v)
  {
    if(v)
    {
      if(v[0])
        delete [] v[0];
      delete [] v;
    }
  }
  void direct_solver(double **a, double *y, double *sol, int n)
  {
    std::vector<int> indx(n);
    double d;
    ludcmp(a, n, indx, d);
    std::copy(y, y+n, sol);
    lubksb(a,n, indx,sol);
  }
  void invmat(double **a, double **y, int n)
  {
    vector<double> col(n);
    vector<int> indx(n);
    double d;
    ludcmp(a, n, indx, d);
    for(int j=0;j<n;j++){
      fill(col.begin(), col.end(), 0.0);
      col[j]=1.0;
      lubksb(a,n,indx, col);
      for(int i=0;i<n;i++) y[i][j]=col[i];
    }
  }
  template <typename IntArray>
  void ludcmp(double **a, int n, IntArray& indx, double& d)
  {
    int i, imax, j,k;
    imax = 0; // initialize imax to zero to quell compiler warnings
    double big,dum,sum,temp;
    std::vector<double> vv(n);
    d=1.0;
    for(i=0; i<n; i++) {
      big=0.0;
      for(j=0;j<n; j++)
        if((temp=abs(a[i][j])) > big) big = temp;
      if(big==0.0) nrerror("Singular matrix in routine LUDCMP");
      vv[i]=1.0/big;
    }
    for(j=0; j<n;j++){
      for(i=0;i<j;i++){
        sum=a[i][j];
        for(k=0;k<i;k++) sum -= a[i][k]*a[k][j];
        a[i][j] = sum;
      }
      big=0.0;
      for(i=j; i<n;i++){
        sum=a[i][j];
        for(k=0;k<j;k++)
          sum -= a[i][k]*a[k][j];
        a[i][j]=sum;
        if((dum=vv[i]*abs(sum)) >= big){
          big=dum;
          imax=i;
        }
      }
      if(j != imax){
        for(k=0; k <n; k++){
          dum=a[imax][k];
          a[imax][k]=a[j][k];
          a[j][k]=dum;
        }
        d = -(d);
        vv[imax]=vv[j];
      }
      indx[j]=imax;
      if( a[j][j] == 0.0) a[j][j] = TINY;
      if(j !=n){
        dum=1.0/(a[j][j]);
        for(i=j+1;i<n;i++) a[i][j] *= dum;
      }
    }
  }
  void nrerror(const char* error)
  {
    fprintf(stderr,"error\n");
    fprintf(stderr,"%s\n",error);
    std::exit(1);
  }
  template <typename IntArray, typename DoubleArray>
  void lubksb(double **a, int n, IntArray& indx, DoubleArray& b)
  {
    int i,ii=-1,ip,j;
    double sum;
    for(i=0; i< n; i++){
      ip=indx[i];
      sum=b[ip];
      b[ip]=b[i];
      if(ii != -1)
        for(j=ii; j <= i-1; j++) sum -= a[i][j]*b[j];
      else if (sum) ii=i;
      b[i]=sum;
    }
    for(i=n-1; i>=0; i--){
      sum=b[i];
      for(j=i+1;j<n;j++) sum -= a[i][j]*b[j];
      b[i]=sum/a[i][i];
    }
  }
}
