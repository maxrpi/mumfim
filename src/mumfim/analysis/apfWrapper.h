#ifndef APF_WRAPPER_H_
#define APF_WRAPPER_H_
#include <apf.h>
namespace amsi
{
  /**
   * Get the values (scalar, vector, or matrix) of all DOFs
   *  effecting the element in canonical order.
   */
  void getDofValues(apf::Field * fld, apf::Element * elmt, double * vls);
  /**
   * Get the mean value of the field values (scalar, vector, or matrix)
   *  on all nodes on all mesh entities in range [bgn,nd).
   * @note Collective on AMSI_COMM_SCALE
   */
  template <class I>
    void getAvgFieldValue(apf::Field * fld, I bgn, I nd, double * avg);
  class TagIterator
  {
  private:
  public:
    bool operator++();
    bool operator==(const TagIterator & o);
  };
};
#include "apfWrapper_impl.h"
#endif
