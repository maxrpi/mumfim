#ifndef APF_MEASURE_IMPL_H_
#define APF_MEASURE_IMPL_H_
namespace amsi
{
  template <typename I>
    double measureDisplacedModelEntities(I b, I e, apf::Field * u)
  {
    double v = 0.0;
    for(auto i = b; i != e; ++i)
      v += measureDisplacedModelEntity(*i,u);
    return v;
  }
}
#endif
