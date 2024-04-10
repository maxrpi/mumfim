#include "amsiFields.h"
#include <pystring.h>
#include <algorithm>
#include <cassert>
namespace amsi
{
  /*
  static char const * fld_units[] =
  {
    "untyped",
    "displacement"
  };
  static char const * fld_tps[] =
  {
    "full",
    "delta"
  };
  */
  char const * fieldUnitString(FieldUnit nt)
  {
    assert(nt >= 0 && nt < FieldUnit::num_field_units);
    return FieldUnits[nt];
  }
  char const * fieldTypeString(FieldType tp)
  {
    assert(tp >= 0 && tp < FieldType::num_field_types);
    return FieldTypes[tp];
  }
  int decodeFieldUnit(const std::string & fld_nm)
  {
    std::vector<std::string> tks;
    pystring::split(fld_nm,tks,std::string("_"));
    int tp = std::distance(&FieldUnits[0],std::find(&FieldUnits[0],&FieldUnits[FieldUnit::num_field_units],tks[0]));
    assert(tp < FieldUnit::num_field_units);
    return tp;
  }
  int decodeFieldType(const std::string & fld_nm)
  {
    std::vector<std::string> tks;
    pystring::split(fld_nm,tks,std::string("_"));
    int tp = std::distance(&FieldTypes[0],std::find(&FieldTypes[0],&FieldTypes[FieldType::num_field_types],tks[2]));
    assert(tp < FieldType::num_field_types);
    return tp;
  }
  std::string decodeFieldName(apf::Field * fld)
  {
    std::vector<std::string> tks;
    pystring::split(apf::getName(fld),tks,std::string("_"));
    assert(tks.size() == 5);
    return tks[4];
  }
  std::string composeFieldName(int nt, int tp, const std::string & nm)
  {
    std::vector<std::string> nm_vc(3);
    nm_vc[0] = std::string(fieldUnitString(static_cast<FieldUnit>(nt)));
    nm_vc[1] = std::string(fieldTypeString(static_cast<FieldType>(tp)));
    nm_vc[2] = nm;
    return pystring::join("_",nm_vc);
   }
}
