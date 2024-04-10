#ifndef AMSI_FIELDS_H_
#define AMSI_FIELDS_H_
#include "amsiLAS.h"
#include <amsiEnumOps.h>
#include <apf.h>
#include <vector>
namespace amsi
{
  #define FIELD_UNITS(OP) OP(unitless), OP(displacement), OP(num_field_units)
  #define FIELD_TYPES(OP) OP(full), OP(delta), OP(num_field_types)
  enum FieldUnit{FIELD_UNITS(MAKE_ENUM_OP)};
  const char * const FieldUnits[] = {FIELD_UNITS(MAKE_STRING_OP)};
  enum FieldType{FIELD_TYPES(MAKE_ENUM_OP)};
  const char * const FieldTypes[] = {FIELD_TYPES(MAKE_STRING_OP)};
  // note: modifications here MUST be reflected in amsiFields.cc::fld_strs
  char const * fieldUnitString(FieldUnit nt);
  char const * fieldTypeString(FieldType tp);
  int decodeFieldUnit(const std::string & fld_nm);
  int decodeFieldType(const std::string & fld_nm);
  std::string decodeFieldName(apf::Field * fld);
  std::string composeFieldName(int nt, int tp, const std::string & nm);
}
#include "amsiFields_impl.h"
#endif
