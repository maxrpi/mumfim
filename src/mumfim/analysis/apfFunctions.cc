#include "apfFunctions.h"
#include "ElementalSystem.h"
#include <apfDynamicMatrix.h>
#include <apfGeometry.h>
#include <apfMesh.h>
#include <apfShape.h>
#include <cassert>
namespace amsi
{
  void writePvdFile(const std::string & col_fnm, const std::string & msh_prfx, int sz)
  {
    std::string pvd(fs->getResultsDir() + col_fnm);
    std::fstream pvdf(pvd.c_str(), std::ios::out);
    pvdf << "<VTKFile type=\"Collection\" version=\"0.1\">" << std::endl;
    pvdf << "  <Collection>" << std::endl;
    for(int ii = 0; ii < sz; ii++)
    {
      std::ostringstream oss;
      oss << msh_prfx << ii+1;
      std::string vtu = oss.str();
      pvdf << "    <DataSet timestep=\"" << ii << "\" group=\"\" ";
      pvdf << "part=\"0\" file=\"" << vtu << "/" << vtu;
      pvdf << ".pvtu\"/>" << std::endl;
    }
    pvdf << "  </Collection>" << std::endl;
    pvdf << "</VTKFile>" << std::endl;
  }
  void writePvdFile(const std::string & col_fnm, const std::vector<PvdData> & pvd_data) {
     std::stringstream pvd;
     if(fs)
       pvd<<fs->getResultsDir();
     pvd <<col_fnm;
    //std::string pvd(fs->getResultsDir() + col_fnm);
    std::fstream pvdf(pvd.str().c_str(), std::ios::out);
    pvdf << "<VTKFile type=\"Collection\" version=\"0.1\">" << std::endl;
    pvdf << "  <Collection>" << std::endl;
    for(std::size_t ii = 0; ii < pvd_data.size(); ii++)
    {
      pvdf << "    <DataSet timestep=\"" << pvd_data[ii].timestep << "\" group=\"\" ";
      pvdf << "part=\""<<pvd_data[ii].part<<"\" file=\"" << pvd_data[ii].filename << "/" << pvd_data[ii].filename;
      pvdf << ".pvtu\"/>" << std::endl;
    }
    pvdf << "  </Collection>" << std::endl;
    pvdf << "</VTKFile>" << std::endl;
  }
  bool isFixed(int n)
  {
    return n == -2;
  }
  void faceNormal(apf::Mesh * msh, apf::MeshEntity * fc, apf::Vector3 & n)
  {
    assert(msh);
    assert(fc);
    n.zero();
    apf::Field * crds = msh->getCoordinateField();
    apf::Adjacent vs;
    msh->getAdjacent(fc,0,vs);
    apf::Vector3 vcrds[3];
    apf::getVector(crds,vs[0],0,vcrds[0]);
    apf::getVector(crds,vs[1],0,vcrds[1]);
    apf::getVector(crds,vs[2],0,vcrds[2]);
    apf::Plane p = apf::Plane::fromPoints(vcrds[0],vcrds[1],vcrds[2]);
    n = p.normal;
  }
  void vertexNormal(apf::Mesh * msh, apf::MeshEntity * vrt, apf::Vector3 & n)
  {
    assert(msh);
    assert(vrt);
    n.zero();
    apf::Adjacent adj;
    msh->getAdjacent(vrt,2,adj);
    APF_ITERATE(apf::Adjacent,adj,fc)
    {
      apf::Vector3 fn;
      faceNormal(msh,*fc,fn);
      n += fn;
    }
    n = n.normalize();
  }
  void displaceMesh(apf::Field * displacement_field)
  {
    apf::Mesh * mesh = apf::getMesh(displacement_field);
    apf::Field * apf_coord_field = mesh->getCoordinateField();
    for(int dim = 0; dim <= mesh->getDimension(); dim++)
    {
      apf::FieldShape * fs = apf::getShape(displacement_field);
      if(!fs->hasNodesIn(dim))
        break;
      apf::MeshIterator * it = mesh->begin(dim);
      while(apf::MeshEntity * me = mesh->iterate(it))
      {
        int num_nodes = fs->countNodesOn(mesh->getType(me));
        for(int jj = 0; jj < num_nodes; jj++)
        {
          apf::Vector3 disp, coord;
          apf::getVector(displacement_field,me,jj,disp);
          apf::getVector(apf_coord_field,me,jj,coord);
          apf::Vector3 new_coord = coord + disp;
          apf::setVector(apf_coord_field,me,jj,new_coord);
        }
      }
      mesh->end(it);
    }
  }
  void printNumbering(std::ostream & out, apf::Numbering * numbering)
  {
    apf::Mesh * mesh = apf::getMesh(numbering);
    apf::Field * field = apf::getField(numbering);
    int num_components = apf::countComponents(field);
    apf::FieldShape * field_shape = apf::getShape(field);
    int analysis_dim = mesh->getDimension();
    for(int dim = 0; dim <= analysis_dim; dim++)
    {
      out << dim << " : " << std::endl;
      apf::MeshIterator * it = mesh->begin(dim);
      while(apf::MeshEntity * me = mesh->iterate(it))
      {
        if(mesh->isOwned(me))
        {
          int num_nodes = field_shape->countNodesOn(mesh->getType(me));
          for(int ii = 0; ii < num_nodes; ii++)
          {
            apf::Vector3 coord;
            mesh->getPoint(me,ii,coord);
            for(int jj = 0; jj < num_components; jj++)
            {
              if(apf::isNumbered(numbering,me,ii,jj))
                out << apf::getNumber(numbering,me,ii,jj) << "\t";
              else
                out << "F\t";
              out << coord[jj] << std::endl;
            }
          }
        }
      }
      mesh->end(it);
    }
  }
  void SymmMatrixToVoigtVector(const apf::Matrix3x3 & mat,
                               apf::Vector<6> & vec)
  {
    vec[0] = mat[0][0];
    vec[1] = mat[1][1];
    vec[2] = mat[2][2];
    vec[3] = mat[0][1];
    vec[4] = mat[1][2];
    vec[5] = mat[0][2];
  }
  void VoigtVectorToSymmMatrix(const apf::Vector<6> & vec,
                               apf::Matrix3x3 & mat)
  {
    mat[0][0] = vec[0];
    mat[1][1] = vec[1];
    mat[2][2] = vec[2];
    mat[0][1] = mat[1][0] = vec[3];
    mat[1][2] = mat[2][1] = vec[4];
    mat[0][2] = mat[2][0] = vec[5];
  }
  void VoigtVectorToSymmMatrix(const apf::DynamicVector & vec,
                               apf::Matrix3x3 & mat)
  {
    mat[0][0] = vec[0];
    mat[1][1] = vec[1];
    mat[2][2] = vec[2];
    mat[0][1] = mat[1][0] = vec[3];
    mat[1][2] = mat[2][1] = vec[4];
    mat[0][2] = mat[2][0] = vec[5];
  }
}
