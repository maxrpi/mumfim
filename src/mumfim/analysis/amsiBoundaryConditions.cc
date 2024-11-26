#include "amsiBoundaryConditions.h"
#include <apfShape.h>
#include "amsiNeumannIntegratorsMT.h"
namespace amsi {
  // returns the component of a vector that a scalar represents
  // for example if we get x,y,z as scalar name, they should
  // actually assign to the vector field
  static int ScalarComponent(const std::string& name)
  {
    if (name == "x") {
      return 0;
    }
    if (name == "y") {
      return 1;
    }
    if (name == "z") {
      return 2;
    }
    return -1;
  }
  static int SetScalarValue(const std::string& name, apf::Field* field,
                            apf::Field* delta_field, apf::Numbering* nm,
                            apf::MeshEntity* ent, int nd,
                            std::vector<double>& data,
                            apf::Numbering* already_set, double val)
  {
    int cmp = ScalarComponent(name);
    if (cmp < 1) {
      cmp = 0;
    }
    // if the node has already been set by lower order BC don't
    // set it again
    if (apf::isNumbered(already_set, ent, nd, cmp)) {
      return 0;
    }
    apf::number(already_set, ent, nd, cmp, 1);
    data.resize(apf::countComponents(field));
    apf::getComponents(field, ent, nd, data.data());
    double delta_val = val - data[cmp];
    data[cmp] = val;
    apf::setComponents(field, ent, nd, data.data());
    // set the delta fields
    if (delta_field) {
      apf::getComponents(delta_field, ent, nd, data.data());
      data[cmp] = delta_val;
      apf::setComponents(delta_field, ent, nd, data.data());
    }
    apf::fix(nm, ent, nd, cmp, true);
    return 1;
  }
  static apf::Vector3 GetPosition(apf::Mesh* msh, apf::MeshEntity* ent, int nd)
  {
    apf::Vector3 pt;
    msh->getPoint(ent, nd, pt);
    // map local coord to global coord
    apf::Vector3 xyz;
    apf::MeshElement* melmt = apf::createMeshElement(msh, ent);
    apf::mapLocalToGlobal(melmt, pt, xyz);
    apf::destroyMeshElement(melmt);
    return xyz;
  }
  template <typename MT, typename... Args>
  int SetVectorValue(apf::Field* field, apf::Field* delta_field,
                     apf::Numbering* nm, apf::MeshEntity* ent, int nd,
                     apf::Numbering* already_set, const MT& mt, Args... args)
  {
    int num_fixed = 0;
    size_t num_comp = apf::countComponents(field);
    if (mt.size() != num_comp) {
      std::cerr << "Vector BC with wrong number of components\n";
      exit(1);
    }
    std::vector<double> data(num_comp);
    std::vector<double> delta_data(num_comp);
    apf::getComponents(field, ent, nd, data.data());
    if (delta_field != nullptr) {
      apf::getComponents(delta_field, ent, nd, delta_data.data());
    }
    for (size_t i = 0; i < num_comp; ++i) {
      // if the node has already been set by lower order BC don't
      // set it again
      if (!apf::isNumbered(already_set, ent, nd, i)) {
        apf::number(already_set, ent, nd, i, 1);
        double val = mt(i, args...);
        delta_data[i] = val - data[i];
        data[i] = val;
        apf::fix(nm, ent, nd, i, true);
        num_fixed += 1;
      }
    }
    apf::setComponents(field, ent, nd, data.data());
    // set the delta fields
    if (delta_field != nullptr) {
      apf::setComponents(delta_field, ent, nd, delta_data.data());
    }
    return num_fixed;
  }
  class SetDirichletVisitor : public mt::MTVisitor {
    public:
    SetDirichletVisitor(const std::string& name, apf::Field* field,
                        apf::Field* delta_field, apf::Numbering* nm,
                        double time, apf::Numbering* already_set,
                        int& num_fixed)
        : field_(field)
        , delta_field_(delta_field)
        , nm_(nm)
        , time_(time)
        , already_set_(already_set)
        , name_(name)
        , ent_()
        , nd_()
        , num_fixed_(num_fixed)
    {
    }
    void update(apf::MeshEntity* ent, int nd)
    {
      ent_ = ent;
      nd_ = nd;
    }
    void visit(mt::BoolMT& mt) override
    {
      num_fixed_ += SetScalarValue(name_, field_, delta_field_, nm_, ent_, nd_,
                                   data_, already_set_, mt());
    }
    void visit(mt::ScalarMT& mt) override
    {
      num_fixed_ += SetScalarValue(name_, field_, delta_field_, nm_, ent_, nd_,
                                   data_, already_set_, mt());
    };
    void visit(mt::IntMT& mt) override
    {
      num_fixed_ += SetScalarValue(name_, field_, delta_field_, nm_, ent_, nd_,
                                   data_, already_set_, mt());
    }
    void visit(mt::VectorMT& mt) override
    {
      num_fixed_ += SetVectorValue(field_, delta_field_, nm_, ent_, nd_,
                                   already_set_, mt);
    }
    void visit(mt::StringMT&) override
    {
      std::cerr << "Cannot apply a string boundary condition\n";
      exit(1);
    }
    void visit(mt::MatrixMT&) override
    {
      std::cerr << "Cannot apply a matrix boundary condition\n";
      exit(1);
    }
    // equation types
    // 4 parameters (space and time)
    void visit(mt::BoolFunctionMT<4>& mt) override
    {
      auto xyz = GetPosition(apf::getMesh(field_), ent_, nd_);
      num_fixed_ +=
          SetScalarValue(name_, field_, delta_field_, nm_, ent_, nd_, data_,
                         already_set_, mt(time_, xyz[0], xyz[1], xyz[2]));
    }
    void visit(mt::ScalarFunctionMT<4>& mt) override
    {
      auto xyz = GetPosition(apf::getMesh(field_), ent_, nd_);
      num_fixed_ +=
          SetScalarValue(name_, field_, delta_field_, nm_, ent_, nd_, data_,
                         already_set_, mt(time_, xyz[0], xyz[1], xyz[2]));
    }
    void visit(mt::IntFunctionMT<4>& mt) override
    {
      auto xyz = GetPosition(apf::getMesh(field_), ent_, nd_);
      num_fixed_ +=
          SetScalarValue(name_, field_, delta_field_, nm_, ent_, nd_, data_,
                         already_set_, mt(time_, xyz[0], xyz[1], xyz[2]));
    }
    void visit(mt::VectorFunctionMT<4>& mt) override
    {
      auto xyz = GetPosition(apf::getMesh(field_), ent_, nd_);
      num_fixed_ +=
          SetVectorValue(field_, delta_field_, nm_, ent_, nd_, already_set_, mt,
                         time_, xyz[0], xyz[1], xyz[2]);
    }
    void visit(mt::StringFunctionMT<4>&) override
    {
      std::cerr << "Cannot apply a string boundary condition\n";
      exit(1);
    }
    void visit(mt::MatrixFunctionMT<4>&) override
    {
      std::cerr << "Cannot apply a matrix boundary condition\n";
      exit(1);
    }
    // 3 parameters (3D space)
    void visit(mt::BoolFunctionMT<3>& mt) override
    {
      auto xyz = GetPosition(apf::getMesh(field_), ent_, nd_);
      num_fixed_ +=
          SetScalarValue(name_, field_, delta_field_, nm_, ent_, nd_, data_,
                         already_set_, mt(xyz[0], xyz[1], xyz[2]));
    }
    void visit(mt::ScalarFunctionMT<3>& mt) override
    {
      auto xyz = GetPosition(apf::getMesh(field_), ent_, nd_);
      num_fixed_ +=
          SetScalarValue(name_, field_, delta_field_, nm_, ent_, nd_, data_,
                         already_set_, mt(xyz[0], xyz[1], xyz[2]));
    }
    void visit(mt::IntFunctionMT<3>& mt) override
    {
      auto xyz = GetPosition(apf::getMesh(field_), ent_, nd_);
      num_fixed_ +=
          SetScalarValue(name_, field_, delta_field_, nm_, ent_, nd_, data_,
                         already_set_, mt(xyz[0], xyz[1], xyz[2]));
    }
    void visit(mt::VectorFunctionMT<3>& mt) override
    {
      auto xyz = GetPosition(apf::getMesh(field_), ent_, nd_);
      num_fixed_ += SetVectorValue(field_, delta_field_, nm_, ent_, nd_,
                                   already_set_, mt, xyz[0], xyz[1], xyz[2]);
    }
    void visit(mt::StringFunctionMT<3>&) override
    {
      std::cerr << "Cannot apply a string boundary condition\n";
      exit(1);
    }
    void visit(mt::MatrixFunctionMT<3>&) override
    {
      std::cerr << "Cannot apply a matrix boundary condition\n";
      exit(1);
    }
    // 2 parameters (2D space)
    void visit(mt::BoolFunctionMT<2>& mt) override
    {
      auto xyz = GetPosition(apf::getMesh(field_), ent_, nd_);
      num_fixed_ += SetScalarValue(name_, field_, delta_field_, nm_, ent_, nd_,
                                   data_, already_set_, mt(xyz[0], xyz[1]));
    }
    void visit(mt::ScalarFunctionMT<2>& mt) override
    {
      auto xyz = GetPosition(apf::getMesh(field_), ent_, nd_);
      num_fixed_ += SetScalarValue(name_, field_, delta_field_, nm_, ent_, nd_,
                                   data_, already_set_, mt(xyz[0], xyz[1]));
    }
    void visit(mt::IntFunctionMT<2>& mt) override
    {
      auto xyz = GetPosition(apf::getMesh(field_), ent_, nd_);
      num_fixed_ += SetScalarValue(name_, field_, delta_field_, nm_, ent_, nd_,
                                   data_, already_set_, mt(xyz[0], xyz[1]));
    }
    void visit(mt::VectorFunctionMT<2>& mt) override
    {
      auto xyz = GetPosition(apf::getMesh(field_), ent_, nd_);
      num_fixed_ += SetVectorValue(field_, delta_field_, nm_, ent_, nd_,
                                   already_set_, mt, xyz[0], xyz[1]);
    }
    void visit(mt::StringFunctionMT<2>&) override
    {
      std::cerr << "Cannot apply a string boundary condition\n";
      exit(1);
    }
    void visit(mt::MatrixFunctionMT<2>&) override
    {
      std::cerr << "Cannot apply a matrix boundary condition\n";
      exit(1);
    }
    // 1 parameters (time)
    void visit(mt::BoolFunctionMT<1>& mt) override
    {
      num_fixed_ += SetScalarValue(name_, field_, delta_field_, nm_, ent_, nd_,
                                   data_, already_set_, mt(time_));
    }
    void visit(mt::ScalarFunctionMT<1>& mt) override
    {
      num_fixed_ += SetScalarValue(name_, field_, delta_field_, nm_, ent_, nd_,
                                   data_, already_set_, mt(time_));
    }
    void visit(mt::IntFunctionMT<1>& mt) override
    {
      num_fixed_ += SetScalarValue(name_, field_, delta_field_, nm_, ent_, nd_,
                                   data_, already_set_, mt(time_));
    }
    void visit(mt::VectorFunctionMT<1>& mt) override
    {
      num_fixed_ += SetVectorValue(field_, delta_field_, nm_, ent_, nd_,
                                   already_set_, mt, time_);
    }
    void visit(mt::StringFunctionMT<1>&) override
    {
      std::cerr << "Cannot apply a string boundary condition\n";
      exit(1);
    }
    void visit(mt::MatrixFunctionMT<1>&) override
    {
      std::cerr << "Cannot apply a matrix boundary condition\n";
      exit(1);
    }

    private:
    // set once
    apf::Field* field_;
    apf::Field* delta_field_;
    apf::Numbering* nm_;
    double time_;
    apf::Numbering* already_set_;
    // set before each call
    const std::string& name_;
    apf::MeshEntity* ent_;
    int nd_;
    // don't set
    std::vector<double> data_;
    int& num_fixed_;
  };
  static int setDirichletValue(mt::IModelTrait* bc, const std::string& name,
                               apf::Numbering* nm, apf::Field* field,
                               apf::Field* delta_field, apf::MeshEntity* ent,
                               double t, apf::Numbering* already_set)
  {
    int num_fixed = 0;
    SetDirichletVisitor visitor(name, field, delta_field, nm, t, already_set,
                                num_fixed);
    assert(field != nullptr);
    auto* field_shape = apf::getShape(field);
    int num_nds = field_shape->countNodesOn(apf::getMesh(field)->getType(ent));
    for (int nd = 0; nd < num_nds; ++nd) {
      visitor.update(ent, nd);
      bc->accept(visitor);
    }
    return num_fixed;
  }
  int applyDirichletBCs(apf::Numbering* nm, apf::Field* delta_field,
                        const mt::AssociatedModelTraits<mt::DimIdGeometry>& bcs,
                        const std::vector<DirichletBCEntry>& bc_paths, double t)
  {
    int num_fixed = 0;
    auto* mesh = apf::getMesh(nm);
    auto* field = apf::getField(nm);
    // a numbering so we can keep track of if we already set the components of
    // the field
    auto* already_set =
        apf::createNumbering(mesh, "bc_already_set", apf::getShape(field),
                             apf::countComponents(field));
    for (int dimension = 0; dimension < mesh->getDimension(); ++dimension) {
      auto* it = mesh->begin(dimension);
      apf::MeshEntity* e = nullptr;
      while ((e = mesh->iterate(it))) {
        if (!mesh->isOwned(e)) {
          continue;
        }
        auto* geom_entity = mesh->toModel(e);
        int model_dim = mesh->getModelType(geom_entity);
        int model_tag = mesh->getModelTag(geom_entity);
        const auto* traits = bcs.Find({model_dim, model_tag});
        if (traits == nullptr) {
          continue;
        }
        for (const auto& path : bc_paths) {
          const mt::AssociatedCategoryNode* nd = traits;
          for (auto& category : path.categories) {
            nd = nd->FindCategoryByType(category);
            if (nd == nullptr) {
              break;
            }
          }
          if (nd == nullptr) {
            continue;
          }
          const auto& bc_name = path.mt_name;
          auto* bc = mt::GetCategoryModelTraitByType(nd, bc_name);
          if (bc != nullptr) {
            num_fixed +=
                setDirichletValue(const_cast<mt::IModelTrait*>(bc), bc_name, nm,
                                  field, delta_field, e, t, already_set);
            // apply the bcs to any entities classified on the model entity
            // of a lower dimension. This makes sure we apply the dirichlet bc
            // to all entities on the closure of the geometry
            for (int lower_dim = 0; lower_dim < dimension; ++lower_dim) {
              apf::Adjacent adjacent;
              mesh->getAdjacent(e, lower_dim, adjacent);
              for (auto& adj : adjacent) {
                num_fixed += setDirichletValue(const_cast<mt::IModelTrait*>(bc),
                                               bc_name, nm, field, delta_field,
                                               adj, t, already_set);
              }
            }
          }
        }
      }
      mesh->end(it);
    }
    apf::destroyNumbering(already_set);
    return num_fixed;
  }
  void applyNeumannBCs(LAS* las, apf::Numbering* nm,
                       const mt::AssociatedModelTraits<mt::DimIdGeometry>& bcs,
                       const std::vector<NeumannBCEntry>& bc_paths, double t)
  {
    std::unordered_map<NeumannBCType, std::unique_ptr<NeumannIntegratorMT>>
        integrators;
    apf::Field* fld = apf::getField(nm);
    apf::Mesh* mesh = apf::getMesh(nm);
    for (int dimension = 0; dimension < mesh->getDimension(); ++dimension) {
      auto* it = mesh->begin(dimension);
      apf::MeshEntity* e = nullptr;
      while ((e = mesh->iterate(it))) {
        if (!mesh->isOwned(e)) {
          continue;
        }
        auto* geom_entity = mesh->toModel(e);
        int model_dim = mesh->getModelType(geom_entity);
        // only integrate when the model dimension is of the same dimension
        // as the mesh entity dimension
        if (model_dim != dimension) {
          continue;
        }
        int model_tag = mesh->getModelTag(geom_entity);
        const auto* traits = bcs.Find({model_dim, model_tag});
        if (traits == nullptr) {
          continue;
        }
        for (const auto& path : bc_paths) {
          const mt::AssociatedCategoryNode* nd = traits;
          for (auto& category : path.categories) {
            nd = nd->FindCategoryByType(category);
            if (nd == nullptr) {
              break;
            }
          }
          if (nd == nullptr) {
            continue;
          }
          const auto& bc_name = path.mt_name;
          auto* bc = mt::GetCategoryModelTraitByType(nd, bc_name);
          if (bc != nullptr) {
            NeumannIntegratorMT* integrator;
            auto* mnt = apf::createMeshElement(mesh, e);
            auto integrator_type = path.mt_type;
            auto rslt = integrators.find(integrator_type);
            if (rslt == integrators.end()) {
              auto r = integrators.emplace(
                  integrator_type,
                  createNeumannIntegrator(las, fld, bc, 1, t, integrator_type));
              integrator = r.first->second.get();
            }
            else {
              integrator = rslt->second.get();
            }
            integrator->process(mnt);
            apf::NewArray<int> dofs;
            apf::getElementNumbers(nm, e, dofs);
            las->AddToVector(integrator->getnedofs(), &dofs[0],
                             &integrator->getFe()[0]);
            las->AddToMatrix(integrator->getnedofs(), &dofs[0],
                             integrator->getnedofs(), &dofs[0],
                             &(integrator->getKe()(0,0)) );
            apf::destroyMeshElement(mnt);
          }
        }
      }
    }
  }
}  // namespace amsi
