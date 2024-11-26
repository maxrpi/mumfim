#ifndef AMSI_AMSINEUMANNINTEGRATORSMT_H
#define AMSI_AMSINEUMANNINTEGRATORSMT_H
#include <apf.h>
#include <apfDynamicVector.h>
#include <apfDynamicMatrix.h>
#include <model_traits/ModelTrait.h>
#include "amsiLAS.h"
#include "apfFunctions.h"
namespace amsi {
  enum class NeumannBCType { pressure, traction, robin};

  class MTEvaluator : public mt::MTVisitor {
    template <typename MT, typename... Args>
    void SetVector(const MT& mt, Args... args)
    {
      vals_->reserve(mt.size());
      for (size_t i = 0; i < mt.size(); ++i) {
        vals_->push_back(mt(i, args...));
      }
    }
    public:
    MTEvaluator(const mt::IModelTrait* mt) : vals_(nullptr), model_trait_(mt) {}
    void operator()(double t, double x, double y, double z,
                    std::vector<double>& out_vls)
    {
      vals_ = &out_vls;
      vals_->clear();
      if (model_trait_ != nullptr) {
        x_ = x;
        y_ = y;
        z_ = z;
        t_ = t;
        // the visitor doesn't modify the mt, but there is currently not ConstVisitor
        const_cast<mt::IModelTrait*>(model_trait_)->accept(*this);
      }
    }
    void visit(mt::BoolMT& mt) final { vals_->push_back(mt()); }
    void visit(mt::ScalarMT& mt) final { vals_->push_back(mt()); };
    void visit(mt::IntMT& mt) final { vals_->push_back(mt()); }
    void visit(mt::VectorMT& mt) final { SetVector(mt); }
    void visit(mt::StringMT&) final
    {
      std::cerr << "Cannot apply a string boundary condition\n";
      exit(1);
    }
    void visit(mt::MatrixMT&) final
    {
      std::cerr << "Cannot apply a matrix boundary condition\n";
      exit(1);
    }
    // equation types
    // 4 parameters (space and time)
    void visit(mt::BoolFunctionMT<4>& mt) final
    {
      vals_->push_back(mt(t_, x_, y_, z_));
    }
    void visit(mt::ScalarFunctionMT<4>& mt) final
    {
      vals_->push_back(mt(t_, x_, y_, z_));
    }
    void visit(mt::IntFunctionMT<4>& mt) final
    {
      vals_->push_back(mt(t_, x_, y_, z_));
    }
    void visit(mt::VectorFunctionMT<4>& mt) final
    {
      SetVector(mt, t_, x_, y_, z_);
    }
    void visit(mt::StringFunctionMT<4>&) final
    {
      std::cerr << "Cannot apply a string boundary condition\n";
      exit(1);
    }
    void visit(mt::MatrixFunctionMT<4>&) final
    {
      std::cerr << "Cannot apply a matrix boundary condition\n";
      exit(1);
    }
    // 3 parameters (3D space)
    void visit(mt::BoolFunctionMT<3>& mt) final
    {
      vals_->push_back(mt(x_, y_, z_));
    }
    void visit(mt::ScalarFunctionMT<3>& mt) final
    {
      vals_->push_back(mt(x_, y_, z_));
    }
    void visit(mt::IntFunctionMT<3>& mt) final
    {
      vals_->push_back(mt(x_, y_, z_));
    }
    void visit(mt::VectorFunctionMT<3>& mt) final { SetVector(mt, x_, y_, z_); }
    void visit(mt::StringFunctionMT<3>&) final
    {
      std::cerr << "Cannot apply a string boundary condition\n";
      exit(1);
    }
    void visit(mt::MatrixFunctionMT<3>&) final
    {
      std::cerr << "Cannot apply a matrix boundary condition\n";
      exit(1);
    }
    // 2 parameters (2D space)
    void visit(mt::BoolFunctionMT<2>& mt) final
    {
      vals_->push_back(mt(x_, y_));
    }
    void visit(mt::ScalarFunctionMT<2>& mt) final
    {
      vals_->push_back(mt(x_, y_));
    }
    void visit(mt::IntFunctionMT<2>& mt) final { vals_->push_back(mt(x_, y_)); }
    void visit(mt::VectorFunctionMT<2>& mt) final { SetVector(mt, x_, y_); }
    void visit(mt::StringFunctionMT<2>&) final
    {
      std::cerr << "Cannot apply a string boundary condition\n";
      exit(1);
    }
    void visit(mt::MatrixFunctionMT<2>&) final
    {
      std::cerr << "Cannot apply a matrix boundary condition\n";
      exit(1);
    }
    // 1 parameters (time)
    void visit(mt::BoolFunctionMT<1>& mt) final { vals_->push_back(mt(t_)); }
    void visit(mt::ScalarFunctionMT<1>& mt) final { vals_->push_back(mt(t_)); }
    void visit(mt::IntFunctionMT<1>& mt) final { vals_->push_back(mt(t_)); }
    void visit(mt::VectorFunctionMT<1>& mt) final { SetVector(mt, t_); }
    void visit(mt::StringFunctionMT<1>&) final
    {
      std::cerr << "Cannot apply a string boundary condition\n";
      exit(1);
    }
    void visit(mt::MatrixFunctionMT<1>&) final
    {
      std::cerr << "Cannot apply a matrix boundary condition\n";
      exit(1);
    }
    private:
    std::vector<double>* vals_;
    const mt::IModelTrait* model_trait_;
    double x_{};
    double y_{};
    double z_{};
    double t_{};
  };

  class NeumannIntegratorMT : public apf::Integrator {
    protected:
    LAS* las;
    apf::DynamicMatrix Ke;
    apf::DynamicVector fe;
    apf::Field* fld;
    apf::MeshElement* me;
    apf::Element* e;
    int nedofs;
    int nenodes;
    int nfcmps;
    double tm;
    std::vector<int> dofs;
    MTEvaluator mt_evaluator;
    public:
    NeumannIntegratorMT(LAS* l, apf::Field* f, const mt::IModelTrait* model_trait,
                        int o, double t = 0.0);
    void setTime(double t);
    int getnedofs();
    apf::DynamicVector& getFe() { return fe; }
    apf::DynamicMatrix& getKe() { return Ke; }
    void inElement(apf::MeshElement* m) override;
    void outElement() override;
  };

  class SurfaceTractionMT : public NeumannIntegratorMT {
    public:
    SurfaceTractionMT(LAS* l, apf::Field* f, const mt::IModelTrait* mt, int o,
                      double t = 0.0);
    void atPoint(apf::Vector3 const& p, double w, double dV) final;
  };

  class PressureMT : public NeumannIntegratorMT {
    private:
    apf::Mesh* msh;
    apf::MeshEntity* ent;
    public:
    PressureMT(LAS* l, apf::Field* f, const mt::IModelTrait* mt, int o, double t);
    void inElement(apf::MeshElement* m) final;
    void atPoint(apf::Vector3 const& p, double w, double dV) final;
  };

  class RobinMT : public NeumannIntegratorMT {
    private:
    apf::Mesh* msh;
    apf::MeshEntity* ent;
    public:
    RobinMT(LAS* l, apf::Field* f, const mt::IModelTrait* mt, int o, double t);
    void inElement(apf::MeshElement* m) final;
    void atPoint(apf::Vector3 const& p, double w, double dV) final;
  };

  std::unique_ptr<NeumannIntegratorMT> createNeumannIntegrator(LAS* las,
                                                               apf::Field* fld,
                                                               const mt::IModelTrait* mt,
                                                               int o,
                                                               double t,
      NeumannBCType tp);
}  // namespace amsi
#endif  // AMSI_AMSINEUMANNINTEGRATORSMT_H
