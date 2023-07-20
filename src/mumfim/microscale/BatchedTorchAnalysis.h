#ifndef MUMFIM_SRC_MUMFIM_MICROSCALE_BATCHEDTORCHANALYSIS_H
#define MUMFIM_SRC_MUMFIM_MICROSCALE_BATCHEDTORCHANALYSIS_H
#include <fmt/format.h>
#include <mumfim/exceptions.h>
#include <mumfim/microscale/BatchedRVEAnalysis.h>
#include <mumfim/microscale/ContinuumMechanics.h>
#include <mumfim/torch/TorchUtilities.h>
#include <torch/script.h>

#include <Kokkos_Core.hpp>

namespace mumfim
{

  template <typename Scalar = mumfim::Scalar,
            typename LocalOrdinal = mumfim::LocalOrdinal,
            typename ExeSpace = Kokkos::DefaultExecutionSpace>
  class BatchedTorchAnalysis
      : public BatchedRVEAnalysis<Scalar, LocalOrdinal, ExeSpace>
  {
    public:
    using exe_space = ExeSpace;
    using memory_space = typename ExeSpace::memory_space;

    BatchedTorchAnalysis(const std::string & model_path, int nrves)
        : BatchedRVEAnalysis<Scalar, LocalOrdinal, ExeSpace>(nrves),
          F_trial_("F_trial", nrves), stiffness_trial_("stiffness_trial", nrves),
          F_current_("F_current", nrves)
    {
      FillIdentity(F_current_);
      try
      {
        model_ = torch::jit::load(model_path);
        deformation_gradient_tensor_ = torch::empty({nrves, 3, 3});
        right_cauchy_green_.emplace_back(torch::empty({nrves, 3, 3}));
      }
      // forward torch exception as
      catch (const torch::Error & e)
      {
        throw material_error(fmt::format("mumfim TorchMaterial: {}", e.what()));
      }
    }

    bool run(Kokkos::DualView<Scalar * [3][3], ExeSpace> Fincrement,
             Kokkos::DualView<Scalar * [6], ExeSpace> sigma,
             bool update_coords = true) final
    {
      KOKKOS_ASSERT(Fincrement.extent(0) == this->num_rves_);
      KOKKOS_ASSERT(sigma.extent(0) == this->num_rves_);
      // 1. update the deformation gradient Finc@F
      Fincrement.template sync<ExeSpace>();
      UpdateDeformationGradient(Fincrement.template view<ExeSpace>(),
                                F_current_, F_trial_);
      auto right_cauchy_green = ComputeRightCauchyGreenDeformation(F_trial_);
      // convert kokkos arrays to torch tensors

      right_cauchy_green_[0].toTensor().set_requires_grad(false);
      ml::KokkosViewToTorchArray(right_cauchy_green,
                             right_cauchy_green_[0].toTensor());
      ml::KokkosViewToTorchArray(F_trial_, deformation_gradient_tensor_);
      right_cauchy_green_[0].toTensor().set_requires_grad(true);
      //  2. convert kokkos deformation gradients to torch tensors
      auto energy = model_.forward(right_cauchy_green_).toTensor();
      auto [cauchy_stress, stiffness] =
          ml::ComputeCauchyStressStiffness(energy, deformation_gradient_tensor_,
                                           right_cauchy_green_[0].toTensor());
      ml::TorchArrayToKokkosView(cauchy_stress, sigma.template view<ExeSpace>());
      ml::TorchArrayToKokkosView(stiffness, stiffness_trial_);
      sigma.template modify<ExeSpace>();

      return true;
    }

    // Accept the state from the last run. For this type of analysis, this
    // corresponds to updating the deformation current deformation gradient to
    // the previous state
    void accept() final { Kokkos::deep_copy(F_current_, F_trial_); }

    // computes the material stiffness tensor at the current deformation state
    // this can be computed by F@F@D@F@F where F is deformation tensor and
    // D= dPK2/dE.
    void computeMaterialStiffness(
        Kokkos::DualView<Scalar * [6][6], ExeSpace> C) final
    {
      KOKKOS_ASSERT(C.extent(0) == this->num_rves_);
      Kokkos::deep_copy(C.template view<ExeSpace>(), stiffness_trial_);
    }

    private:
    using LO = LocalOrdinal;
    // the tiral deformatio that was used in the last run. This is stored so
    // that it can update the current when the stat is accepted
    Kokkos::View<Scalar * [3][3], ExeSpace> F_trial_;
    // material stiffness tensor corresponding to the trial deformation state
    Kokkos::View<Scalar * [6][6], ExeSpace> stiffness_trial_;

    // current deformation gradient that is updated after each accept
    Kokkos::View<Scalar * [3][3], ExeSpace> F_current_;

    torch::jit::script::Module model_;

    torch::Tensor deformation_gradient_tensor_;
    // inputs are the right cauchy_green deformation_tensor
    std::vector<torch::jit::IValue> right_cauchy_green_;
  };

}  // namespace mumfim

#endif  // MUMFIM_SRC_MUMFIM_MICROSCALE_BATCHEDTORCHANALYSIS_H
