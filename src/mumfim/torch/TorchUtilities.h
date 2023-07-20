#ifndef MUMFIM_SRC_MUMFIM_PYTORCH_TORCHUTILITIES_H
#define MUMFIM_SRC_MUMFIM_PYTORCH_TORCHUTILITIES_H
#include <apfDynamicMatrix.h>
#include <apfMatrix.h>
#include <torch/script.h>

#include <Kokkos_core.hpp>

namespace mumfim
{
  namespace ml
  {
    void ComputeRightCauchyGreen(const torch::Tensor & deformation_gradient,
                                 torch::Tensor & right_cauchy_green);

    [[nodiscard]] at::Tensor ComputeRightCauchyGreen(
        at::Tensor & deformation_gradient);

    [[nodiscard]] int VoigtIndex(int i, int j);

    [[nodiscard]] std::pair<int, int> Tensor4Index(int i);

    [[nodiscard]] torch::Tensor Tensor4ToVoigt(const torch::Tensor & tensor);

    [[nodiscard]] torch::Tensor VoigtToTensor4(const torch::Tensor & voigt);

    [[nodiscard]] torch::Tensor TotalLagrangianStiffnessToULStiffness(
        const torch::Tensor & material_stiffness,
        const torch::Tensor & deformation_gradient);

    [[nodiscard]] at::Tensor MandelProbeInv(int batch_size);

    [[nodiscard]] at::Tensor ProbeDirections(int batch_size);

    void MandelToVoigt(torch::Tensor & tensor);

    void VoigtToMandel(torch::Tensor & tensor);

    [[nodiscard]] torch::Tensor ComputeMaterialStiffness(
        const torch::Tensor & PK2,
        const torch::Tensor & right_cauchy_green);

    [[nodiscard]] torch::Tensor PK2ToCauchy(const torch::Tensor & PK2,
                                            const torch::Tensor & F);

    struct DerivativeResult
    {
      at::Tensor cauchy_stress;
      at::Tensor material_stiffness;
    };

    [[nodiscard]] DerivativeResult ComputeCauchyStressStiffness(
        const at::Tensor & energy,
        const at::Tensor & deformation_gradient,
        const at::Tensor & right_cauchy_green);

    void ApfMatrixToTorchTensor(const apf::Matrix3x3 & apf_matrix,
                                torch::Tensor & tensor);

    [[nodiscard]] apf::DynamicMatrix StressToApfDynamicMatrix(
        const torch::Tensor & cauchy_stress);

    [[nodiscard]] apf::DynamicMatrix TorchTensorToApfDynamicMatrix(
        const torch::Tensor & tensor_in);

    template <typename ViewT>
    void CheckTensorViewCompatibility(const torch::Tensor & tensor,
                                      const ViewT & view)
    {
      KOKKOS_ASSERT(ViewT::rank == tensor.dim());
      for (int i = 0; i < ViewT::rank; ++i)
      {
        KOKKOS_ASSERT(view.extent(i) == tensor.sizes()[i]);
      }
      static_assert(
          Kokkos::SpaceAccessibility<typename ViewT::memory_space,
                                     Kokkos::HostSpace>::accessible,
          "torch array conversion not currently implemented on GPU\n");
    }

    template <typename ViewT>
    void KokkosViewToTorchArray(const ViewT & view, torch::Tensor tensor)
    {
      CheckTensorViewCompatibility(tensor, view);
      static_assert(ViewT::rank == 2 || ViewT::rank == 3,
                    "view must be rank 2 or 3");
      if constexpr (ViewT::rank == 2)
      {
        auto policy = Kokkos::MDRangePolicy<Kokkos::Rank<ViewT::rank>>(
            {0, 0}, {view.extent(1), view.extent(2)});
        Kokkos::parallel_for(
            "KokkosViewToTorchArray", policy,
            KOKKOS_LAMBDA(int i, int j) {
              tensor.index({i,j}) = view(i, j);
            });
      }
      if constexpr (ViewT::rank == 3)
      {
        auto policy = Kokkos::MDRangePolicy<Kokkos::Rank<ViewT::rank>>(
            {0, 0, 0}, {view.extent(0), view.extent(1), view.extent(2)});
        Kokkos::parallel_for(
            "KokkosViewToTorchArray", policy,
            KOKKOS_LAMBDA(int i, int j, int k) {
              tensor.index({i,j,k}) = view(i, j, k);
            });
      }

    }

    template <typename ViewT>
    void TorchArrayToKokkosView(const torch::Tensor & tensor, ViewT & view)
    {
      CheckTensorViewCompatibility(tensor, view);
      static_assert(ViewT::rank == 2 || ViewT::rank == 3,
                    "view must be rank 2 or 3");
      // this should by the torch array type which is float by default
      auto tensor_accessor = tensor.accessor<float, ViewT::rank>();
      if constexpr (ViewT::rank == 2)
      {
        auto policy = Kokkos::MDRangePolicy<Kokkos::Rank<ViewT::rank>>(
            {0, 0}, {view.extent(0), view.extent(1)});
        Kokkos::parallel_for(
            "KokkosViewToTorchArray", policy,
            KOKKOS_LAMBDA(const int i, const int j) {
              view(i, j) = tensor_accessor[i][j];
            });
      }
      else if constexpr (ViewT::rank == 3)
      {
        auto policy = Kokkos::MDRangePolicy<Kokkos::Rank<ViewT::rank>>(
            {0, 0, 0}, {view.extent(0), view.extent(1), view.extent(2)});
        Kokkos::parallel_for(
            "KokkosViewToTorchArray", policy,
            KOKKOS_LAMBDA(const int i, const int j, const int k) {
              view(i, j, k) = tensor_accessor[i][j][k];
            });
      }
    }
  }  // namespace ml
}  // namespace mumfim

#endif  // MUMFIM_SRC_MUMFIM_PYTORCH_TORCHUTILITIES_H
