#ifndef MUMFIM_SRC_MUMFIM_PYTORCH_TORCHUTILITIES_H
#define MUMFIM_SRC_MUMFIM_PYTORCH_TORCHUTILITIES_H
#include <torch/script.h>
#include <apfDynamicMatrix.h>
#include <apfMatrix.h>

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

  }  // namespace torch
}  // namespace mumfim

#endif  // MUMFIM_SRC_MUMFIM_PYTORCH_TORCHUTILITIES_H
