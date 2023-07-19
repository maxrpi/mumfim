#include <apfDynamicMatrix.h>
#include <apfMatrixUtil.h>
#include <fmt/format.h>
#include <mumfim/exceptions.h>
#include <mumfim/macroscale/materials/TorchMaterial.h>
#include <torch/script.h>
#include <torch/torch.h>
#include <mumfim/torch/TorchUtilities.h>

#include <iostream>

namespace mumfim
{
  TorchMaterial::TorchMaterial(const torch::jit::Module & Model) : model(Model)
  {
    try
    {
      deformation_gradient_tensor = torch::empty({1, 3, 3});
      inputs.emplace_back(torch::empty({1, 3, 3}));
    }
    // forward torch exception as
    catch (const torch::Error & e)
    {
      throw material_error(fmt::format("mumfim TorchMaterial: {}", e.what()));
    }
  }

  TorchMaterial::TorchMaterial(const std::string & model_path)
  {
    try
    {
      model = torch::jit::load(model_path);
      deformation_gradient_tensor = torch::empty({1, 3, 3});
      inputs.emplace_back(torch::empty({1, 3, 3}));
    }
    // forward torch exception as
    catch (const torch::Error & e)
    {
      throw material_error(fmt::format("mumfim TorchMaterial: {}", e.what()));
    }
  }

  MaterialResult TorchMaterial::operator()(
      const apf::Matrix3x3 & deformation_gradient,
      apf::MeshEntity * /* element */,
      int /* integration_point */)
  {
    try
    {
      ml::ApfMatrixToTorchTensor(deformation_gradient, deformation_gradient_tensor);
      // we set requires_grad to false because we do in place tensor operations
      inputs[0].toTensor().set_requires_grad(false);
      ml::ComputeRightCauchyGreen(deformation_gradient_tensor,
                              inputs[0].toTensor());
      // set requires_grad to true so we can take derivatives w.r.t. right
      // cauchy green deformation tensor
      inputs[0].toTensor().set_requires_grad(true);
      auto energy = model.forward(inputs).toTensor();
      auto [cauchy_stress, stiffness] = ml::ComputeCauchyStressStiffness(
          energy, deformation_gradient_tensor, inputs[0].toTensor());

      MaterialResult result;
      // stress is 6x1 for symmetric components, but material result
      // expects the full 3x3 matrix
      result.cauchy_stress = ml::StressToApfDynamicMatrix(cauchy_stress);
      result.material_stiffness = ml::TorchTensorToApfDynamicMatrix(stiffness);
      return result;
    }
    // forward torch exception as
    catch (const torch::Error & e)
    {
      throw material_error(fmt::format("mumfim TorchMaterial: {}", e.what()));
    }
  }

}  // namespace mumfim
