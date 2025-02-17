#ifndef MUMFIM_SRC_MUMFIM_MICROSCALE_BATCHEDNEOHOOKEANANALYSIS_H
#define MUMFIM_SRC_MUMFIM_MICROSCALE_BATCHEDNEOHOOKEANANALYSIS_H
#include <Kokkos_Core.hpp>
#include "mumfim/microscale/BatchedRVEAnalysis.h"
#include "mumfim/microscale/MicroTypeDefinitions.h"
#include "mumfim/microscale/TensorUtilities.h"
#include "mumfim/microscale/ContinuumMechanics.h"
namespace mumfim
{
  template <typename Scalar = mumfim::Scalar,
            typename LocalOrdinal = mumfim::LocalOrdinal,
            typename ExeSpace = Kokkos::DefaultExecutionSpace>
  class BatchedNeohookeanAnalysis
      : public BatchedRVEAnalysis<Scalar, LocalOrdinal, ExeSpace>
  {
    public:
    BatchedNeohookeanAnalysis(LocalOrdinal num_rves,
                              Scalar youngs_modulus,
                              Scalar poissons_ratio)
        : BatchedRVEAnalysis<Scalar, LocalOrdinal, ExeSpace>(num_rves)
        , F("f", num_rves)
        , F_updated("f_updated", num_rves)
        , shear_modulus_(youngs_modulus / (2.0 * (1.0 + poissons_ratio)))
        , lambda_((2.0 * shear_modulus_ * poissons_ratio) /
                  (1.0 - 2.0 * poissons_ratio))
    {
      // initialize original deformation gradient
      FillIdentity(F);
    }
    bool run(Kokkos::DualView<Scalar * [3][3], ExeSpace> dfmGrds,
             Kokkos::DualView<Scalar * [6], ExeSpace> sigma,
             bool update_coords = true) final
    {
      KOKKOS_EXPECTS((F.extent(0) == dfmGrds.extent(0)) ||
                     (dfmGrds.extent(0) == 1));
      KOKKOS_EXPECTS(sigma.extent(0) == F.extent(0));
      dfmGrds.sync_device();
      dfmGrds.sync_host();
      Kokkos::View<Scalar * [3][3]> green_strain("green_strain", F.extent(0));
      dfmGrds.sync_device();
      auto dfmGrds_d = dfmGrds.d_view;
      static_assert(std::is_same_v<decltype(dfmGrds_d), decltype(F)>);
      UpdateDeformationGradient(dfmGrds_d, F, F_updated);
      auto left_cauchy_green = ComputeLeftCauchyGreenDeformation(F_updated);

      // to enable copy to device (w/o this ptr)
      auto shear_modulus = shear_modulus_;
      auto lambda = lambda_;
      Kokkos::parallel_for(
          F.extent(0), KOKKOS_LAMBDA(const int i) {
            auto sigma_i = Kokkos::subview(sigma.d_view, i, Kokkos::ALL());
            auto F_i =
                Kokkos::subview(F_updated, i, Kokkos::ALL(), Kokkos::ALL());
            auto J = Determinant<3>(F_i);
            // printf("FUpdated: %f %f %f %f %f %f %f %f %f \n", F_i(0, 0),
            //        F_i(0, 1), F_i(0, 2), F_i(1, 0), F_i(1, 1), F_i(1, 2),
            //        F_i(2, 0), F_i(2, 1), F_i(2, 2));
            sigma_i(0) = (left_cauchy_green(i, 0, 0) - 1) * shear_modulus / J +
                         lambda / J * log(J);
            sigma_i(1) = (left_cauchy_green(i, 1, 1) - 1) * shear_modulus / J +
                         lambda / J * log(J);
            sigma_i(2) = (left_cauchy_green(i, 2, 2) - 1) * shear_modulus / J +
                         lambda / J * log(J);
            sigma_i(3) = (left_cauchy_green(i, 1, 2)) * shear_modulus / J;
            sigma_i(4) = (left_cauchy_green(i, 0, 2)) * shear_modulus / J;
            sigma_i(5) = (left_cauchy_green(i, 0, 1)) * shear_modulus / J;
            // printf("sigma %f %f %f %f %f %f : % f %f %f\n", sigma_i(0),
            //        sigma_i(1), sigma_i(2), sigma_i(3), sigma_i(4),
            //        sigma_i(5), lambda, shear_modulus, J);
          });
      sigma.modify_device();
      if (update_coords)
      {
        Kokkos::deep_copy(F, F_updated);
        // This is needed for the finite difference calculation of the stress
        //Kokkos::deep_copy(this->current_stress_.template view<ExeSpace>(),
        //                  sigma.template view<ExeSpace>());
        //this->current_stress_.template modify<ExeSpace>();
      }
      return true;
    }
    void computeMaterialStiffness(
        Kokkos::DualView<Scalar * [6][6], ExeSpace> C) final
    {
      C.sync_device();
      auto lambda = lambda_;
      auto shear_modulus = shear_modulus_;
      auto F_local = F;
      Kokkos::parallel_for(
          C.extent(0), KOKKOS_LAMBDA(int i) {
            auto C_i =
                Kokkos::subview(C.d_view, i, Kokkos::ALL(), Kokkos::ALL());
            auto F_i = Kokkos::subview(F_local, i, Kokkos::ALL(), Kokkos::ALL());
            auto J = Determinant<3>(F_i);
            auto lambda_prime = lambda / J;
            auto shear_modulus_prime = (shear_modulus - lambda * log(J)) / J;
            C_i(0, 0) = lambda_prime + 2 * shear_modulus_prime;
            C_i(1, 1) = lambda_prime + 2 * shear_modulus_prime;
            C_i(2, 2) = lambda_prime + 2 * shear_modulus_prime;
            C_i(0, 1) = lambda_prime;
            C_i(0, 2) = lambda_prime;
            C_i(1, 0) = lambda_prime;
            C_i(1, 2) = lambda_prime;
            C_i(2, 0) = lambda_prime;
            C_i(2, 1) = lambda_prime;
            C_i(3, 3) = shear_modulus_prime;
            C_i(4, 4) = shear_modulus_prime;
            C_i(5, 5) = shear_modulus_prime;
          });
    }

    private:
    // explicitly specify voids in F for template matching. Eventually functions
    // should be fixed to not require this...
    Kokkos::View<Scalar * [3][3], ExeSpace, void, void> F;
    Kokkos::View<Scalar * [3][3], ExeSpace, void, void> F_updated;
    Scalar shear_modulus_;
    Scalar lambda_;
  };
}  // namespace mumfim
#endif  // MUMFIM_SRC_MUMFIM_MICROSCALE_BATCHEDNEOHOOKEANANALYSIS_H
