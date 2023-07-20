#include <apf.h>
#include <mpi.h>
#include <iostream>
#include "mumfim/microscale/MicroFOParams.h"
#include "mumfim/microscale/FiberNetworkLibrary.h"
#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
//#include <mumfim/microscale/BatchedFiberRVEAnalysisExplicit.h>
#include <memory>
#include <PCU.h>
#include <mumfim/microscale/BatchedTorchAnalysis.h>
template <typename T>
void stressToMat(int idx, T stress_view, apf::Matrix3x3 & stress)
{
  stress[0][0] = stress_view(idx,0);
  stress[1][1] = stress_view(idx,1);
  stress[2][2] = stress_view(idx,2);
  stress[1][2] = stress_view(idx,3);
  stress[2][1] = stress_view(idx,3);
  stress[0][2] = stress_view(idx,4);
  stress[2][0] = stress_view(idx,4);
  stress[0][1] = stress_view(idx,5);
  stress[1][0] = stress_view(idx,5);
}
int main(int argc, char * argv[])
{
  amsi::MPI mpi{argc, argv, MPI_COMM_WORLD};
#ifdef MICRO_USING_PETSC
  las::initPETScLAS(&argc, &argv, MPI_COMM_WORLD);
#endif
  Kokkos::ScopeGuard kokkos(argc, argv);
  if(argc != 3)
  {
    std::cerr<<"Usage: "<<argv[0]<<" job.yaml num_rves"<<std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  PCU_Switch_Comm(MPI_COMM_SELF);
  auto BatchNum = std::atoi(argv[2]);
  std::vector<mumfim::MicroCase> cases;
  mumfim::loadMicroFOFromYamlFile(argv[1], cases);
  mumfim::printMicroFOCase(cases[0]);
  std::string file_name = cases[0].pd.meshFile;

  using ExeSpace = Kokkos::DefaultExecutionSpace;
  // using ExeSpace = Kokkos::Serial;
  // using Scalar = float;//mumfim::Scalar;
  using Scalar = mumfim::Scalar;
  using Ordinal = mumfim::LocalOrdinal;
  mumfim::BatchedTorchAnalysis<Scalar,Ordinal,Kokkos::HostSpace> batched_analysis(file_name, BatchNum);
  Kokkos::DualView<Scalar*[3][3],ExeSpace> deformation_gradient("deformation gradients",BatchNum);
  Kokkos::DualView<Scalar*[6][6],ExeSpace> stiffness("stiffness",BatchNum);
  Kokkos::DualView<Scalar*[6], ExeSpace> stress("stress",BatchNum);
  auto deformation_gradient_h = deformation_gradient.h_view;
  for(int i=0; i<BatchNum; ++i)
  {
    for(int ei=0; ei<3; ++ei)
    {
      for(int ej=0; ej<3; ++ej)
      {
        deformation_gradient_h(i,ei,ej) = cases[0].pd.deformationGradient[ei*3+ej];
      }
    }
  }
  deformation_gradient.modify<Kokkos::HostSpace>();
  Kokkos::Timer timer;
  double time1 = timer.seconds();
  auto success = !batched_analysis.run(deformation_gradient, stress);
  stress.sync<Kokkos::HostSpace>();
  auto stress_h = stress.h_view;
  std::cout<<"Stress :";
  for (int i = 0; i < 6; ++i)
  {
    std::cout<<stress_h(0,i)<<" ";
  }
  std::cout<<"\n";
  batched_analysis.computeMaterialStiffness(stiffness);
  stiffness.sync<Kokkos::HostSpace>();
  auto stiffness_h = stiffness.h_view;
  std::cout<<"Stiffness :\n";
  for (int i = 0; i < 6; ++i)
  {
    for (int j = 0; j < 6; ++j)
    {
      std::cout << stiffness_h(0, i, j) << " ";
    }
    std::cout << "\n";
  }
  std::cout << std::endl;
  double time2 = timer.seconds();
  std::cout << "Took: " << time2 - time1 << " seconds." << std::endl;
  return success;
}
