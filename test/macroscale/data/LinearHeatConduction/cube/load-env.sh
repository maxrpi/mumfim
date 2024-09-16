if ( test -f /etc/redhat-release && grep -q "7\." "/etc/redhat-release" ) ; then
  echo "RHEL7 not supported"
  return 1
elif ( test -f /etc/redhat-release && grep -q "9\." "/etc/redhat-release" ) ; then
  export MODULEPATH=/opt/scorec/spack/rhel9/v0201_4/lmod/linux-rhel9-x86_64/Core/:/opt/scorec/modules 

  module load gcc/12.3.0-iil3lno mpich/4.1.1-xpoyz4t
  module load simmetrix-simmodsuite/2024.0-240119dev-7abimo4
  module load simmetrix/simModeler/2024.0-240119-dev
  module load parasolid
  module load pumi/develop-simmodsuite-2024.0-240119dev-int32-shared-re4vh42
  module load gdb
  export LD_LIBRARY_PATH=/net/common/meshSim/rhel7/SimModeler2023.0-230121-dev/bin/:${LD_LIBRARY_PATH}
else 
  echo "Could not determine OS"
  return 1
fi

