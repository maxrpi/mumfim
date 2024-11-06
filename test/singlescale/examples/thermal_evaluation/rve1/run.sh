#!/usr/bin/env bash

BASE_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
ENV_SCRIPT=${BASE_DIR}/load-env.sh


source ${ENV_SCRIPT}  &> /dev/null
if [ $? -ne 0 ]; then
  echo 'Could not load appropriate modules. See "${ENV_SCRIPT}"'
  exit 1
fi

EXEDIR=/lore/bloomm2/mumfim/mumfim/build-linux-rhel9-skylake-3tts4ek/spack-build-3tts4ek/src/
EXE=${EXEDIR}/mumfim_thermal_evaluation
MODELFILE=${BASE_DIR}/rve1.dmg
MESHFILE=${BASE_DIR}/rve1.smb
ATTRIBUTEFILE=${BASE_DIR}/rve1.yaml
CASENAME=rve1
AMSI_OPTIONS=${BASE_DIR}/amsi.yaml
GDB="/opt/scorec/spack/rhel9/v0201_4/install/linux-rhel9-x86_64/gcc-7.4.0/gdb-13.1-z4rrapkqramzepkhrfj4tkcmvbimrsik/bin/gdb -args"

DEBUG_RUNCOMMAND=(
  mpirun -np 1  -hostfile hosts.txt \
  ${GDB} \
  ${EXE} \
  -g ${MODELFILE}  -m ${MESHFILE} -b ${ATTRIBUTEFILE} -c ${CASENAME} -a ${AMSI_OPTIONS}
  )

RUNCOMMAND=(
  $EXE -np 1 \
  -g ${MODELFILE}  -m ${MESHFILE} -b ${ATTRIBUTEFILE} -c ${CASENAME} -a ${AMSI_OPTIONS}
)

echo "Running command:
  ${RUNCOMMAND[@]}
  --------------------------"

${RUNCOMMAND[@]}
