#!/usr/bin/env bash

BASE_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
ENV_SCRIPT=${BASE_DIR}/load-env.sh

source ${ENV_SCRIPT}  &> /dev/null
if [ $? -ne 0 ]; then
  echo 'Could not load appropriate modules. See "${ENV_SCRIPT}"'
  exit 1
fi

EXEDIR=/lore/maxim/mumfim/mumfim/spack-build-kgwxwib/src
EXE=${EXEDIR}/mumfim_singlescale 
MODELFILE=${BASE_DIR}/simplecube.dmg
MESHFILE=${BASE_DIR}/simplecube.smb
ATTRIBUTEFILE=${BASE_DIR}/simplecube.yaml
CASENAME=simplecube
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
