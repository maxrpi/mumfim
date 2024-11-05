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
MODELFILE=${BASE_DIR}/bar.dmg
MESHFILE=${BASE_DIR}/bar.smb
ATTRIBUTEFILE=${BASE_DIR}/bar.yaml
CASENAME=thermal
AMSI_OPTIONS=${BASE_DIR}/amsi.yaml
KAPPA_TAG_FILE=${BASE_DIR}/bar.map

RUNCOMMAND=(
  $EXE -np 1 \
  -g ${MODELFILE}  -m ${MESHFILE} -b ${ATTRIBUTEFILE} -c ${CASENAME} -a ${AMSI_OPTIONS} -k ${KAPPA_TAG_FILE}
)

${RUNCOMMAND[@]}
