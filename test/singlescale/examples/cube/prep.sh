#!/usr/bin/env bash

MODEL_TRAITS_INSTALL='/lore/maxim/mt/bin'

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
ENV_SCRIPT=${SCRIPT_DIR}/load-env.sh

source ${ENV_SCRIPT}  &> /dev/null
if [ $? -ne 0 ]; then
  echo 'Could not load appropriate modules. See "${ENV_SCRIPT}"'
  exit 1
fi

convert  cube.smd cube-case1.sms cube.smb #--native-model=cube_nat.x_t
simTranslate cube_nat.x_t cube.smd simplecube.smd
mdlConvert cube.smd cube.dmg
${MODEL_TRAITS_INSTALL}/smd2yaml simplecube.smd > simplecube.yaml
