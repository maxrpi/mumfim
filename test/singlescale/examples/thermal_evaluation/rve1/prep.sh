#!/usr/bin/env bash

MODEL_TRAITS_INSTALL='/lore/maxim/mt/bin'

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
ENV_SCRIPT=${SCRIPT_DIR}/load-env.sh

source ${ENV_SCRIPT}  &> /dev/null
if [ $? -ne 0 ]; then
  echo 'Could not load appropriate modules. See "${ENV_SCRIPT}"'
  exit 1
fi

convert  rve1.smd rve1-case1.sms rve1.smb --native-model=rve1_nat.x_t
#simTranslate rve1_nat.x_t rve1.smd simplebar.smd
mdlConvert rve1.smd rve1.dmg
${MODEL_TRAITS_INSTALL}/smd2yaml rve1.smd > rve-tags.yaml
