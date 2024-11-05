#!/usr/bin/env bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
ENV_SCRIPT=${SCRIPT_DIR}/load-env.sh

source ${ENV_SCRIPT}  &> /dev/null
if [ $? -ne 0 ]; then
  echo 'Could not load appropriate modules. See "${ENV_SCRIPT}"'
  exit 1
fi

convert  bar.smd bar.sms bar.smb --native-model=bar.xmt_txt
mdlConvert bar.smd bar.dmg
