#!/bin/bash
# properties = {properties}

source /g/data1/xe2/.profile


export TMPDIR=$PBS_JOBFS

. raijin/modules.sh

set -ueo pipefail
{exec_job}
