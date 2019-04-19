#!/bin/bash
# properties = {properties}

source /g/data1/xe2/.profile

set -ueo pipefail

export TMPDIR=$PBS_JOBFS

. raijin/modules.sh

{exec_job}
