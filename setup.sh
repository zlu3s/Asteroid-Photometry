#!/usr/bin/sh

export PIPELINE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $PIPELINE/venv/bin/activate
alias id_light="${PIPELINE}/bin/id_light.sh $!"
alias reduce="${PIPELINE}/bin/reduce_images.sh $!"
