# SPDX-License-Identifier: MIT

#!/bin/env python3

import sys
import READutils as READut
from diag_Landau import plot_Landau_damping

if len(sys.argv) != 1:
    print('Usage: {}'.format(sys.argv[0]), file=sys.stderr)
    sys.exit(1)

if not READut.valid_directory('.'):
    print('==> There are no VOICEXX results in this directory')
    sys.exit(1)

vxx_res = READut.Read_VOICEXX_results('.')

plot_Landau_damping(vxx_res.electrostatic_potential,
                    vxx_res.time_saved,
                    dirname='.',
                    plotfig=True,
                    validate=True)
