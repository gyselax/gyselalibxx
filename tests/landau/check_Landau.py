# SPDX-License-Identifier: MIT

#!/bin/env python3

import sys
import READutils as READut
from diag_VOICE_XVx import plot_Phi_growthrate_frequency

if len(sys.argv) != 1:
    print('Usage: {}'.format(sys.argv[0]), file=sys.stderr)
    sys.exit(1)

if not READut.valid_directory('.'):
    print('==> There are no VOICEXX results in this directory')
    sys.exit(1)

vxx_res = READut.Read_VOICEXX_results('.')

plot_Phi_growthrate_frequency (vxx_res.electrostatic_potential,
                               vxx_res.time_saved,
                               dirname='.',
                               plotfig=True,
                               validate=True,
                               growth_rate_theory=-0.153,
                               frequency_theory=1.4156)
