#for nice module access

from modules.general_simulation_modules import simulation
from modules.preset_simulation_modules import emin, md, sd
from modules.ti_modules import TI_sampling

from approaches.hvap_calculation.hvap_calculation import Hvap_calculation

from hpc_queuing import submission_systems