import socket

from pygromos.hpc_queuing.submission_systems.dummy import DUMMY
from pygromos.hpc_queuing.submission_systems.local import LOCAL
from pygromos.hpc_queuing.submission_systems.lsf import LSF


def get_submission_system(testing: bool = False):
    if (testing):
        return DUMMY
    if ("eu" in socket.gethostname()):
        return LSF
    else:
        return LOCAL