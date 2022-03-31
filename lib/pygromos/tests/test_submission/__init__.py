import tempfile
from pygromos.tests.out_testresults import test_dirs

out_test_root_dir = tempfile.mkdtemp(dir=test_dirs, prefix="tmp_hpc_queuing_")
