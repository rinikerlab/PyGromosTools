from pygromos.files.gromos_system.gromos_system import Gromos_System

import importlib

if importlib.util.find_spec("openforcefield") is not None:
    from pygromos.files.gromos_system.ff.serenityff import serenityff
    from pygromos.files.gromos_system.ff.openforcefield2gromos import openforcefield2gromos
