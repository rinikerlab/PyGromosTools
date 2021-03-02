from pygromos.files.gromos_system.gromos_system import gromos_system

import importlib
if(importlib.util.find_spec("openforcefield") != None):
    from pygromos.files.gromos_system.serenityff import *
    from pygromos.files.gromos_system.openforcefield2gromos import *
