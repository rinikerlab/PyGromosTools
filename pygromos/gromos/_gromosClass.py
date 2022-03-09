from typing import Union
from pygromos.utils.compiledProgram import _compiled_program


class _gromosClass(_compiled_program):
    def __init__(self, in_bin_dir: str, dummy: bool = False) -> Union[str, None]:

        pass
        """
        for func in dir(self):
            if callable(getattr(self, func)) and not func.startswith("__"):
                setattr(self, func, self._gromosTypeConverter(getattr(self, func)))
                pass
        """
        super().__init__(in_bin_dir, dummy)
