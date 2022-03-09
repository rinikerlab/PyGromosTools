from typing import Union
from pygromos.utils.compiledProgram import _compiled_program

import os
import functools
from pygromos.files._basics import _general_gromos_file


class _gromosClass(_compiled_program):
    def __init__(self, in_bin_dir: str, dummy: bool = False) -> Union[str, None]:

        """
        for func in dir(self):
            if callable(getattr(self, func)) and not func.startswith("__") and not func.startswith("_"):
                setattr(self, func, self._gromosTypeConverter(getattr(self, func)))
                pass
        """
        super().__init__(in_bin_dir, dummy)

    def _gromosTypeConverter(func) -> callable:
        """
            This decorator can be used to automatically convert
        Parameters
        ----------
        func

        Returns
        -------
        callable

        """

        @functools.wraps(func)
        def convert_pyGromos_types(self, *args, **kwargs):
            print("Converter1: ", func.__name__, args, kwargs)

            # no key-word parameters
            nargs = []
            for v in args:
                if isinstance(v, _general_gromos_file._general_gromos_file):
                    if v.path is None or not os.path.exists(v.path):
                        raise IOError(
                            "please write out the "
                            + str(v.__name__)
                            + " first to use the function "
                            + str(func.__name__)
                            + ""
                        )
                    v = v.path
                nargs.append(v)

            nargs = list(filter(lambda x: x != self, nargs))  # avoid double selfing

            # key-value - parameters
            for k, v in kwargs.items():
                if isinstance(v, _general_gromos_file._general_gromos_file):
                    if v.path is None or not os.path.exists(v.path):
                        raise IOError(
                            "please write out the " + k + " first to use the function " + str(func.__name__) + ""
                        )
                    kwargs[k] = v.path

            print("Converter2: ", func.__name__, self, nargs, kwargs)
            return func(self, *nargs, **kwargs)

        return convert_pyGromos_types
