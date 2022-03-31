import os
import functools
from pygromos.files._basics import _general_gromos_file
from pygromos.utils.typing import Callable
from pygromos.utils.compiledProgram import _compiled_program


class _gromosClass(_compiled_program):
    def __init__(self, in_bin_dir: str, dummy: bool = False, _check_binary_paths: bool = True) -> str:
        """
          This parent class contains wrappers for gromos functionalities.
          E.g. gromosTypeConverter converts a passed gromos obj (Cnf, Top etc.) to a string path, such it can be passed to the comand line tools.

        Parameters
        ----------
        in_bin_dir : Union[str, None]
            directory containing binaries of the gromos program.
        dummy : bool, optional
            For dummy executions, will not throw errors on binary check fails, by default False
        _dont_check_binary : bool, optional
            This flag removes the checks of the binary presence for this obj. This can make sense if system access is slow!, by default False - checks will be made

        Returns
        -------
        Union[str, None]
            _description_
        """
        super().__init__(in_bin_dir, dummy, _check_binary_paths=_check_binary_paths)

    @staticmethod
    def _gromosTypeConverter(func: Callable) -> Callable:
        """
            This decorator can be used to automatically convert gromos files to the str path, where this obj, was written to.

        Parameters
        ----------
        func: Callable
            function to be decorated

        Returns
        -------
        Callable
            decorated function

        """

        @functools.wraps(func)
        def convert_pyGromos_types(self, *args, **kwargs):
            # print("Converter1: ", func.__name__, args, kwargs)

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

            # print("Converter2: ", func.__name__, self, nargs, kwargs)
            return func(self, *nargs, **kwargs)

        return convert_pyGromos_types
