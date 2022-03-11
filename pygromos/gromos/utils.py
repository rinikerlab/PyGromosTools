import os
import functools
from pygromos.files._basics import _general_gromos_file


def gromosTypeConverter(func) -> callable:
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

        # key-value - parameters
        for k, v in kwargs.items():
            if isinstance(v, _general_gromos_file._general_gromos_file):
                if v.path is None or not os.path.exists(v.path):
                    raise IOError("please write out the " + k + " first to use the function " + str(func.__name__) + "")
                kwargs[k] = v.path

        print("Converter2: ", func.__name__, args, kwargs)
        # *nargs,
        return func(self, **kwargs)

    return convert_pyGromos_types
