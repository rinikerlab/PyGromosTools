"""
    This abstract parent class, should provide functionality for  checking if binaries are present.
"""
import inspect
import functools
from pygromos.utils.typing import Union, Dict, Callable
from pygromos.utils import bash


class _compiled_program:
    """
    This class contains functions, that are required to detect if a binary program is accesible.
    It tracks the presence of the binary dir and wraps the gromos functionality with wrappers, that check on function call if the binary is present.

    """

    _bin: str
    _force_bin_present: bool
    _check_binary_paths: bool

    _found_binary_dir: Dict[str, bool]  # found? binary dir
    _found_binary: Dict[str, bool]  # found? binary
    _found_binary_paths: Dict[str, str]  # the found binary paths.

    def __init__(
        self, in_bin_dir: Union[str, None], _force_bin_present: bool = True, _check_binary_paths: bool = True
    ) -> Union[str, None]:
        """
            The  _compiled_program parent class can be used, to ensure on runtime, that certain binaries are present.

        Parameters
        ----------
        in_bin_dir : Union[str, None]
            directory that should contain the binaries. If None, the assumption is made, that the path is part of the PATH variable.
        _force_bin_present : bool, optional
            if True, the check_binary or check_binary_folder will throw errors, if they don't find the targets, by default True
        _dont_check_binary : bool, optional
            This is a kill switch for all checks of this class., by default False

        Returns
        -------
        Union[str, None]
            _description_
        """
        # init structures
        self._force_bin_present = _force_bin_present
        self._found_binary_dir = {}
        self._found_binary = {}
        self._found_binary_paths = {}
        self._function_binary = {}
        self._check_binary_paths = _check_binary_paths

        # Check initial status of binaries
        if in_bin_dir is None or in_bin_dir == "None" or in_bin_dir == "":
            self._bin = ""
        elif in_bin_dir.endswith("/"):
            self._bin = in_bin_dir
        else:
            self._bin = in_bin_dir + "/"

        self._check_binary_dir(in_bin_dir=self._bin)
        self._check_all_binaries()
        self.__wrap_programms_with_binary_checks()

    def __getstate__(self):
        """
        preperation for pickling:
        remove the non trivial pickling parts
        """
        return {
            "_bin": self._bin,
            "_check_binary_paths": self._check_binary_paths,
            "_force_bin_present": self._force_bin_present,
        }

    def __setstate__(self, state):
        self.__dict__ = state
        self._found_binary_dir = {}
        self._found_binary = {}
        self._found_binary_paths = {}
        self._function_binary = {}

        # Check initial status of binaries
        self._check_all_binaries()
        self.__wrap_programms_with_binary_checks()

    """
        properties
    """

    @property
    def bin(self) -> str:
        """
        This attribute is the binary directory.
        """
        if not hasattr(self, "_bin") or self._bin == "":
            return None
        else:
            return self._bin

    @bin.setter
    def bin(self, in_bin_dir: str):
        if in_bin_dir is None or in_bin_dir == "None" or in_bin_dir == "":
            self._bin = ""
        elif in_bin_dir.endswith("/"):
            self._bin = in_bin_dir
        else:
            self._bin = in_bin_dir + "/"

        self._check_binary_dir(in_bin_dir=in_bin_dir)
        self._check_all_binaries()
        self.__wrap_programms_with_binary_checks()

    """
        checking functions:
    """

    def _check_binary(self, test_program: str) -> bool:
        """
            this function checks if the binary exists in the given binary dir. (if binary_dir == None -> it will check the PATH variable.)

        Parameters
        ----------
        test_program : str
            name of the binary.

        Returns
        -------
        bool
            if the binary was found

        Raises
        ------
        IOError
            If the binary dir was not found and _dont_check_bin was False (default: False)
        """

        if (test_program in self._found_binary and self._found_binary[test_program]) or not self._check_binary_paths:
            return True

        elif self.bin is not None and bash.command_exists(self._bin + test_program):
            self._found_binary[test_program] = True
            self._found_binary_paths[test_program] = self._bin + test_program
            return self._bin + test_program

        elif self.bin is None and bash.command_exists(test_program):
            self._found_binary[test_program] = True
            self._found_binary_paths[test_program] = test_program
            return test_program

        else:
            self._found_binary[test_program] = False
            if self._force_bin_present:
                raise IOError(
                    "No binary could be found! Please make sure, the program was compiled and the path were passed to this obj."
                    + " provided binary path: "
                    + str(self.bin)
                    + "\tbinary name: "
                    + str(test_program)
                )

    def _check_binary_dir(self, in_bin_dir: str) -> bool:
        """
            this fuction checks and records if the provided binary dir was found.

        Parameters
        ----------
        in_bin_dir : str
            the expected binary dir.

        Returns
        -------
        bool
            Does the binary dir exist?

        Raises
        ------
        IOError
            If the binary dir was not found and _dont_check_bin was False (default: False)
        """

        if (
            in_bin_dir in self._found_binary_dir and self._found_binary_dir[in_bin_dir]
        ) or not self._check_binary_paths:  # did we already check this
            return True
        elif isinstance(in_bin_dir, str) and in_bin_dir != "" and bash.directory_exists(in_bin_dir):
            self._found_binary_dir[in_bin_dir] = True
            if not in_bin_dir.endswith("/"):
                in_bin_dir += "/"
            return True

        elif not self._force_bin_present and (in_bin_dir is None or in_bin_dir == "" or in_bin_dir == "None"):
            self._found_binary_dir[in_bin_dir] = True
            return True
        else:
            self._found_binary_dir[in_bin_dir] = False
            if self._force_bin_present:
                raise IOError(
                    "No binary directory could be found! Please make sure the directory exists! "
                    + " and either pass the path to the binary directory or set the PATH variable. The given folder path was: "
                    + str(in_bin_dir)
                )
            return False

    def _check_all_binaries(self, _force_bin_present: bool = False) -> bool:
        """
        This function checks all present programs of this class, if the binary can be found and executed.
        It does not trigger an Exception, if a binary cannot be found, except force_present is True.
        However it records, if a binary not found in the _function_binary dict.
        The recorded information is used to add restrictive wrappers to the progams, that can throw errors on executions (if _dont_check_bin attribute is False)

        Parameters
        ----------
        force_present : bool, optional
            the function can be restrictive and throw errors, if a function was not found, by default False

        Returns
        -------
        bool
            all binaries present?
        """

        funcs = {key: getattr(self, key) for key in dir(self) if (not key.startswith("_") and key != "bin")}
        tmp_dont_check_bin = self._force_bin_present
        for key, f in funcs.items():
            binary = inspect.signature(f).parameters["_binary_name"].default
            self._check_binary(test_program=binary)
            self._function_binary[binary] = key
        self._force_bin_present = tmp_dont_check_bin

        return all(self._found_binary.values())

    """
        Utils for the binary wrapping to check binary on the fly.
    """

    def _check_binaries_decorator(self, func: Callable) -> Callable:
        """
            This function wrapper adds a binary check before, the function is executed.

        Parameters
        ----------
        func : callable
            function using a binary and having the parameter _binary_name:str

        Returns
        -------
        callable
            wrapped function

        Raises
        ------
        Exception
            if no binary can be found in the function call or definition
        """

        @functools.wraps(func)
        def control_binary(*args, **kwargs) -> any:
            # print("binaryChecker", func.__name__, args, kwargs)

            func_signature = inspect.signature(func)
            if "_binary_name" in kwargs:
                self._check_binary(self._bin + kwargs["_binary_name"])
            elif "_binary_name" in func_signature.parameters:
                self._check_binary(self._bin + func_signature.parameters["_binary_name"].default)
            else:
                raise Exception(
                    "Could not find Binary name in function signature: " + str(func) + "\n found: " + str(kwargs)
                )

            args = list(filter(lambda x: x != self, args))  # avoid double selfing
            return func(self, *args, **kwargs)

        return control_binary

    def __wrap_programms_with_binary_checks(self, remove=False):
        """
            Wraps all functions with the __check_binaries_decorator

        Parameters
        ----------
        remove : bool, optional
            remove all wrappers?, by default False
        """
        if self._check_binary_paths:
            pass

        else:
            v = {}
            for binary, func in self._function_binary.items():
                if remove:  # or self._found_binary[binary]:
                    v[func] = getattr(self.__class__, func)
                else:
                    v[func] = self._check_binaries_decorator(getattr(self, func))
            self.__dict__.update(v)
