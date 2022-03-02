"""
    This abstract parent class, should provide functionality for  checking if binaries are present.
"""
import inspect
import functools
from typing import Union, Dict
from pygromos.utils import bash


class _compiled_programm:

    _bin: str = ""
    __dont_check_bin: bool = False

    _found_binary_dir: Dict[str, bool]
    _found_binary: Dict[str, bool]

    @property
    def bin(self) -> Union[str, None]:
        if not hasattr(self, "_bin") or self._bin == "":
            return None
        else:
            return self._bin

    @bin.setter
    def bin(self, in_bin_dir: str):
        self._bin = self._check_binary_dir(in_bin_dir=in_bin_dir)

    def _check_binary_dir_constructor(self, in_bin_dir: str) -> Union[str, None]:
        self._found_binary_dir = {}
        self.__wrap_programms_with_binary_checks()

        if not in_bin_dir is None:
            return self._check_binary_dir(in_bin_dir=in_bin_dir)
        else:
            return None

    def _check_binary(self, in_bin_dir: str, test_program: str) -> str:
        self._check_binary_dir(in_bin_dir=in_bin_dir)
        if test_program in self._found_binary and self._found_binary[test_program]:
            return test_program if (in_bin_dir is None) else in_bin_dir + "/" + test_program
        elif bash.command_exists(f"{in_bin_dir}/" + test_program):
            self._found_binary[test_program] = True
            return in_bin_dir + "/" + test_program
        elif self.__dont_check_bin or (in_bin_dir is None and bash.command_exists(test_program)):
            self._found_binary[test_program] = True
            return test_program
        else:
            self._found_binary[test_program] = False
            raise IOError(
                "No binary could be found! Please make sure, the program was compiled and the path were passed to this obj."
                + " provided binary path: "
                + str(in_bin_dir)
                + "\tbinary name: "
                + str(test_program)
            )

    def _check_binary_dir(self, in_bin_dir: str) -> str:
        if in_bin_dir in self._found_binary_dir and self._found_binary_dir[in_bin_dir]:  # did we already check this
            return "" if (in_bin_dir is None) else in_bin_dir + "/"
        elif isinstance(in_bin_dir, str) and bash.directory_exists(in_bin_dir):
            self._found_binary_dir[in_bin_dir] = True
            return in_bin_dir + "/"
        elif self.__dont_check_bin or in_bin_dir is None:
            self._found_binary_dir[in_bin_dir] = True
            return ""
        else:
            self._found_binary_dir[in_bin_dir] = False
            raise IOError(
                "No binary directory could be found! Please make sure the directory exists! "
                + " and either pass the path to the binary directory or set the PATH variable. The given folder path was: "
                + str(in_bin_dir)
            )

    def __check_binaries_decorator(self, func: callable) -> callable:
        @functools.wraps(func)
        def control_binary(*args, **kwargs) -> any:
            print("WUHU!")
            func_signature = inspect.signature(func)
            if "_binary_name" in func_signature.parameters:
                bin_name = func_signature.parameters["_binary_name"].default
                self._check_binary(self._bin, bin_name)
            else:
                raise Exception(
                    "Could not find Binary name in function signature: "
                    + str(func)
                    + "\n found: "
                    + str(func_signature)
                )
            return func(self, *args, **kwargs)

        return control_binary

    def __wrap_programms_with_binary_checks(self, remove=False):
        class_funcs = self.__class__
        func = [k for k in dir(class_funcs) if (not k.startswith("_") and k != "bin")]
        self._found_binary = {f: False for f in func}

        v = {}
        if remove:
            v = {f: getattr(class_funcs, f) for f in func}
        else:
            v = {f: self.__check_binaries_decorator(getattr(class_funcs, f)) for f in func}
        self.__dict__.update(v)
