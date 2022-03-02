"""_summary_
"""

from pygromos.utils import bash


class _compiled_programm:

    _bin: str = ""
    __dont_check_bin: bool = False

    def _check_binary_dir(self, in_bin_dir: str, test_programm: str):
        if (
            isinstance(in_bin_dir, str)
            and bash.directory_exists(in_bin_dir)
            and bash.command_exists(f"{in_bin_dir}/pdb2g96")
        ):
            return in_bin_dir + "/"
        elif self.__dont_check_bin or (in_bin_dir is None and bash.command_exists("pdb2g96")):
            return ""
        else:
            raise IOError(
                "No "
                + __name__
                + " binary directory could be found! Please make sure you compiled "
                + __name__
                + " and either pass the path to the binary directory or set the PATH variable"
            )
