from typing import List, Dict


class _fileManagment_base_class:
    def __init__(self):
        pass

    def __str__(self, offset: str = "") -> str:
        msg = offset + self.__class__.__name__ + "\n"
        for key, value in vars(self).items():
            if not key.startswith("_") and not issubclass(value.__class__, _fileManagment_base_class):
                if isinstance(value, List):
                    msg += offset + "\t" + key + ":  " + "\n\t".join(map(str, value)) + "\n"
                else:
                    msg += offset + "\t" + key + ":  " + str(value) + "\n"
            elif issubclass(value.__class__, _fileManagment_base_class):
                msg += offset + "\t" + key + ": \n" + value.__str__(offset=offset + "\t\t") + "\n"
        return msg

    def get_script_generation_command(self, var_name: str = None, var_prefixes="system") -> str:

        if isinstance(var_name, type(None)):
            var_name = var_prefixes + self.__class__.__name__

        gen_cmd = "#Generate " + self.__class__.__name__ + "\n"
        gen_cmd += (
            "from "
            + self.__module__
            + " import "
            + self.__class__.__name__
            + " as "
            + self.__class__.__name__
            + "_obj"
            + "\n"
        )
        constr_str = var_name + " = " + self.__class__.__name__ + "_obj( "
        for key, value in vars(self).items():
            if not key.startswith("_") and not issubclass(value.__class__, _fileManagment_base_class):

                if isinstance(value, List):
                    gen_cmd += var_prefixes + "_" + key + ' = [ "' + '",\n\t"'.join(map(str, value)) + '" ]\n'
                else:
                    gen_cmd += var_prefixes + "_" + key + ' = "' + str(value) + '"\n'
                constr_str += key + "=" + var_prefixes + "_" + key + ", "

            elif issubclass(value.__class__, _fileManagment_base_class):
                gen_cmd += (
                    "\n"
                    + value.get_script_generation_command(
                        var_name=var_prefixes + "_" + value.__class__.__name__, var_prefixes=var_prefixes
                    )
                    + "\n"
                )
                constr_str += key + "=" + var_prefixes + "_" + value.__class__.__name__ + ", "

        constr_str += ")\n"
        cmd = gen_cmd + "\n\n" + constr_str + "\n"

        return cmd
