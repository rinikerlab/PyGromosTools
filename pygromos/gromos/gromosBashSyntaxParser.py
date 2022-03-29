from pygromos.utils.typing import Union, List


class gromosBashSyntaxParser:
    """
    Helper class to parse general gromos bash syntax

    all methods should be static
    """

    @staticmethod
    def multiplyArgumentParser(args: Union[str, List[str]], multiplier: Union[int, List[int]] = 1) -> str:
        """
        Parser for multiplier syntax to gromos scripts

        example: com_top @topo 1:Protein 2:Na 2:Cl 100:SPC .....

        Parameters
        ----------
        args : str or list(str)
            The actual argument which should be multiplied (ex. top)
        multiplier : int or list(int)
            the multiplier for each argument provided in args
        """
        command = ""
        if multiplier != 1:
            if type(args) == list and len(args) >= 1:
                if len(args) != len(multiplier):
                    raise ValueError("multiplier does not match the number of arguments provided!")
                else:
                    for mult, topo in zip(multiplier, args):
                        command += str(mult) + ":" + topo + " "
            else:
                command = str(multiplier) + ":" + args
        else:
            if type(args) == list and len(args) >= 1:
                command = " ".join(args)
            else:
                command = args
        return command

    @staticmethod
    def atomSliceParser():
        raise NotImplementedError("WIP")

    @staticmethod
    def moleculeSliceParser():
        raise NotImplementedError("WIP")
