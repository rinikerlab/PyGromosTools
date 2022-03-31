import __main__
from pygromos.gromos.gromosPP import GromosPP
from pygromos.gromos.gromosXX import GromosXX
from pygromos.gromos.gromosBashSyntaxParser import gromosBashSyntaxParser

# Possilbe but not super helpful ;)
from pygromos.gromos.compile_gromos import install_gromos
from pygromos.utils.utils import dynamic_parser


if __name__ == "__main__":
    # INPUT JUGGELING
    kargs = dynamic_parser(install_gromos, title="Install Gromos")
    print(kargs)
    # Install Command
    install_gromos(**vars(kargs))
