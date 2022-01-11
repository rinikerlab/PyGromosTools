#from pygromos.files.gromos_system.gromos_system import Gromos_System
#from pygromos.simulations.modules.preset_simulation_modules import emin, md
from pygromos.files.qmmm.qmmm import QMMM

def simulation():
    # binaries
    gromosPP = "/home/fpultar/bin/gromos++/bin"
    gromosXX = "/home/fpultar/bin/gromosXX/bin"

    # folders and title
    system_name = "menthol-dmf"
    work_folder = f"/home/fpultar/Documents/calc/pygromos/{system_name}"

    # files
    in_cnf_path = f"{work_folder}/{system_name}-all-atom_54a7.cnf"
    in_top_path = f"{work_folder}/{system_name}-all-atom_54a7.top"

    # system
    system = Gromos_System(
        work_folder, 
        system_name, 
        in_top_path=in_top_path, 
        in_cnf_path=in_cnf_path,
        in_gromosPP_bin_dir=gromosPP,
        in_gromosXX_bin_dir=gromosXX
        )


    new_system = emin(system)
    new_system.imd = None
    new_system.top = in_top_path
    new_system.cnf = "/home/fpultar/Documents/calc/pygromos/menthol-dmf/emin/analysis/data/emin.cnf"
    
    md(new_system, equilibration_runs=2, simulation_runs=5)

def qmmm():
    QMMM("/home/fpultar/Documents/calc/mdfptools-test/qmmm/test.qmmm")

def main():
    # simulation()
    qmmm()

if __name__ == "__main__":
    main()