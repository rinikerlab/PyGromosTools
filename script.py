from pygromos.files.gromos_system.gromos_system import Gromos_System
from pygromos.simulations.modules.preset_simulation_modules import emin, md
from pygromos.files.qmmm.qmmm import QMMM
from pygromos.files.simulation_parameters.imd import Imd

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

    print(system)
    print(new_system)
    #new_system.imd = None
    #new_system.top = in_top_path
    #new_system.cnf = "/home/fpultar/Documents/calc/pygromos/menthol-dmf/emin/analysis/data/emin.cnf"
    
    #md(new_system, equilibration_runs=2, simulation_runs=5)

def qmmm():
    qmmm_file = QMMM("/home/fpultar/Documents/calc/mdfptools-test/qmmm/xphos-methanol-dmf.qmmm")

def imd():
    imd_path = "/home/fpultar/Documents/calc/pygromos-qmmm/md.imd"
    imd_file = Imd(imd_path)
    print(imd_file.QMMM.NTQMMM)
    print(imd_file.QMMM.NTQMSW)
    print(imd_file.QMMM.MMSCAL)
    print()

    imd_file.QMMM.NTQMSW = 4

    print(imd_file.QMMM.NTQMMM)
    print(imd_file.QMMM.NTQMSW)
    print(imd_file.QMMM.MMSCAL)
    print()

    imd_file.TITLE.content = "Felix File"
    print(imd_file.TITLE.content)
    print()

    print(imd_file)

def main():
    simulation()
    #qmmm()
    # imd()

if __name__ == "__main__":
    main()