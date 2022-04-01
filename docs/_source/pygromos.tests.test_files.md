---
date: '2022-04-01T13:55:18.499Z'
docname: _source/pygromos.tests.test_files
images: {}
path: /source-pygromos-tests-test-files
title: pygromos.tests.test_files package
---

# pygromos.tests.test_files package

## Submodules

## pygromos.tests.test_files.general_file_functions module


### _class_ pygromos.tests.test_files.general_file_functions.general_file_tests(methodName='runTest')
Bases: `unittest.case.TestCase`


#### test_copy()

#### test_deepcopy()

#### test_equalAfterCopy()

#### test_equality()

#### test_pickle()
## pygromos.tests.test_files.test_cnf_funcs module


### _class_ pygromos.tests.test_files.test_cnf_funcs.test_cnf(methodName='runTest')
Bases: `pygromos.tests.test_files.general_file_functions.general_file_tests`


#### class_type()
alias of [`pygromos.files.coord.cnf.Cnf`](#pygromos.files.coord.cnf.Cnf)


#### in_file_path(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/in_testfiles/cnf/in_cnf1.cnf_ )

#### root_out(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/out_testresults/tmp_test_files_3j7dxdxv/cnf_5mmpd33o_ )

#### test_clean_posiResNumsByName()

#### test_delete_res()

#### test_parse()

#### test_visualize()

#### test_write_out()
## pygromos.tests.test_files.test_disres module


### _class_ pygromos.tests.test_files.test_disres.test_disres(methodName='runTest')
Bases: `pygromos.tests.test_files.general_file_functions.general_file_tests`


#### class_type()
alias of [`pygromos.files.topology.disres.Disres`](#pygromos.files.topology.disres.Disres)


#### in_file_path(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/in_testfiles/top/disres5.disres_ )

#### root_out(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/out_testresults/tmp_test_files_3j7dxdxv/top_i7aifk30_ )

#### test_parsing_test_file()
## pygromos.tests.test_files.test_disres_file module


### _class_ pygromos.tests.test_files.test_disres_file.test_disres(methodName='runTest')
Bases: `pygromos.tests.test_files.general_file_functions.general_file_tests`


#### class_type()
alias of [`pygromos.files.topology.disres.Disres`](#pygromos.files.topology.disres.Disres)


#### in_file_path(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/in_testfiles/disres/in_disres.dat_ )

#### root_out(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/out_testresults/tmp_test_files_3j7dxdxv/disres_qqiyri7a_ )

#### test_parsing_test_file()
## pygromos.tests.test_files.test_gromosSystem module


### _class_ pygromos.tests.test_files.test_gromosSystem.test_gromos_system(methodName='runTest')
Bases: `unittest.case.TestCase`


#### class_type()
alias of [`pygromos.files.gromos_system.gromos_system.Gromos_System`](#pygromos.files.gromos_system.gromos_system.Gromos_System)


#### input_cnf_path(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/in_testfiles/small_system/6J29.cnf_ )

#### input_top_path(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/in_testfiles/small_system/6J29.top_ )

#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_construct_empty()

#### test_construct_files()

#### test_load()

#### test_rebase()

#### test_write()

#### verbose(_ = Tru_ )
## pygromos.tests.test_files.test_gromosSystem_forcefields module


### _class_ pygromos.tests.test_files.test_gromosSystem_forcefields.test_gromos_system_2016H66(methodName='runTest')
Bases: `pygromos.tests.test_files.test_gromosSystem_forcefields.test_gromos_system_forcefields`


#### ff(_ = <pygromos.files.forcefield.gromos.gromosff.GromosFF object_ )

#### top_residue_list(_ = ['MTL'_ )

### _class_ pygromos.tests.test_files.test_gromosSystem_forcefields.test_gromos_system_54A7(methodName='runTest')
Bases: `pygromos.tests.test_files.test_gromosSystem_forcefields.test_gromos_system_forcefields`


#### ff(_ = <pygromos.files.forcefield.gromos.gromosff.GromosFF object_ )

#### top_residue_list(_ = ['CH3OH'_ )

### _class_ pygromos.tests.test_files.test_gromosSystem_forcefields.test_gromos_system_forcefields(methodName='runTest')
Bases: `unittest.case.TestCase`


#### ff(_ = <pygromos.files.forcefield.gromos.gromosff.GromosFF object_ )

#### file_class()
alias of [`pygromos.files.gromos_system.gromos_system.Gromos_System`](#pygromos.files.gromos_system.gromos_system.Gromos_System)


#### smiles(_ = 'CO_ )

#### test_construct_empty()

#### test_construct_top_from_ff()

#### top_residue_list(_ = ['MTL'_ )

#### verbose(_ = Tru_ )
## pygromos.tests.test_files.test_imd module


### _class_ pygromos.tests.test_files.test_imd.test_imd(methodName='runTest')
Bases: `pygromos.tests.test_files.general_file_functions.general_file_tests`


#### class_type()
alias of [`pygromos.files.simulation_parameters.imd.Imd`](#pygromos.files.simulation_parameters.imd.Imd)


#### in_file_path(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/in_testfiles/imd/in_REEDS1.imd_ )

#### root_out(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/out_testresults/tmp_test_files_3j7dxdxv/imd_4aqplr3p_ )

#### test_edit_REEDS()

#### test_parsing_test_file()

#### test_to_string()

#### test_write_out()
## pygromos.tests.test_files.test_imd_block module


### _class_ pygromos.tests.test_files.test_imd_block.test_imd_block(methodName='runTest')
Bases: `unittest.case.TestCase`


#### test_some_blocks()
## pygromos.tests.test_files.test_mtb module


### _class_ pygromos.tests.test_files.test_mtb.test_mtb(methodName='runTest')
Bases: `pygromos.tests.test_files.general_file_functions.general_file_tests`


#### class_type()
alias of [`pygromos.files.topology.mtb.Mtb`](#pygromos.files.topology.mtb.Mtb)


#### in_file_path(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/data/ff/Gromos2016H66/2016H66.mtb_ )

#### root_out(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/out_testresults/tmp_test_files_3j7dxdxv/mtb_f75_lk_0_ )

#### test_parsing_test_file()

#### test_write()

### _class_ pygromos.tests.test_files.test_mtb.test_mtb_g54a7(methodName='runTest')
Bases: `pygromos.tests.test_files.test_mtb.test_mtb`


#### in_file_path(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/data/ff/Gromos54A7/54a7.mtb_ )

### _class_ pygromos.tests.test_files.test_mtb.test_mtb_orga(methodName='runTest')
Bases: `pygromos.tests.test_files.test_mtb.test_mtb`


#### in_file_path(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/data/ff/Gromos2016H66/2016H66_orga.mtb_ )

#### test_CHE_read()

#### test_all_mtb_solutes_read()
## pygromos.tests.test_files.test_ptp module


### _class_ pygromos.tests.test_files.test_ptp.test_ptp(methodName='runTest')
Bases: `pygromos.tests.test_files.general_file_functions.general_file_tests`


#### class_type()
alias of [`pygromos.files.topology.ptp.Pertubation_topology`](#pygromos.files.topology.ptp.Pertubation_topology)


#### in_file_path(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/in_testfiles/ptp/eds.ptp_ )

#### root_out(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/out_testresults/tmp_test_files_3j7dxdxv/ptp_ji58jh_b_ )

#### test_IO()

#### test_add_eds_state()

#### test_delete_eds_state()

#### test_eds_gen_all_states_ptp()

#### test_get_eds_states()

#### test_new_eds_ptp_from_scratch()

#### test_remove_eds_atoms()
## pygromos.tests.test_files.test_qmmm module


### _class_ pygromos.tests.test_files.test_qmmm.test_qmmm(methodName='runTest')
Bases: `pygromos.tests.test_files.general_file_functions.general_file_tests`


#### class_type()
alias of [`pygromos.files.qmmm.qmmm.QMMM`](#pygromos.files.qmmm.qmmm.QMMM)


#### in_file_path(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/in_testfiles/qmmm/menthol-methanol-dmf.qmmm_ )

#### root_out(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/out_testresults/tmp_test_files_3j7dxdxv/qmmm_kmeqcg93_ )

#### test_parsing_test_file()

#### test_to_string()

#### test_write_out()

### _class_ pygromos.tests.test_files.test_qmmm.test_qmmm_imd(methodName='runTest')
Bases: `pygromos.tests.test_files.general_file_functions.general_file_tests`


#### class_type()
alias of [`pygromos.files.simulation_parameters.imd.Imd`](#pygromos.files.simulation_parameters.imd.Imd)


#### in_file_path(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/in_testfiles/qmmm/md.imd_ )

#### root_out(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/out_testresults/tmp_test_files_3j7dxdxv/qmmm_kmeqcg93_ )

#### test_parsing_test_file()

#### test_to_string()

#### test_write_out()
## pygromos.tests.test_files.test_repdat module


### _class_ pygromos.tests.test_files.test_repdat.test_repdat(methodName='runTest')
Bases: `unittest.case.TestCase`


#### class_type()
alias of [`pygromos.files.otherfiles.repdat.Repdat`](#pygromos.files.otherfiles.repdat.Repdat)


#### in_file_path(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/in_testfiles/repdat/in_REEDS_repdat2_short.dat_ )

#### root_out(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/out_testresults/tmp_test_files_3j7dxdxv/repdat_4zzp0dy__ )

#### test_cleaning_concat()

#### test_get_transitions()

#### test_parsing_new_file()

#### test_parsing_test_file()

#### test_test_transitions()

#### test_write_out()
## pygromos.tests.test_files.test_resnlib module


### _class_ pygromos.tests.test_files.test_resnlib.test_resnlib(methodName='runTest')
Bases: `pygromos.tests.test_files.general_file_functions.general_file_tests`


#### class_type()
alias of [`pygromos.files.otherfiles.residue_library.residue_library`](#pygromos.files.otherfiles.residue_library.residue_library)


#### in_file_path(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/data/pdb2g96.lib_ )

#### root_out(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/out_testresults/tmp_test_files_3j7dxdxv/resnLib_uoo_oji5_ )

#### test_IO()
## pygromos.tests.test_files.test_top module


### _class_ pygromos.tests.test_files.test_top.test_top(methodName='runTest')
Bases: `pygromos.tests.test_files.general_file_functions.general_file_tests`


#### class_type()
alias of [`pygromos.files.topology.top.Top`](#pygromos.files.topology.top.Top)


#### in_file_path(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/in_testfiles/top/test.top_ )

#### root_out(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/out_testresults/tmp_test_files_3j7dxdxv/top_r_x5so39_ )

#### test_parsing_test_file()

### _class_ pygromos.tests.test_files.test_top.test_top_simple(methodName='runTest')
Bases: `pygromos.tests.test_files.general_file_functions.general_file_tests`


#### class_type()
alias of [`pygromos.files.topology.top.Top`](#pygromos.files.topology.top.Top)


#### in_file_path(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/in_testfiles/top/simpleTest.top_ )

#### root_out(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/out_testresults/tmp_test_files_3j7dxdxv/top_r_x5so39_ )

#### test_add_eq_mul()

#### test_additon()

#### test_multiplication()

#### test_parsing_test_file()
## pygromos.tests.test_files.test_trajectory module


### _class_ pygromos.tests.test_files.test_trajectory.test_trc(methodName='runTest')
Bases: `unittest.case.TestCase`


#### class_name()
alias of [`pygromos.files.trajectory.trc.Trc`](#pygromos.files.trajectory.trc.Trc)


#### help_class(_ = TITLE gch found 22 hydrogen atoms in /home/bschroed/Documents/projects/PyGromosTools_release3/examples/example_files/Tutorial_System/a_build_initial_files/peptide.cnf 0 were within 0.1% of minimum energy bond length 22 were assigned new coordinates based on geometry  	>>> Generated with PyGromosTools (riniker group) <<< END POSITION # 	      1 VAL   H1         1    1.241783665    1.501556791    1.518273147     1 VAL   H2         2    1.221597569    1.373896735    1.618130891     1 VAL   N          3    1.196200000    1.413300000    1.529800000     1 VAL   H3         4    1.097093311    1.426247772    1.526603233     1 VAL   CA         5    1.237300000    1.322400000    1.421900000     1 VAL   CB         6    1.116100000    1.268700000    1.345400000     1 VAL   CG1        7    1.041700000    1.369500000    1.257600000     1 VAL   CG2        8    1.157800000    1.148100000    1.261100000     1 VAL   C          9    1.335700000    1.393300000    1.328600000     1 VAL   O         10    1.447800000    1.344200000    1.316600000     2 TYR   N         11    1.297800000    1.514300000    1.288500000     2 TYR   H         12    1.210575283    1.551469583    1.320284758     2 TYR   CA        13    1.381800000    1.594400000    1.198300000     2 TYR   CB        14    1.320800000    1.731300000    1.167700000     2 TYR   CG        15    1.316700000    1.828700000    1.285500000     2 TYR   CD1       16    1.209300000    1.823700000    1.373600000     2 TYR   HD1       17    1.130290031    1.749796485    1.360307326     2 TYR   CD2       18    1.420400000    1.919600000    1.303000000     2 TYR   HD2       19    1.504990885    1.920076950    1.234260240     2 TYR   CE1       20    1.204100000    1.914300000    1.478900000     2 TYR   HE1       21    1.120330759    1.912712834    1.548622271     2 TYR   CE2       22    1.415700000    2.009900000    1.408600000     2 TYR   HE2       23    1.496221894    2.081795752    1.423707128     2 TYR   CZ        24    1.306500000    2.006900000    1.494600000     2 TYR   OH        25    1.303000000    2.103700000    1.590000000     2 TYR   HH        26    1.385890905    2.159332924    1.584162563     2 TYR   C         27    1.517900000    1.617700000    1.264200000     2 TYR   O         28    1.624900000    1.586100000    1.212400000     3 ARG   N         29    1.514400000    1.657100000    1.391100000     3 ARG   H         30    1.427383863    1.670744779    1.438449888     3 ARG   CA        31    1.642400000    1.680000000    1.459700000     3 ARG   CB        32    1.617000000    1.721100000    1.604800000     3 ARG   CG        33    1.539900000    1.617300000    1.686500000     3 ARG   CD        34    1.608900000    1.606000000    1.822500000     3 ARG   NE        35    1.527400000    1.525500000    1.914700000     3 ARG   HE        36    1.438096810    1.493764282    1.882797265     3 ARG   CZ        37    1.565300000    1.492000000    2.038800000     3 ARG   NH1       38    1.686400000    1.525700000    2.085100000     3 ARG   HH11      39    1.713360963    1.500046328    2.177917001     3 ARG   HH12      40    1.749840898    1.576510923    2.026846223     3 ARG   NH2       41    1.485200000    1.418600000    2.117200000     3 ARG   HH21      42    1.514453755    1.394156276    2.209648484     3 ARG   HH22      43    1.396160354    1.388258647    2.083268543     3 ARG   C         44    1.729900000    1.554500000    1.457100000     3 ARG   O         45    1.845400000    1.559500000    1.415000000     4 LYSH  N         46    1.668200000    1.443600000    1.497000000     4 LYSH  H         47    1.571910359    1.447946461    1.523634813     4 LYSH  CA        48    1.738800000    1.314800000    1.502700000     4 LYSH  CB        49    1.646500000    1.196800000    1.533800000     4 LYSH  CG        50    1.572300000    1.213000000    1.666700000     4 LYSH  CD        51    1.497100000    1.082100000    1.691400000     4 LYSH  CE        52    1.396500000    1.097200000    1.805600000     4 LYSH  NZ        53    1.321000000    0.974500000    1.834700000     4 LYSH  HZ1       54    1.350271838    0.901778367    1.772613579     4 LYSH  HZ2       55    1.223278482    0.991693073    1.822253991     4 LYSH  HZ3       56    1.337909275    0.946523380    1.929206006     4 LYSH  C         57    1.798300000    1.287600000    1.364400000     4 LYSH  O         58    1.917500000    1.258100000    1.359000000     5 GLN   N         59    1.716100000    1.293600000    1.260000000     5 GLN   H         60    1.618525916    1.310613276    1.273778479     5 GLN   CA        61    1.769400000    1.275100000    1.124300000     5 GLN   CB        62    1.656700000    1.300200000    1.023800000     5 GLN   CG        63    1.627000000    1.446800000    0.991700000     5 GLN   CD        64    1.523100000    1.473500000    0.882700000     5 GLN   OE1       65    1.501500000    1.390900000    0.794200000     5 GLN   NE2       66    1.460300000    1.590500000    0.891200000     5 GLN   HE21      67    1.390971831    1.615073255    0.823452231     5 GLN   HE22      68    1.482431587    1.653856690    0.965335839     5 GLN   C         69    1.892900000    1.359900000    1.093300000     5 GLN   O1        70    1.991900000    1.296800000    1.050300000     5 GLN   O2        71    1.885100000    1.482600000    1.116000000 END GENBOX     0     0.000000000    0.000000000    0.000000000     0.000000000    0.000000000    0.000000000     0.000000000    0.000000000    0.000000000     0.000000000    0.000000000    0.000000000 END_ )

#### in_file_path(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/in_testfiles/trc/in_test.trc_ )

#### in_file_path_h5(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/in_testfiles/trc/in_test_trc.h5_ )

#### in_file_w_genbox_cnf_path(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/in_testfiles/trc/in_test_genbox.cnf_ )

#### in_file_w_genbox_path(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/in_testfiles/trc/in_test_genbox.trc_ )

#### outpath(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/out_testresults/tmp_test_files_3j7dxdxv/trajs_f1is2xc_/out_trc1.h5_ )

#### test_constructor_empty()

#### test_constructor_trc_file_noTop_path()

#### test_constructor_trc_file_path()

#### test_constructor_trc_h5_file_path()

#### test_to_conf()

#### test_to_trc_file()

#### test_trc_with_boxes_traj()

#### test_write()

#### trc_outpath(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/out_testresults/tmp_test_files_3j7dxdxv/trajs_f1is2xc_/out_.trc.gz_ )

### _class_ pygromos.tests.test_files.test_trajectory.test_tre(methodName='runTest')
Bases: `pygromos.tests.test_files.test_trajectory.traj_standard_tests`


#### class_name()
alias of [`pygromos.files.trajectory.tre.Tre`](#pygromos.files.trajectory.tre.Tre)


#### in_file2_path(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/in_testfiles/tre/in_tre2.tre_ )

#### in_file_eds_path(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/in_testfiles/tre/in_eds.tre_ )

#### in_file_h5_path(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/in_testfiles/tre/in_tre1.tre.h5_ )

#### in_file_lam_path(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/in_testfiles/tre/in_lam.tre_ )

#### in_file_path(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/in_testfiles/tre/in_tre1.tre_ )

#### outpath(_: st_ _ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/out_testresults/tmp_test_files_3j7dxdxv/trajs_f1is2xc_/out_tre1.tre.h5_ )

#### test_get_eds()

#### test_get_lam()

#### test_get_totals()

### _class_ pygromos.tests.test_files.test_trajectory.test_trg(methodName='runTest')
Bases: `pygromos.tests.test_files.test_trajectory.traj_standard_tests`


#### class_name()
alias of [`pygromos.files.trajectory.trg.Trg`](#pygromos.files.trajectory.trg.Trg)


#### in_file_h5_path(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/in_testfiles/trg/test.trg.h5_ )

#### in_file_path(_ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/in_testfiles/trg/test.trg_ )

#### outpath(_: st_ _ = '/home/bschroed/Documents/projects/PyGromosTools/pygromos/tests/out_testresults/tmp_test_files_3j7dxdxv/trajs_f1is2xc_/out_tre1.tre.h5_ )

#### test_get_lambdas()

#### test_get_precalclam()

#### test_get_totals()

### _class_ pygromos.tests.test_files.test_trajectory.traj_standard_tests(methodName='runTest')
Bases: `unittest.case.TestCase`


#### class_name()
alias of `pygromos.files.trajectory._general_trajectory._General_Trajectory`


#### in_file_h5_path(_ = Non_ )

#### in_file_path(_ = Non_ )

#### outpath(_: st_ )

#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_add()

#### test_constructor_empty()

#### test_constructor_trg_file_path()

#### test_constructor_trg_h5_file_path()

#### test_write()
## Module contents
