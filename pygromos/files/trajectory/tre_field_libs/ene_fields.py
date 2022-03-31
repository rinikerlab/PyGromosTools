from pygromos.utils.typing import List


class gromos_tre_block_names_table:
    totals_subblock_names: List[str]
    eds_subblock_names_singleState: List[str]
    lam_subblock_names_singleLam: List[str]


class gromos_2015_tre_block_names_table(gromos_tre_block_names_table):
    totals_subblock_names = [
        "totene",
        "totkin",
        "totpot",
        "totcov",
        "totbond",
        "totangle",
        "totimproper",
        "totdihedral",
        "totcrossdihedral",
        "totnonbonded",
        "totlj",
        "totcrf",
        "totls",
        "totlspair",
        "totlsreal",
        "totlsk",
        "totlsa",
        "totlsself",
        "totlssurf",
        "totpolself",
        "totspecial",
        "totsasa",
        "totsasavol",
        "totconstraint",
        "totdisres",
        "totdisfieldres",
        "totdihres",
        "totposres",
        "totjval",
        "totxray",
        "totle",
        "totorder",
        "totsymm",
        "eds_vr,entropy",
        "totqm",
        "totbsleus",
        "totrdc",
        "wip1",
    ]

    eds_subblock_names_singleState = ["total", "nonbonded", "special"]
    eds_subblock_names = (
        None  # is generated on the fly in get_eds of TRE - depends on num_states -> simulation specific
    )

    lam_subblock_names_singleLam = [
        "A_e_lj",
        "B_e_lj",
        "A_e_crf",
        "B_e_crf",
        "AB_kinetic",
        "AB_bond",
        "AB_angle",
        "AB_improper",
        "AB_disres",
        "AB_dihres",
        "AB_disfld",
    ]
    lam_subblock_names = (
        None  # is generated on the fly in get_eds of TRE - depends on num_states -> simulation specific
    )


class gromos_2020_tre_block_names_table(gromos_tre_block_names_table):
    totals_subblock_names = [
        "totene",
        "totkin",
        "totpot",
        "totcov",
        "totbond",
        "totangle",
        "totimproper",
        "totdihedral",
        "totcrossdihedral",
        "totnonbonded",
        "totlj",
        "totcrf",
        "totls",
        "totlspair",
        "totlsreal",
        "totlsk",
        "totlsa",
        "totlsself",
        "totlssurf",
        "totpolself",
        "totspecial",
        "totsasa",
        "totsasavol",
        "totconstraint",
        "totdisres",
        "totdisfieldres",
        "totdihres",
        "totposres",
        "totjval",
        "totxray",
        "totle",
        "totorder",
        "totsymm",
        "eds_vr,entropy",
        "totqm",
        "totbsleus",
        "totrdc",
        "wip1",
        "wip2",
        "wip3",
        "wip4",
        "wip5",
        "wip6",
        "wip7",
    ]

    eds_subblock_names_singleState = ["total", "nonbonded", "special", "offset"]
    eds_subblock_names = (
        None  # is generated on the fly in get_eds of TRE - depends on num_states -> simulation specific
    )

    lam_subblock_names_singleLam = [
        "A_e_lj",
        "B_e_lj",
        "A_e_crf",
        "B_e_crf",
        "AB_kinetic",
        "AB_bond",
        "AB_angle",
        "AB_improper",
        "AB_disres",
        "AB_dihres",
        "AB_disfld",
    ]
    lam_subblock_names = (
        None  # is generated on the fly in get_eds of TRE - depends on num_states -> simulation specific
    )


class gromos_2021_tre_block_names_table(gromos_tre_block_names_table):
    totals_subblock_names = [
        "totene",
        "totkin",
        "totpot",
        "totcov",
        "totbond",
        "totangle",
        "totimproper",
        "totdihedral",
        "totcrossdihedral",
        "totnonbonded",
        "totlj",
        "totcrf",
        "totls",
        "totlspair",
        "totlsreal",
        "totlsk",
        "totlsa",
        "totlsself",
        "totlssurf",
        "totpolself",
        "totspecial",
        "totsasa",
        "totsasavol",
        "totconstraint",
        "totdisres",
        "totdisfieldres",
        "totdihres",
        "totposres",
        "totjval",
        "totxray",
        "totle",
        "totorder",
        "totsymm",
        "eds_vmix",
        "eds_vr",
        "eds_emax",
        "eds_emin",
        "eds_globmin",
        "eds_globminfluc",
        "entropy",
        "totqm",
        "totbsleus",
        "totrdc",
        "totangres",
        "wip1",
        "wip2",
        "wip3",
        "wip4",
        "wip5",
        "wip6",
    ]

    eds_subblock_names_singleState = ["total", "nonbonded", "special", "offset"]
    eds_subblock_names = (
        None  # is generated on the fly in get_eds of TRE - depends on num_states -> simulation specific
    )

    lam_subblock_names_singleLam = [
        "A_e_lj",
        "B_e_lj",
        "A_e_crf",
        "B_e_crf",
        "AB_kinetic",
        "AB_bond",
        "AB_angle",
        "AB_improper",
        "AB_disres",
        "AB_dihres",
        "AB_disfld",
    ]
    lam_subblock_names = (
        None  # is generated on the fly in get_eds of TRE - depends on num_states -> simulation specific
    )
