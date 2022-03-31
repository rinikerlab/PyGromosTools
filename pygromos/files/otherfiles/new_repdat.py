"""
FUNCTIONLIB:            Repdat File
Description:
    From a Replica-Exchange simulation a repdat file will be created, that gives insight on the replica exchanges of the simulation
Author: Benjamin Schroeder
"""

import pandas as pd

from pygromos.files._basics import parser
from pygromos.files.blocks import replica_exchange_blocks as blocks
from pygromos.utils.typing import Dict, List, Union, Repdat_Type


class Repdat(pd.DataFrame):  #
    """Replica exchange statistic file
    This class is a representation for all transition information during a replica exchange run. it adds some useful functionality.


    """

    _gromos_file_ending = "restat"

    # pandas specific parameters:
    @property
    def _constructor(self):
        return Repdat

    _metadata = ["added_property"]
    added_property = 1  # This will be passed to copies

    # gromos fun
    SYSTEM: blocks.repex_system
    DATA: pd.DataFrame

    # transition_traces[replica][["trials", "positions". "state_pot"]]
    transition_traces: Dict[int, Dict[str, List[float]]] = None

    # count_state_per_position[replicaposition][["tot_nup", "tot_ndown", "states_index", "dt", "dt_nup", "dt_ndown"]]
    count_state_per_position: Dict[int, Dict[str, Union[List[int], int]]] = None

    # count_state_per_position[replica]
    replica_round_trips: Dict[int, int] = None

    def __init__(self, input_path: str):
        """Repdat_Constructor

        Parameters
        ----------
        input_path :    str
             path to gromos repdadat file
        """

        if type(input_path) is str:
            system, df = parser.read_repdat(input_path, Vj_header=True)
            self.system = system
            self.DATA = df  # future data field a pandas frame!
            self.path = input_path
        else:
            raise NotImplementedError("Not correct yet!")

    def _clean_replica_round_trips(self, replica_round_trips: Dict[int, int]) -> Dict[int, int]:
        """_clean_replica_round_trips - privat

            This function cleans up so that the minimal rountrip number in a roundtrip dict is 0

        Parameters
        ----------
        replica_round_trips : Dict[int:int]
            a dictionary containing all replica roundtrip counts

        Returns
        -------
        Dict[int:int]
            a dictionary containing all replica roundtrip counts with lowest value 0

        """

        # clean up indices that are -1
        clean_replica_round_trips = {}
        for key, item in replica_round_trips.items():
            if item != -1:
                clean_replica_round_trips.update({key: item})
            else:
                clean_replica_round_trips.update({key: 0})
        return clean_replica_round_trips

    def _caculate_transition_traces(self) -> None:
        """
        TODO: refactor code!
        ..autofunction: _caculate_transition_traces
            calculates the transition traces for all replicas from raw data and stores them in self.transition_traces.
            In the end you recieve the trace a replica coord system moved through s dist
            format: {replicaID: {[trials...], [position...], [PotE...]}}

        :return: None
        :rtype: None
        """

        replicas = len(self.system.s)

        # follow transitions of one state
        transition_dict = {
            x: x for x in range(1, replicas + 1)
        }  # keeps track of unique id and current replica position.
        tmp_dict = {x: x for x in range(1, replicas + 1)}
        transition_result_dict = {
            x: {"trial": [], "position": [], "state_pot": []} for x in range(1, replicas + 1)
        }  # init transition dicts following one replica with inital start

        # go through repda and count
        tmp_run = 1
        for index, row in self.DATA_NEW.iterrows():
            if tmp_run != row.run:  # new trial
                transition_dict = tmp_dict
                tmp_dict = {x: x for x in range(1, replicas + 1)}
                tmp_run = row.run

            # Exchange Replica
            replica = int(transition_dict[int(row.ID)])  # get the replica unique id
            # record Exchange
            if row.s == 1:  # only hit when exchange and not partner already exchangeds
                # new_pos
                transition_result_dict[replica]["trial"].append(int(row.run))
                transition_result_dict[replica]["position"].append(int(row.partner))
                transition_result_dict[replica]["state_pot"].append(row.state_potentials)
                # exchange reps
                tmp_dict[int(row.partner)] = replica

            else:
                transition_result_dict[replica]["trial"].append(int(row.run))
                transition_result_dict[replica]["position"].append(int(row.ID))
                transition_result_dict[replica]["state_pot"].append(row.state_potentials)
                tmp_dict[int(row.ID)] = replica

            # if (replica == 2 and row.run < 10):
            #    print("trial ", row.run, "ID ", row.ID, "partner ", row.partner)
            #    print("replica: ", transition_dict[row.ID], "vs.", transition_dict[row.partner])
            #    print("transd", transition_result_dict[2]["position"])

        traces = {x: pd.DataFrame(transition_result_dict[x]) for x in transition_result_dict}
        [df.insert(0, "replicaID", replicaID) for replicaID, df in traces.items()]
        self.transition_traces = pd.concat(traces)

    def _calculate_ndowns_nups_for_each_state(
        self, time_stride: int = -1, min_state_potential_treshold: float = None, verbose: bool = False
    ):
        """_calculate_ndowns_nups_for_each_state

                        calculates the visit counts for each replicaID position (Temperature or s_value).
            It splits into substates depending on the state potentials, to destinguish which state passed by.
            Up and Downs are counted from top to bottom.
            Additionally also a time dependend series is generated. This can be binned by the  argument time_window.
            If the min_state_potential_treshold is given, than a minimal state is also dependent on the other states, if they are below the threshold, the state is undefined.

            In the end you recieve the position state visit counts in a dict:
            format: {{replicaposition:{"tot_nup":[], "tot_ndown":[], "dt":float, "dt_nup":[], "dt_ndown":[]}}

        Parameters
        ----------
        time_stride :  int, optional
            determines the window bin size of flow trajectory for each replicaID. This there are total_transitions/time_window bins containing time_window many flow values. default -1 counts all frames
        min_state_potential_treshold :  float, optional
             a threshold, defining if a state is governing a system at a time point t
        verbose :   bool, optional

        Returns
        -------
        None
        """

        # define needed stuff for calc:
        replica_traces = self.get_replica_traces()
        num_states = len(self.system.state_eir)
        num_replica = len(self.system.s)

        if time_stride < 1:
            time_stride = 1  # arbitrary window size value, that seems reasonable!  len(replica_traces[list(replica_traces.keys())[0]]["trial"]) * 0.01

        extreme_positions = (1, num_replica)  # gives the extreme values of the replicaID dist
        replica_extreme_position_memory = {
            replica: -1 for replica in range(1, num_replica + 1)
        }  # which replicaID visited which extreme?
        replica_extreme_position_memory.update(
            {1: extreme_positions[0], num_replica: extreme_positions[1]}
        )  # init pos1 and last one

        # result vars
        # for easier keeping track of state indices
        state_index = {key: key for key in range(num_states)}
        if min_state_potential_treshold is not None:  # add an undef state if multiple residues are below threshold.
            state_index.update({"undefined": num_states})
            num_states += 1

        count_state_perpos = {
            positionID: {
                "tot_nup": [0 for state in state_index],
                "tot_ndown": [0 for state in state_index],
                "dt": time_stride,
                "pot_tresh": min_state_potential_treshold,
                "states_index": state_index,
                "dt_nup": [[0 for state in state_index]],
                "dt_ndown": [[0 for state in state_index]],
            }
            for positionID in range(1, num_replica + 1)
        }

        if verbose:
            print("general: ", extreme_positions)
        if verbose:
            print("time_window_size: ", time_stride)
        # if verbose: print("INITIAL")
        # if verbose: print("Initial count_per_repPos\n", count_state_per_position)
        # if verbose: print("Initial current_extremePos_replica\n", replica_extreme_position_memory)

        # as side product easily the round trips can be calculated!
        replica_round_trips = {replica: -1 for replica in range(1, num_replica + 1)}
        replica_round_trips[1] = 0
        replica_round_trips[num_replica] = 0

        for index, (replicaID, trial, position, pot_energies) in replica_traces[::time_stride].iterrows():
            count_state_perpos[position]["dt_ndown"].append([0 for state in state_index])
            count_state_perpos[position]["dt_nup"].append([0 for state in state_index])

            if position in extreme_positions and replica_extreme_position_memory[replicaID] != position:
                replica_extreme_position_memory.update({replicaID: position})
                replica_round_trips[replicaID] += 1

            # This replicaID has already seen an extreme pos
            if replica_extreme_position_memory[replicaID] in extreme_positions:
                # who is the active state?
                if (
                    min_state_potential_treshold is not None
                ):  # NEW shall no other state be in an undersampling situation?
                    undersampling_state_energies = [
                        float(val) for val in list(pot_energies.values()) if (float(val) < min_state_potential_treshold)
                    ]
                    if 1 == len(
                        undersampling_state_energies
                    ):  # clean active states - only one state at a time underSampling
                        active_state = undersampling_state_energies.index(min(undersampling_state_energies))
                    else:  # no clear state presen skip
                        active_state = state_index["undefined"]
                else:
                    undersampling_state_energies = [
                        float(val) for val in list(pot_energies.values())
                    ]  # if(float(val) < 200)]
                    active_state = undersampling_state_energies.index(min(undersampling_state_energies))

                # determine if replicaID comes from top or bottom and add +1 to stat
                if replica_extreme_position_memory[replicaID] == extreme_positions[0]:  # coming from top
                    count_state_perpos[position]["tot_ndown"][active_state] += 1
                    count_state_perpos[position]["dt_ndown"][-1][active_state] += 1
                elif replica_extreme_position_memory[replicaID] == extreme_positions[1]:  # coming_from bottom
                    count_state_perpos[position]["tot_nup"][active_state] += 1
                    count_state_perpos[position]["dt_nup"][-1][active_state] += 1
                else:  # NEW has never seen any thing
                    raise ValueError("A replicaID has never seen a extreme position should not reach this code!")
            else:
                continue

        if verbose:
            print("\nFINAL")
        if verbose:
            print("Final extreme_positions", replica_extreme_position_memory)
        if verbose:
            print(
                "Final position counts totup/totdown: ",
                [
                    (count_state_perpos[pos]["tot_nup"], count_state_perpos[pos]["tot_ndown"])
                    for pos in count_state_perpos
                ],
            )
        if verbose:
            print("Final positoin counts keys", count_state_perpos[1].keys())
        if verbose:
            print("counted rountrips per replicaID!: ", replica_round_trips)

        # store trajs in pd.dataframes.
        tmp_counts = self._clean_replica_round_trips(replica_round_trips)
        for x, data in tmp_counts.items():
            column_names = ["state_" + str(data["states_index"][x]) for x in sorted(data["states_index"])]

            tmp = pd.DataFrame(data["dt_nup"])
            tmp.columns = column_names
            data["dt_nup"] = tmp

            tmp = pd.DataFrame(data["dt_ndown"])
            tmp.columns = column_names
            data["dt_ndown"] = tmp

        self.count_state_per_position = count_state_perpos
        self.replica_round_trips = tmp_counts

    def _calculate_replica_roundtrips(self):
        """
        ..autofunction: _calculate_replica_roundtrips
            This function is calculating the roundtrips over all replica positions for each replica.

        :return: None
        :rtype: None
        """
        # define needed stuff for calc:
        replica_traces = self.get_replica_traces()
        _ = len(self.system.state_eir)
        num_replica = len(self.system.s)

        extreme_positions = (1, num_replica)  # gives the extreme values of the replica dist
        replica_extreme_position_memory = {
            replica: -1 for replica in range(1, num_replica + 1)
        }  # which replica visited which extreme?
        replica_extreme_position_memory.update(
            {1: extreme_positions[0], num_replica: extreme_positions[1]}
        )  # init pos1 and last one

        # as side product easily the round trips can be calculated!
        replica_round_trips = {replica: -1 for replica in range(1, num_replica + 1)}
        replica_round_trips[1] = 0
        replica_round_trips[num_replica] = 0

        # only go over_extreme postitions.
        extreme_position_trace = replica_traces.loc[replica_traces.position.isin(extreme_positions)].sort_values(
            "trial"
        )
        for index, (replicaID, trial, position, pot_energies) in extreme_position_trace.sort_values(
            "trial"
        ).iterrows():  # go through each replica trace
            # print(trial, position, pot_energies)
            if position in extreme_positions and replica_extreme_position_memory[replicaID] != position:
                replica_extreme_position_memory.update({replicaID: position})
                replica_round_trips[replicaID] += 1
            else:
                continue
        self.replica_round_trips = self._clean_replica_round_trips(replica_round_trips)

    def append(self, repdat: Union[List[Repdat_Type], Repdat_Type]):
        """append

            This function concatenates two repdat files into the executing obj.

        Parameters
        ----------
        repdat :    List[Repdat] or Repdat
            one or multiple Repdat files.

        Returns
        -------
        None

        """
        # if(self.system != repdat.system):
        #    raise ValueError("The two repdats seem not to come from the same simulation, as the system settings are different!")

        if not isinstance(repdat, List):
            repdat = [repdat]

        self.DATA_NEW = pd.concat([self.DATA_NEW, *map(lambda x: x.data2, repdat)], ignore_index=True)
        self.DATA_NEW.run = pd.Series(map(lambda i: 1 + int(i) // len(self.system.s), self.DATA_NEW.index))

    def clean_file_runs(self):
        """clean_file
        DEAPRECEATED - ABOUT TO BE REMOVED!

                Updates the run numbers to be continous sequential. (for example needed for concatenation)
        """

        tmp_run = 0
        offset = 0
        for line in self.DATA:
            if int(tmp_run) > int(line.run):
                tmp_run = int(line.run)
            elif int(tmp_run) < int(line.run):
                tmp_run += line.run
                offset += 1
            line.run = offset

    def get_replica_traces(self, recalculate: bool = False) -> pd.DataFrame:
        """get_replica_traces
                        returns a replica_traces dictionary.
        Parameters
        ----------
        recalculate :   bool, optional
            shall the dict be recalculated, if already present?

        Returns
        -------
        Dict[int, Dict[str,List[float]]]
            dictionary containing all individual replica_traces
        """

        if not isinstance(self.transition_traces, pd.DataFrame) or recalculate:
            self._caculate_transition_traces()
        return self.transition_traces

    def get_replicaPosition_dependend_nup_ndown_for_each_state(
        self, time_window_size: int = -1, potential_treshold: float = None, recalculate: bool = False
    ) -> Dict[int, Dict[str, Union[List or float]]]:
        """
        ..autofunction: get_replicaPosition_dependend_nup_ndown_for_each_state
            This function is returning the replica position visit counts by each simulation state, per state.

        :param time_window_size: how many timesteps shall be binned
        :type time_window_size: int
        :param potential_treshold: if defined, and there is a time window, in which multiple states are below this threshold, the count is in an undefined state.
        :type potential_treshold: float
        :param recalculate: shall the dict be recalculated?
        :type recalculate: bool
        :return: returns a dict for all replica positions and their state visit counts.
        :rtype: Dict[int, Dict[str, Union[List or float]]]
        """
        if not isinstance(self.count_state_per_position, pd.DataFrame) or recalculate:
            self._calculate_ndowns_nups_for_each_state(
                time_stride=time_window_size, min_state_potential_treshold=potential_treshold
            )
        else:
            if not all(
                [
                    self.count_state_per_position[1]["dt"] == time_window_size,
                    self.count_state_per_position[1]["pot_tresh"] == potential_treshold,
                ]
            ):
                self._calculate_ndowns_nups_for_each_state(
                    time_stride=time_window_size, min_state_potential_treshold=potential_treshold
                )
        return self.count_state_per_position

    def get_replicaPosition_dependend_nup_ndown(
        self, time_window_size: int = -1, potential_treshold: float = None, recalculate: bool = False
    ) -> Dict[int, Dict[str, Union[List, pd.DataFrame, dict, float]]]:
        """
        ..autofunction: get_replicaPosition_dependend_nup_ndown
            This function is returning the replica position visit counts by all simulation state.

        :param time_window_size: how many timesteps shall be binned
        :type time_window_size: int
        :param potential_treshold: if defined, and there is a time window, in which multiple states are below this threshold, the count is in an undefined state.
        :type potential_treshold: float
        :param recalculate: shall the dict be recalculated?
        :type recalculate: bool
        :return: returns a dict for all replica positions and the visit counts.
        :rtype: Dict[int, Dict[str, Union[List or float]]]
        """
        if not isinstance(self.replicas_pos_visit_counts, Dict):
            replicas_pos_visit_counts = {}
            for replica, statistics in self.get_replicaPosition_dependend_nup_ndown_for_each_state(
                time_window_size=time_window_size, potential_treshold=potential_treshold, recalculate=recalculate
            ).items():
                replica_pos_visit_counts = {
                    replica: {
                        "tot_nup": sum(statistics["tot_nup"]),
                        "tot_ndown": sum(statistics["tot_ndown"]),
                        "dt": statistics["dt"],
                        "dt_nup": list(map(lambda x: sum(x), statistics["dt_nup"])),
                        "dt_ndown": list(map(lambda x: sum(x), statistics["dt_ndown"])),
                    }
                }

                replicas_pos_visit_counts.update(replica_pos_visit_counts)
        return replicas_pos_visit_counts

    def get_replica_roundtrips(self, recalculate: bool = False) -> Dict[int, int]:
        """
        ..autofunction: get_replica_roundtrips
            This function is returning the count of rountrips (RT) for each replica.
        :param recalculate: shall the dict be recalculated?
        :type recalculate: bool
        :return: returns a dict for all replica and their rountrip counts.
        :rtype: Dict[int,int]
        """
        if not isinstance(self.replica_round_trips, pd.DataFrame) or recalculate:
            self._calculate_replica_roundtrips()
        return self.replica_round_trips

    def write(self, out_path: str) -> str:
        """
        ..autofunction:
            Write out a repdat file to the outpath.
        :param out_path:    determines the output path
        :type out_path: str
        :return: out_path
        :rtype:str
        """
        # content_dict{"system":..., "header":..., "data":....}
        file = open(out_path, "w")
        file.write(str(self.system))

        file.write("\n\n")
        # print header:
        file.write(
            "#"
            + "\t".join(
                [
                    "ID",
                    "partner",
                    "run",
                    "li",
                    "Ti",
                    "Epoti",
                    "lj",
                    "Tj",
                    "Epotj",
                    "p",
                    "s",
                    "si",
                    "sj",
                    " ".join(["V" + str(i) for i in range(1, len(self.system.state_eir) + 1)]),
                ]
            )
            + "\n"
        )

        # file.write("\n#"+"\t".join(self.content["header"])+"\n")
        for line in self.DATA:
            file.write(str(line))
        file.close()
        return out_path
