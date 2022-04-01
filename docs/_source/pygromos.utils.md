---
date: '2022-04-01T13:55:18.499Z'
docname: _source/pygromos.utils
images: {}
path: /source-pygromos-utils
title: pygromos.utils package
---

# pygromos.utils package

## Submodules

## pygromos.utils.amino_acids module

Amino Acid Library

### Notes

This module contains all natural amino acids as named tuples.
Also contained are dictionaries for three letter code and so on.
Todo: add number of atoms in full atomistic case and other useful properties.


### _class_ pygromos.utils.amino_acids.amino_acid(name, oneLetter, threeLetter, numUnitedAtoms, numFullAtomistic)
Bases: `tuple`


#### name()
Alias for field number 0


#### numFullAtomistic()
Alias for field number 4


#### numUnitedAtoms()
Alias for field number 3


#### oneLetter()
Alias for field number 1


#### threeLetter()
Alias for field number 2

## pygromos.utils.bash module

bash-wrapper

### Description

> This lib contains bash-wrappers for general use in a workflow, with nice error messages and other features.
> Additionally there are some convenience functions for error managment and dependencies check


* **author**

    Benjamin Schroeder



### pygromos.utils.bash.check_path_dependencies(check_required_paths: Union[Dict[any, str], List[str]], check_warn_paths: Union[str, List[str]] = [], verbose: bool = True)
> checks a list of dependencies if each path is present or not. throws an
> IOError, if an Path is not existing.


* **Parameters**

    
    * **check_required_paths** (*Union[t.Dict[any, str], List[str]]*) – if path does not exist, an error wil be raised


    * **check_warn_paths** (*Union[t.Dict[any, str], List[str]]*) – if path does not exist, a warning will be written out.


    * **verbose** (*bool, optional*)



* **Raises**

    **IOERROR** – if a file is not existing



* **Returns**

    in_dependencies



* **Return type**

    Union[t.Dict[any, str], List[str]]



### pygromos.utils.bash.command_exists(command: str)
Does the command exists in the current system / path?

Args:

    command (str): Name of the command.

Returns:

    bool: Does command exists or not.


### pygromos.utils.bash.compress_gzip(in_path: str, out_path: Optional[str] = None, extract: bool = False, force: bool = True, verbose: bool = False)
> compress a file or directory with tar via the OS


* **Parameters**

    
    * **in_path** (*str*)


    * **out_path** (*str*)


    * **gunzip_compression** (*bool, optional*)


    * **remove_in_file_afterwards** (*bool, optional*)


    * **verbose** (*bool, optional*)



* **Raises**

    **OSError** – if the bas command fails an OSError will be raised.



* **Returns**

    
    * *str* – out_path


    * *FileIO* – process return log.




### pygromos.utils.bash.compress_tar(in_path: str, out_path: Optional[str] = None, gunzip_compression: bool = False, remove_in_file_afterwards: bool = False, remove_in_dir_afterwards: bool = False, verbose: bool = False)
> compress a file or directory with tar via the OS


* **Parameters**

    
    * **in_path** (*str*)


    * **out_path** (*str*)


    * **gunzip_compression** (*bool, optional*)


    * **remove_in_file_afterwards** (*bool, optional*)


    * **verbose** (*bool, optional*)



* **Raises**

    **OSError** – if the bas command fails an OSError will be raised.



* **Returns**

    
    * *str* – out_path


    * *FileIO* – process return log.




### pygromos.utils.bash.concatenate_text_files(in_file_paths: List[str], out_file_path: str, verbose: bool = False)
> concatenated multiple text files via the OS
> TODO: refactor to use execute command!


* **Parameters**

    
    * **in_file_paths** (*List[str]*)


    * **out_file_path** (*str*)


    * **verbose** (*bool, optional*)



* **Raises**

    
    * **OSError** – if the bash command fails an OSError will be raised.


    * **IOError** – if the the directory for the output file cann  not be found or not all in_paths exist.



* **Return type**

    None



### pygromos.utils.bash.copy_file(in_file_path: str, out_file_path: str, copy_a_directory: bool = False, additional_option: str = '')
> copy files via the OS
> TODO: refactor to use execute command!


* **Parameters**

    
    * **in_file_path** (*str*)


    * **out_file_path** (*str*)


    * **copy_a_directory** (*bool, optional*)


    * **additional_option** (*str, optional*)



* **Raises**

    **OSError** – if the bas command fails an OSError will be raised.



* **Returns**

    out_file_path



* **Return type**

    str



### pygromos.utils.bash.directory_exists(path: str)
Tests whether the provided path is valid and
also is a directory. Returns false if either condition
is not fullfilled.

Args:

    path (str): Path of the file or directory.

Returns:

    bool: Is a directory with a valid path or not.


### pygromos.utils.bash.execute(command: Union[str, List[str]], verbose: bool = False, catch_STD: Union[bool, str] = False, env: Optional[dict] = None)

### pygromos.utils.bash.execute_old(command: Union[str, List[str]], verbose: bool = False, ignore_return_code: bool = False, wait_fail=False, out_cnf_path: Optional[str] = None)
execute

    this command executes a command on the os-layer. (e.g.: on linux a bash command)


* **Parameters**

    
    * **command** (*str, Lists[str]*)


    * **verbose** (*bool, optional*)



* **Raises**

    **OSError** – if bash command fails an OSError is raised



* **Returns**

    returns the execution log of the process.



* **Return type**

    io.FileIO



### pygromos.utils.bash.execute_os(command: Union[str, List[str]], verbose: bool = False)
execute

> DEAPRECIATED

> sadly I can only use os.popen, as otherwise euler refuses.
> this command executes a command on the os-layer. (e.g.: on linux a bash command)
> therefore it uses os.popen


* **Parameters**

    
    * **command** (*str, Lists[str]*)


    * **verbose** (*bool, optional*)



* **Raises**

    **OSError** – if bash command fails an OSError is raised



* **Returns**

    returns the execution log of the process.



* **Return type**

    io.FileIO



### pygromos.utils.bash.execute_subprocess(command: Union[str, List[str]], catch_STD: Union[bool, str] = False, env: Optional[dict] = None, verbose: bool = False)
This command starts a subprocess, that is executing the str command in bash.


* **Parameters**

    
    * **command** (*str*) – bash command


    * **catch_STD** – if bool: catch the output and past it into the command line
    if str: catch output and write it into this file


    * **env** (*dict*) – environment


    * **verbose** (*bool*) – verbosity



* **Returns**

    return the executed process obj. (from subprocess)



* **Return type**

    CompletedProcess



### pygromos.utils.bash.extract_tar(in_path: str, out_path: Optional[str] = None, gunzip_compression: bool = False, remove_tar_afterwards: bool = False, verbose: bool = False)
> this wrapper helps you to unpack a tar file via the OS
> – since I always used to forget the command…


* **Parameters**

    
    * **in_path** (*str*)


    * **out_path** (*str*) – optional output, to which the untared in is moved.


    * **gunzip_compression** (*bool, optional*) – is the tar file in the gunzip compression format? (-z added)


    * **remove_tar_afterwards** (*bool, optional*)


    * **verbose** (*bool, optional*)



* **Raises**

    **OSError** – if the bas command fails an OSError will be raised.



* **Returns**

    
    * *str* – out_path


    * *FileIO* – process return log.




### pygromos.utils.bash.file_exists(path: str)
Tests whether the provided path is valid and
also is a file. Returns false if either condition
is not fullfilled.

Args:

    path (str): Path of the file or directory.

Returns:

    bool: Is a file with a valid path or not


### pygromos.utils.bash.is_directory(path: str)
Is the provided path a directory.

Args:

    path (str): Path of the file or directory.

Returns:

    bool: Is a directory or not.


### pygromos.utils.bash.is_file(path: str)
Is the provided path a file.

Args:

    path (str): Path of the file or directory.

Returns:

    bool: Is a file or not.


### pygromos.utils.bash.link_folder(in_directory_path: str, out_link_path: str, additional_options: str = '')
> linking a folder via the OS


* **Parameters**

    
    * **in_directory_path** (*str*)


    * **out_link_path** (*str*)


    * **additional_options** (*str, optional*)



* **Raises**

    
    * **IOError** – If already a file with the link exists.


    * **OSError** – If the bash command failed.



* **Returns**

    out_link_path



* **Return type**

    str



### pygromos.utils.bash.make_folder(in_directory_path: str, additional_option: str = '', verbose: bool = False)
> make a new folder via the os.


* **Parameters**

    
    * **in_directory_path** (*str*)


    * **additional_option** (*str, optional*)


    * **verbose** (*bool, optional*)



* **Raises**

    **OSError** – If the bash command failed.


**WARNING**: if folder already existed, a warning will be printed.


* **Returns**

    directory_path



* **Return type**

    str



### pygromos.utils.bash.move_file(in_file_path: str, out_file_path: str, additional_options: str = '', verbose: bool = False)
> copy files via the OS
> TODO: refactor to use execute command!


* **Parameters**

    
    * **in_file_path** (*str*)


    * **out_file_path** (*str*)


    * **additional_options** (*str, optional*)



* **Raises**

    **OSError** – if the bash command fails an OSError will be raised.



* **Returns**

    out_file_path



* **Return type**

    str



### pygromos.utils.bash.path_exists(path: str)
Does the provided path exists? Gives no information on whether
it is a directory or a file.

Args:

    path (str): Path of the file or directory.

Returns:

    bool: Does path exists or not.


### pygromos.utils.bash.remove_file(in_file_path: str, recursive: bool = False, additional_options: str = '')
> delete a file via the OS.
> TODO: refactor to use execute command!


* **Parameters**

    
    * **in_file_path** (*str*)


    * **recursive** (*bool, optional*) – remove recursivly! CAREFULL is DANGEROUS!


    * **additional_options** (*str, optional*)



* **Raises**

    **OSError** – if the bash command fails an OSError will be raised.



* **Return type**

    None



### pygromos.utils.bash.remove_folder(in_directory_path: str, additional_options: str = '', verbose: bool = False)
> make a new folder via the os.


* **Parameters**

    
    * **in_directory_path** (*str*)


    * **additional_options** (*str, optional*)


    * **verbose** (*bool, optional*)



* **Raises**

    **OSError** – If the bash command failed.


**WARNING**: if folder did not exist, a warning will be printed.


* **Return type**

    None



### pygromos.utils.bash.replace_text_in_text_file(in_file_path: str, find_pattern: str, replace_pattern: str, out_file_path: Optional[str] = None)
> this file replaces a regex pattern in a text file


* **Parameters**

    
    * **in_file_path** (*str*)


    * **find_pattern** (*str*)


    * **replace_pattern** (*str*)


    * **out_file_path** (*str, optional*) – if no out_file_path given, the text is replaced in place.



* **Raises**

    **OSError** – If the bash command failed.



* **Returns**

    out_file



* **Return type**

    str



### pygromos.utils.bash.save_make_folder(in_directory_path: str, additional_options: str = '')
This function checks if a folder already exists and adds an integer to it, so that there will be a guaranteed unique folder.


* **Parameters**

    
    * **in_directory_path** (*str*)


    * **additional_options** (*str, optional*)



* **Raises**

    **OSError** – If the bash command failed.



* **Returns**

    in_directory_path



* **Return type**

    str



### pygromos.utils.bash.wait_for_fileSystem(check_paths: Union[str, List[str]], regex_mode: bool = False, max_waiting_iterations: int = 1000, verbose: bool = False)
This function can be used to circumvent lsf lag times.


* **Parameters**

    
    * **check_path** (*str - Path to file to check if existant*)


    * **max_waiting_iterations** (*int - maximal iteration Time*)


    * **verbose** (*bool - tell me if found*)



* **Returns**

    on success



* **Return type**

    True



### pygromos.utils.bash.which(command: str)
Finds the full path of a command.

Args:

    command (str): Name of the command.

Returns:

    str: Full path of the command.

## pygromos.utils.compiledProgram module

This abstract parent class, should provide functionality for  checking if binaries are present.


### _class_ pygromos.utils.compiledProgram._compiled_program(in_bin_dir: Optional[str], _force_bin_present: bool = True, _check_binary_paths: bool = True)
Bases: `object`

This class contains functions, that are required to detect if a binary program is accesible.
It tracks the presence of the binary dir and wraps the gromos functionality with wrappers, that check on function call if the binary is present.


#### \__init__(in_bin_dir: Optional[str], _force_bin_present: bool = True, _check_binary_paths: bool = True)
> The  _compiled_program parent class can be used, to ensure on runtime, that certain binaries are present.


* **Parameters**

    
    * **in_bin_dir** (*Union[str, None]*) – directory that should contain the binaries. If None, the assumption is made, that the path is part of the PATH variable.


    * **_force_bin_present** (*bool, optional*) – if True, the check_binary or check_binary_folder will throw errors, if they don’t find the targets, by default True


    * **_dont_check_binary** (*bool, optional*) – This is a kill switch for all checks of this class., by default False



* **Returns**

    _description_



* **Return type**

    Union[str, None]



#### _check_all_binaries(_force_bin_present: bool = False)
This function checks all present programs of this class, if the binary can be found and executed.
It does not trigger an Exception, if a binary cannot be found, except force_present is True.
However it records, if a binary not found in the _function_binary dict.
The recorded information is used to add restrictive wrappers to the progams, that can throw errors on executions (if _dont_check_bin attribute is False)


* **Parameters**

    **force_present** (*bool, optional*) – the function can be restrictive and throw errors, if a function was not found, by default False



* **Returns**

    all binaries present?



* **Return type**

    bool



#### _check_binaries_decorator(func: Callable)
> This function wrapper adds a binary check before, the function is executed.


* **Parameters**

    **func** (*callable*) – function using a binary and having the parameter _binary_name:str



* **Returns**

    wrapped function



* **Return type**

    callable



* **Raises**

    **Exception** – if no binary can be found in the function call or definition



#### _check_binary(test_program: str)
> this function checks if the binary exists in the given binary dir. (if binary_dir == None -> it will check the PATH variable.)


* **Parameters**

    **test_program** (*str*) – name of the binary.



* **Returns**

    if the binary was found



* **Return type**

    bool



* **Raises**

    **IOError** – If the binary dir was not found and _dont_check_bin was False (default: False)



#### _check_binary_dir(in_bin_dir: str)
> this fuction checks and records if the provided binary dir was found.


* **Parameters**

    **in_bin_dir** (*str*) – the expected binary dir.



* **Returns**

    Does the binary dir exist?



* **Return type**

    bool



* **Raises**

    **IOError** – If the binary dir was not found and _dont_check_bin was False (default: False)



#### _compiled_program__wrap_programms_with_binary_checks(remove=False)
> Wraps all functions with the __check_binaries_decorator


* **Parameters**

    **remove** (*bool, optional*) – remove all wrappers?, by default False



#### _property_ bin(_: st_ )
This attribute is the binary directory.

## pygromos.utils.pdb module

PDB
.. warning:: This module is very specific and not recomended to be used straigth away!

### Notes

contains some translation and reordering scripts for some pdbs.


### pygromos.utils.pdb.check_ATOM_columns(lines: list) -> (<class 'dict'>, <class 'list'>)

### pygromos.utils.pdb.consecutivley_renumber(header: dict, lines: list)

### pygromos.utils.pdb.filter_atoms_from_residue(filter_out: str, residue_name: str, lines: list) -> (<class 'list'>, <class 'list'>)

### pygromos.utils.pdb.form_columns(header: dict, lines: list)

### pygromos.utils.pdb.read_pdb_simple(path) -> (<class 'list'>, <class 'list'>, <class 'list'>)

### pygromos.utils.pdb.rename_atom_attribute(pattern: str, replace: str, lines: list)

### pygromos.utils.pdb.rename_atom_types(lines, header, in_ff=False, out_ff=False, translation_unknown=False)

### pygromos.utils.pdb.reorder_lines(header: dict, lines: list)
order_dict = {  # charmm_gromos54

    “DPPC”: {

        “N”: 4, # Choline “C11”: 1, “C14”: 2, “C15”: 3, “C12”: 5, “C11”: 6,
        “P”: 8, # Phosphate “O32”:7, “O33”:9, “O34”: 10, “O31”: 11, “C3”:
        12, # Glycerol “C2”: 13, “O21”: 14, “C21”: 15, “O22”: 16, “C22”: 17,
        “C1”: 32, “O11”: 33, “C11”: 34, “O12”: 35, “C12”: 36,

        “C23”, “C24”,”: 18, # Fatty acid chain1 “C25”,”: 19, “C26”,”: 20,
        “C27”,”: 21, “C28”,”: 22, “C29”,”: 23, “C210”: 24, “C211”: 25,
        “C212”: 26, “C213”: 27, “C214”: 28, “C215”: 29, “C216”: 30,

        “C13”: 31, “C14”, # Fatty acid chain2 “C15”, “C2C”: 37, “C16”,
        “C2D”: 38, “C17”, “C2E”: 39, “C18”, “C2F”: 40, “C19”, “C2G”: 41,
        “C110””C2H”: 42, “C111””C2I”: 43, “C112””C2J”: 44, “C113””C2K”: 45,
        “C114””C2L”: 46, “C115””C2M”: 47, “C116””C2N”: 48, “C2O”: 49, “C2P”:
        50},

    “SOLV”: {

    > > “OW”:1,

    > “HW1”:2, “HW2”:3

    }}

new_lines = [”” for x in range(len(lines))] residue = 0 offset = 0

for index.rst, line in enumerate(lines):

    res_index = header[“resnum”] resname_index = header[“resname”]
    atomname_index = header[“atomname”] atomname = line[atomname_index]
    resname = line[resname_index]

    if(line[res_index]!=residue):

        residue = line[res_index] offset = index.rst
        #int(line[header[“atomnum”]])-1

    if(resname in order_dict):

        relativ_pos = order_dict[resname][atomname]
        new_lines[relativ_pos+offset] = line

    else:

    > new_lines[offset] = line

#print(new_lines[len(new_lines)-1]) return new_lines

Args:

    header (dict):
    lines (list):

## pygromos.utils.typing module

## pygromos.utils.utils module

utils

This module should contain small usefull functions.


### pygromos.utils.utils._inline_dict(in_dict: Dict, prefix: str = '\\t')
> translate dictionary to one code line. can be used for meta-scripting


* **Parameters**

    
    * **in_dict** (*Dict*) – analysis control dict


    * **prefix** (*str, optional*) – prfix symbol to dict write out.



* **Returns**

    code line.



* **Return type**

    str



### pygromos.utils.utils.dict_to_nice_string(control_dict: Dict)
> Converts a dictionary of options (like template_control_dict)

>     to a more human readable format. Which can then be printed to a text file,
>     which can be manually modified before submiting analysis jobs.


* **Parameters**

    **control_dict** (*Dict*) – analysis control dictonary



* **Returns**

    nice formatting of the control dictionary for printing.



* **Return type**

    str



### pygromos.utils.utils.dynamic_parser(func: callable, title: str)
> This function builds dynamically a parser obj for any function, that has parameters with annotated types.
> Result is beeing able to parse any function dynamically via bash.


* **Parameters**

    
    * **func** (*callable*) – the function that should be parsed


    * **title** (*str*) – title for the parser



* **Returns**

    The parsed arguments



* **Return type**

    args



* **Raises**

    **IOError** – error if a parsing arg is unknown



### pygromos.utils.utils.nice_s_vals(svals: List[float])
> This helper function formats s-vals for RE-EDS to nice readable values.
> Is/was used in RE-EDS applications. Main functionality is rounding the different s-vals to their significance digits.


* **Parameters**

    **svals** (*List[float]*) – smoothing parameters of a reeds approach



* **Returns**

    nicely rounded s-values, such that only significant digits for each number is around.



* **Return type**

    List[float]



### pygromos.utils.utils.str2bool(v)

### pygromos.utils.utils.write_job_script(out_script_path: str, target_function: callable, variable_dict: dict, python_cmd: str = 'python3', verbose: bool = False)
> this function writes submission commands into a file. The command will be started from a bash env into python.


* **Parameters**

    
    * **out_script_path** (*str*) – path of the output script.


    * **target_function** (*callable*) – the function, that shall be submitted


    * **variable_dict** (*dict*) – variables for this function


    * **python_cmd** (*str, optional*) – which python command shall be supplied


    * **verbose** (*bool, optional*) – c’est la vie



* **Returns**

    returns an out script path.



* **Return type**

    str



* **Raises**

    
    * **IOERROR** – if outpath is not possible


    * **ValueError** – if required variable from the var-dict for the function is missing


## Module contents
