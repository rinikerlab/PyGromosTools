from pygromos.utils import bash
import os

# test for basic bash functionality


def test_which():
    # there should be an installation of bash
    cmd = "bash"
    expected_path_1 = "/bin/bash"
    expected_path_2 = "/usr/bin/bash"
    returned_path = bash.which(cmd)

    assert isinstance(returned_path, str)
    assert (expected_path_1 == returned_path) or (expected_path_2 == returned_path)


def test_command_exists():
    # there should be an installation of python
    cmd = "python"
    expected_command_exists = True
    returned_command_exists = bash.command_exists(cmd)

    assert isinstance(expected_command_exists, bool)
    assert expected_command_exists == returned_command_exists


def test_path_exists():
    expected_outcome = True
    returned_outcome = bash.path_exists(os.getcwd())

    assert isinstance(returned_outcome, bool)
    assert expected_outcome == returned_outcome


def test_is_directory():
    expected_outcome = True
    returned_outcome = bash.is_directory(os.getcwd())

    assert isinstance(returned_outcome, bool)
    assert expected_outcome == returned_outcome


def test_is_file():
    expected_outcome = True
    returned_outcome = bash.is_file(__file__)

    assert isinstance(returned_outcome, bool)
    assert expected_outcome == returned_outcome


def test_directory_exists():
    expected_outcome = True
    returned_outcome = bash.directory_exists(os.getcwd())

    assert isinstance(returned_outcome, bool)
    assert expected_outcome == returned_outcome


def test_file_exists():
    expected_outcome = True
    returned_outcome = bash.file_exists(__file__)

    assert isinstance(returned_outcome, bool)
    assert expected_outcome == returned_outcome
