import os
from contextlib import contextmanager

import pytest
from os.path import expanduser

from pyroSAR._dev_config import ExamineExe, make_dir, write_config_file
from pyroSAR.snap.auxil import ExamineSnap, get_etc_from_config


@contextmanager
def not_raises(ExpectedException):
    try:
        yield

    except ExpectedException:
        raise AssertionError(
            "Did raise exception {0} when it should not!".format(
                repr(ExpectedException)
            )
        )

    except Exception:
        raise AssertionError(
            "An unexpected exception {0} raised.".format(repr(Exception))
        )


class TestExamineExe:
    def test_exception(self):
        with pytest.warns(UserWarning):
            ExamineExe.examine('some_exe_file.exe')

    def test_not_exception(self):
        SNAP_EXECUTABLE = ['snap64.exe', 'snap32.exe', 'snap.exe', 'snap']
        with pytest.warns(None) as record:
            ExamineExe.examine(SNAP_EXECUTABLE)
        # assert len(record) == 0


class TestExamineSnap:
    def test_exception(self):
        with pytest.warns(UserWarning):
            ExamineSnap(snap_executable='some_exe_file.exe')

    # Got some problems with TRAVIS CI because this function needs a input from user.
    def test_not_exception(self):
        with pytest.warns(None) as record:
            ExamineSnap()
        # assert len(record) == 0


class TestSnapExeAuxil:
    def test_make_dir(self):
        path = expanduser("~")
        make_dir(path)
        dir = os.listdir(path)
        assert '.pyrosar' in dir

    def test_write_config(self):
        path = os.path.join(expanduser("~"), '.pyrosar')
        write_config_file('test_data', header=None, path=path)

        dir = os.listdir(path)

        assert 'config.txt' in dir

    def read_config(self):
        path = expanduser("~")
        make_dir(path)

        write_config_file(data='test_path', header='snap_etc:', path=os.path.join(expanduser("~"), '.pyrosar'))

        content = get_etc_from_config()

        assert content == 'test_path'
