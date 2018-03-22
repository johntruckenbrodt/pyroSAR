from contextlib import contextmanager

import pytest

from pyroSAR._dev_config import ExamineExe
from pyroSAR.snap.auxil import ExamineSnap


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


class TestExemineExe:
    def test_exception(self):
        with pytest.warns(UserWarning):
            ExamineExe.examine('some_exe_file.exe')

    def test_not_exception(self):
        SNAP_EXECUTABLE = ['snap64.exe', 'snap32.exe', 'snap.exe', 'snap']
        with pytest.warns(None) as record:
            ExamineExe.examine(SNAP_EXECUTABLE)
        assert len(record) == 0


class TestExamineSnap:
    def test_exception(self):
        with pytest.warns(UserWarning):
            ExamineSnap(snap_executable='some_exe_file.exe')

    def test_not_exception(self):
        with pytest.warns(None) as record:
            ExamineSnap()
        assert len(record) == 0
