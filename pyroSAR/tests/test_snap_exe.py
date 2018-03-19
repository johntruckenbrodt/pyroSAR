import pytest
from contextlib import contextmanager
from pyroSAR._dev_config import ExamineExe
from pyroSAR.snap.auxil import ExamineSnap

@contextmanager
def not_raises(ExpectedException):
    try:
        yield

    except ExpectedException, err:
        raise AssertionError(
            "Did raise exception {0} when it should not!".format(
                repr(ExpectedException)
            )
        )

    except Exception, err:
        raise AssertionError(
            "An unexpected exception {0} raised.".format(repr(err))
        )

class TestExemineExe:
    def test_exception(self):
        with pytest.raises(ValueError):
            ExamineExe.examine('some_exe_file.exe')

    # def test_not_exception(self):
    #     SNAP_EXECUTABLE = ['snap64.exe', 'snap32.exe', 'snap.exe', 'snap']
    #     with not_raises(ValueError):
    #         ExamineExe.examine(SNAP_EXECUTABLE)

# class TestExamineSnap:
#     def test_not_exception(self):
#         with not_raises(AssertionError):
#             test_snap_exe = ExamineSnap()