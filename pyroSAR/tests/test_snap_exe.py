from contextlib import contextmanager

import pytest

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


class TestExamineSnap:
    def test_not_exception(self):
        with pytest.warns(None) as record:
            ExamineSnap()
        assert len(record) == 0
