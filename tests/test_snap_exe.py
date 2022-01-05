from contextlib import contextmanager


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
