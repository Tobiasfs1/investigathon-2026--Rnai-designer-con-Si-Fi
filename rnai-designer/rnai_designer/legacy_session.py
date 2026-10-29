"""Ephemeral per-workflow sharing; never shared across genes or invocations."""
from contextlib import contextmanager
from contextvars import ContextVar
from tempfile import TemporaryDirectory

shared_directory = ContextVar("sifi_shared_directory", default=None)


@contextmanager
def shared_sifi_session():
    with TemporaryDirectory(prefix="sf_") as directory:
        token = shared_directory.set(directory)
        try:
            yield
        finally:
            shared_directory.reset(token)
