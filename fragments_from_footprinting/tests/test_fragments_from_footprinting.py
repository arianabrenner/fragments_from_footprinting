"""
Unit and regression test for the fragments_from_footprinting package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import fragments_from_footprinting


def test_fragments_from_footprinting_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "fragments_from_footprinting" in sys.modules
