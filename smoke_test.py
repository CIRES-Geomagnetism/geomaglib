# The smoke tests is for verifying deployment
import subprocess
import pytest
import os
@pytest.mark.smoke
def test_import():
    try:
        import geomaglib
        from geomaglib import util
        from geomaglib import sh_loader

        from geomaglib import sh_vars
        from geomaglib import magmath
        from geomaglib import legendre
    except ImportError as e:
        assert False, f"Import failed: {e}"