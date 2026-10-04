"""specdb: Extragalactic Galaxy Spectra Database"""

from importlib.metadata import PackageNotFoundError, version as _version

try:
    __version__ = _version("specdb")
except PackageNotFoundError:  # not installed, e.g. running from a source tree
    __version__ = "unknown"
