from importlib.metadata import version

from .qc import qc

__all__ = ["qc"]
__version__ = version("barcodeqc")
