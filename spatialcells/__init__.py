from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("spatialcells")
except PackageNotFoundError:
    __version__ = "(local)"


from . import preprocessing as prep
from . import spatial as spa
from . import measurements as msmt
from . import plotting as plt
from . import utils
