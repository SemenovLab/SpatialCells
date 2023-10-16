from pkg_resources import get_distribution, DistributionNotFound

try:
    __version__ = get_distribution("spatialcells").version
except DistributionNotFound:
    __version__ = "(local)"


from . import preprocessing as prep
from . import spatial as spa
from . import measurements as msmt
from . import plotting as plt
from . import utils
