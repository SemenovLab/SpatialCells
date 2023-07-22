from pkg_resources import get_distribution, DistributionNotFound

try:
    __version__ = get_distribution('spacells').version
except DistributionNotFound:
    __version__ = '(local)'


from . import preprocessing as prep
from . import spatial as spa