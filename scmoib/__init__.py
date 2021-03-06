from . import plotting as pl
from . import preprocessing as pp
from . import tools as tl
# from .utils import check_versions, annotate_doc_types
from ._version import get_versions  # version generated by versioneer

__author__ = ', '.join([
    'Atai Dobrynin',
    'Sheetal Giri',
    'Anna Danese',
])
__email__ = ', '.join([
    'atay.dobrynin@gmail.com',
    'anna.danese@helmholtz-muenchen.de',
    # We don’t need all, the main authors are sufficient.
])

__version__ = get_versions()['version']
