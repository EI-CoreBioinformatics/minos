import pkg_resources

__title__ = "minos"
__author__ = "Christian Schudoma (cschu), Gemy Kaithakottil"
__email__ = "gemy.kaithakottil@earlham.ac.uk"
__license__ = "MIT"
__copyright__ = "Copyright 2019-2023 Earlham Institute"
__version__ = pkg_resources.require("minos")[0].version

DEFAULT_FAKE_PROT = pkg_resources.resource_filename(
    "minos.etc", "fake_prot.fasta"
)
DEFAULT_FAKE_NUC = pkg_resources.resource_filename(
    "minos.etc", "fake_nuc.fasta"
)
