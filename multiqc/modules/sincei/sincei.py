import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

# sincei modules
from .scCountQC import scCountQCMixin
from .scFilterStats import scFilterStatsMixin

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(
    BaseMultiqcModule,
    scFilterStatsMixin,
    scCountQCMixin,
):
    """
    sincei (short for Single Cell Informatics) is a command-line toolkit for
    exploration of single-cell epigenomics data. It accommodates data from a
    wide range of single-cell protocols, such as droplet-based (10x Genomics)
    and plate-based protocols, gene expression (scRNA-seq) and epigenomics
    (scATAC, scCUTnTAG, scBS-seq). sincei can be used for quality control of
    these datasets directly from BAM files (read-level QC), as well as after
    signal aggregation (count-level). Additional functionalities include
    filtering, dimensionality reduction, and plotting of single-cell data.

    The MultiQC module for sincei parses the following text outputs:

    - `scFilterStats` (default output file)
    - `scCountQC --outMetrics` (currently the cell-level metrics are supported)

    Each cell (identified by `Cell_ID` in sincei) is treated as a separate
    sample by MultiQC. Sample/cell names are parsed from the text files
    themselves rather than from file names.
    """

    def __init__(self):
        # Initialise the parent object
        super().__init__(
            name="sincei",
            anchor="sincei",
            target="sincei",
            href="https://sincei.readthedocs.io",
            info="Toolkit for processing and analyzing single-cell (epi)genomics data.",
            doi="10.5281/zenodo.7853375",
        )

        n = dict()

        # scFilterStats
        n["scFilterStats"] = self.parse_scFilterStats()
        if n["scFilterStats"] > 0:
            log.debug(f"Found {n['scFilterStats']} sincei scFilterStats reports")

        # scCountQC
        n["scCountQC"] = self.parse_scCountQC()
        if n["scCountQC"] > 0:
            log.debug(f"Found {n['scCountQC']} sincei scCountQC reports")

        tot = sum(n.values())
        if tot > 0:
            log.info(f"Found {tot} total sincei reports")
        else:
            raise ModuleNoSamplesFound
