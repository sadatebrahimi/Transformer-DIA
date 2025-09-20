"""Mass spectrometry file type input/output operations."""
import collections
import csv
import os
import re
from pathlib import Path
from typing import Any, Dict

from .. import __version__


class MztabWriter:
    """
    Export spectrum identifications to an mzTab file.

    Parameters
    ----------
    filename : str
        The name of the mzTab file.
    """

    def __init__(self, filename: str):
        self.filename = filename

        self.psms = []

    def save(self) -> None:
        """
        Export the spectrum identifications to the mzTab file.
        """
        with open(self.filename, "w") as f:
            writer = csv.writer(f, delimiter="\t", lineterminator=os.linesep)

            writer.writerow(
                [
                    "PSH",
                    "sequence",
                    "PSM_ID",
                    "accession",
                    "unique",
                    "database",
                    "database_version",
                    "search_engine",
                    "search_engine_score[1]",
                    "modifications",
                    "retention_time",
                    "charge",
                    "exp_mass_to_charge",
                    "calc_mass_to_charge",
                    "spectra_ref",
                    "pre",
                    "post",
                    "start",
                    "end",
                    "opt_ms_run[1]_aa_scores",
                ]
            )
            for psm in self.psms:
                writer.writerow(
                    [
                        "PSM",
                        psm[0],  # sequence
                        str(psm[1]),  # PSM_ID
                        "null",  # accession
                        "null",  # unique
                        "null",  # database
                        "null",  # database_version
                        f"[MS, MS, Casanovo_DIA, {__version__}]",
                        psm[2],  # search_engine_score[1]
                        # FIXME: Modifications should be specified as
                        #  controlled vocabulary terms.
                        "null",  # modifications
                        # FIXME: Can we get the retention time from the data
                        #  loader?
                        "null",  # retention_time
                        psm[3],  # charge
                        psm[4],  # exp_mass_to_charge
                        psm[5],  # calc_mass_to_charge
                        f"ms_run[1]:index={psm[1]}",  # spectra_ref
                        "null",  # pre
                        "null",  # post
                        "null",  # start
                        "null",  # end
                        psm[6],  # opt_ms_run[1]_aa_scores
                    ]
                )
