import logging
import csv

from multiqc.plots import violin

from ._helpers import group_median_by_cell_prefix

# Initialise the logger
log = logging.getLogger(__name__)


class scCountQCMixin:
    def parse_scCountQC(self):
        """Find scCountQC output."""
        self.sincei_scCountQC = dict()
        for f in self.find_log_files("sincei/scCountQC", filehandles=True):
            parsed_data = self.parsescCountQCFile(f)
            for k, v in parsed_data.items():
                if k in self.sincei_scCountQC:
                    log.warning(f"Replacing duplicate sample {k}.")
                self.sincei_scCountQC[k] = v
                # Superfluous function call to confirm that it is used in this module
                # Replace None with actual version if it is available
                self.add_software_version(None)
            if len(parsed_data) > 0:
                self.add_data_source(f, section="scCountQC")

        self.sincei_scCountQC = self.ignore_samples(self.sincei_scCountQC)

        if len(self.sincei_scCountQC) > 0:
            # Write data to file
            self.write_data_file(self.sincei_scCountQC, "sincei_count_qc")

            header = dict()
            header["n_genes"] = {
                "title": "# Features",
                "description": "No. of detected features (bins or genes) with non-zero counts (Median of cells)",
                "scale": "RdBu",
                "min": 0,
            }
            header["n_counts_log"] = {
                "title": "# Counts (log1p)",
                "description": "Total counts in features per cell (logN+1 scale, Median of cells)",
                "scale": "OrRd",
                "min": 0,
            }
            header["pct_50"] = {
                "title": "% Counts top 50",
                "description": "Percent of alignments in top 50 features (Median of cells)",
                "scale": "RdYlBu-rev",
                "min": 0,
                "max": 100,
            }
            header["pct_100"] = {
                "title": "% Counts top 100",
                "description": "Percent of alignments in top 100 features (Median of cells)",
                "scale": "BrBG-rev",
                "min": 0,
                "max": 100,
            }
            header["pct_200"] = {
                "title": "% Counts top 200",
                "description": "Percent of alignments in top 200 features (Median of cells)",
                "scale": "RdYlBu-rev",
                "min": 0,
                "max": 100,
            }
            header["pct_500"] = {
                "title": "% Counts top 500",
                "description": "Percent of alignments in top 500 features (Median of cells)",
                "scale": "PuOr-rev",
                "min": 0,
                "max": 100,
            }
            header["gini_coefficient"] = {
                "title": "Gini Coefficient",
                "description": "Gini coefficient of enrichment (inequality) of counts in features (Median of cells)",
                "scale": "BrBG",
                "min": 0,
                "max": 1,
            }
            kys = [
                "n_genes_by_counts",
                "log1p_total_counts",
                "pct_counts_in_top_50_genes",
                "pct_counts_in_top_100_genes",
                "pct_counts_in_top_200_genes",
                "pct_counts_in_top_500_genes",
                "gini_coefficient",
            ]

            test_dict = group_median_by_cell_prefix(self.sincei_scCountQC, kys[0])
            out = {}
            for k in test_dict.keys():
                out[k] = dict.fromkeys(kys)
                for p in kys:
                    dv = group_median_by_cell_prefix(self.sincei_scCountQC, p)
                    out[k].update(dv[k])

            tdata = dict()
            for k, v in out.items():
                tdata[k] = {
                    "SampleName": k,
                    "n_genes": v["n_genes_by_counts"],
                    "n_counts_log": v["log1p_total_counts"],
                    "pct_50": v["pct_counts_in_top_50_genes"],
                    "pct_100": v["pct_counts_in_top_100_genes"],
                    "pct_200": v["pct_counts_in_top_200_genes"],
                    "pct_500": v["pct_counts_in_top_500_genes"],
                    "gini_coefficient": v["gini_coefficient"],
                }

            config = {
                "id": "sincei-scCountQC-plot",
                "title": "sincei: scCountQC counting metrics",
                "namespace": "sincei scCountQC",
            }
            self.add_section(
                name="Counting Metrics",
                anchor="scCountQC",
                description="Statistics of distribution of counts per cells after counting using `scCountQC`",
                plot=violin.plot(tdata, header, config),
            )

        return len(self.sincei_scCountQC)

    def parsescCountQCFile(self, f):
        reader = csv.DictReader(f["f"], delimiter="\t")
        required = {"barcodes", "sample", "n_genes_by_counts", "log1p_n_genes_by_counts"}
        if required.difference(set(reader.fieldnames)):
            log.warning(
                f"{f['fn']} was initially flagged as the tabular output from scCountQC, but that seems to not be the case. Skipping..."
            )
            return dict()

        d = {}
        for row in reader:
            s_name = self.clean_s_name(row["Cell_ID"], f)
            if s_name in d:
                log.debug(f"Replacing duplicate cell_id {s_name}.")
            d[s_name] = {key: row[key] for key in reader.fieldnames}
        return d
