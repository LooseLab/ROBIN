from cnsmeth.subpages.base_analysis import BaseAnalysis
import natsort
from cnsmeth import theme, resources
from cnsmeth import submodules
import pandas as pd
import numpy as np
import os,sys
import time
from nicegui import ui
from io import StringIO
import pysam
import shutil
os.environ["CI"] = "1"
import tempfile
HVPATH = os.path.join(
    os.path.dirname(os.path.abspath(submodules.__file__)), "hv_rapidCNS2"
)

class MGMT_Object(BaseAnalysis):
    def __init__(self, *args, **kwargs):
        self.MGMTbamfile = None
        self.counter = 0
        super().__init__(*args, **kwargs)

    def setup_ui(self):
        with ui.card().style("width: 100%"):
            ui.label("MGMT Methylation").tailwind("drop-shadow", "font-bold")
            self.mgmtplot = ui.row()
            with self.mgmtplot.style("width: 100%"):  # "size-full"):
                ui.label("Plot not yet available.")
            self.mgmtable = ui.row()
            with self.mgmtable:
                ui.label("Table not yet available.")

    def process_bam(self, bamfile, timestamp):
        MGMT_BED = f"{HVPATH}/bin/mgmt_hg38.bed"
        tempbamfile = tempfile.NamedTemporaryFile(suffix=".bam")
        os.system(
            f"bedtools intersect -a {bamfile} -b {MGMT_BED} > {tempbamfile.name}"
        )
        pysam.index(tempbamfile.name)
        if pysam.AlignmentFile(tempbamfile.name, "rb").count(until_eof=True) > 0:
            with self.mgmtplot:
                ui.notify("Running MGMT predictor - MGMT sites found.")
            if not self.MGMTbamfile:
                self.MGMTbamfile = tempfile.NamedTemporaryFile(suffix=".bam")
                shutil.copy2(tempbamfile.name, self.MGMTbamfile.name)
            else:
                tempbamholder = tempfile.NamedTemporaryFile(suffix=".bam")
                pysam.cat("-o",tempbamholder.name, self.MGMTbamfile.name, tempbamfile.name)
                shutil.copy2(tempbamholder.name, self.MGMTbamfile.name)
            tempmgmtdir = tempfile.TemporaryDirectory()
            pysam.sort("-o", os.path.join(tempmgmtdir.name, 'mgmt.bam'), self.MGMTbamfile.name)
            pysam.index(os.path.join(tempmgmtdir.name, 'mgmt.bam'))
            os.system(
                f"modkit pileup -t 4 --filter-threshold 0.73 --combine-mods {os.path.join(tempmgmtdir.name, 'mgmt.bam')} "
                f"{os.path.join(tempmgmtdir.name, 'mgmt.bed')} --suppress-progress  >/dev/null 2>&1 "
            )
            cmd = f"Rscript {HVPATH}/bin/mgmt_pred_v0.3.R --input={os.path.join(tempmgmtdir.name, 'mgmt.bed')} --out_dir={tempmgmtdir.name} --probes={HVPATH}/bin/mgmt_probes.Rdata --model={HVPATH}/bin/mgmt_137sites_mean_model.Rdata --sample=live_analysis"
            os.system(cmd)
            results = pd.read_csv(
                os.path.join(tempmgmtdir.name, "live_analysis_mgmt_status.csv")
            )
            self.counter+=1
            plot_out = f'{self.counter}_mgmt.png'
            os.system(
                f"methylartist locus -i chr10:129466536-129467536 -b {os.path.join(tempmgmtdir.name, 'mgmt.bam')} -o {plot_out}  --motif CG --mods m "
            )
            self.mgmtable.clear()
            with self.mgmtable:
                ui.table.from_pandas(results)

            if os.path.exists(plot_out):
                self.mgmtplot.clear()
                with self.mgmtplot.classes("w-full"):
                    ui.image(plot_out).props("fit=scale-down")
            tempmgmtdir.cleanup()
            with self.mgmtplot:
                ui.notify("MGMT predictor done.")
        self.running=False


def test_me():
    my_connection = None
    with theme.frame("Copy Number Variation Interactive", my_connection):
        TestObject = MGMT_Object(progress=True)
        #path = "tests/static/bam"
        path = "/users/mattloose/datasets/ds1305_Intraop0006_A/20231123_1233_P2S-00770-A_PAS59057_b1e841e7/bam_pass"
        directory = os.fsencode(path)
        for file in os.listdir(directory):
            filename = os.fsdecode(file)
            if filename.endswith(".bam"):
                TestObject.add_bam(os.path.join(path, filename))
    ui.run(port=12345)

# Entrypoint for when GUI is launched by the CLI.
# e.g.: python my_app/my_cli.py
if __name__ in {"__main__", "__mp_main__"}:
    """
    Entrypoint for when GUI is launched by the CLI
    :return: None
    """
    print("GUI launched by auto-reload")
    test_me()