from cnsmeth.subpages.base_analysis import BaseAnalysis
from cnsmeth import theme
from cnsmeth import submodules
import pandas as pd
import os, sys
import asyncio
from nicegui import ui, run
import pysam
import shutil
import click
from pathlib import Path
import natsort

os.environ["CI"] = "1"
import tempfile

HVPATH = os.path.join(
    os.path.dirname(os.path.abspath(submodules.__file__)), "hv_rapidCNS2"
)


def run_methylartist(tempmgmtdir, plot_out):
    os.system(
        f"methylartist locus -i chr10:129466536-129467536 -b {os.path.join(tempmgmtdir, 'mgmt.bam')} -o {plot_out}  --motif CG --mods m > /dev/null 2>&1"
    )

def run_bedtools(bamfile, MGMT_BED, tempbamfile):
    """
    This function runs modkit on a bam file and extracts the methylation data.
    """
    try:
        os.system(f"bedtools intersect -a {bamfile} -b {MGMT_BED} > {tempbamfile}")
        pysam.index(tempbamfile)
    except Exception as e:
        print(e)


def run_modkit(tempmgmtdir, MGMTbamfile, threads):
    pysam.sort("-o", os.path.join(tempmgmtdir, "mgmt.bam"), MGMTbamfile)
    pysam.index(os.path.join(tempmgmtdir, "mgmt.bam"))
    os.system(
        f"modkit pileup -t {threads} --filter-threshold 0.73 --combine-mods {os.path.join(tempmgmtdir, 'mgmt.bam')} "
        f"{os.path.join(tempmgmtdir, 'mgmt.bed')} --suppress-progress  >/dev/null 2>&1 "
    )
    cmd = f"Rscript {HVPATH}/bin/mgmt_pred_v0.3.R --input={os.path.join(tempmgmtdir, 'mgmt.bed')} --out_dir={tempmgmtdir} --probes={HVPATH}/bin/mgmt_probes.Rdata --model={HVPATH}/bin/mgmt_137sites_mean_model.Rdata --sample=live_analysis"
    os.system(cmd)


class MGMT_Object(BaseAnalysis):
    def __init__(self, *args, **kwargs):
        self.MGMTbamfile = None
        self.counter = 0
        super().__init__(*args, **kwargs)

    def setup_ui(self):
        with ui.card().style("width: 100%"):
            ui.label("MGMT Methylation").style('color: #6E93D6; font-size: 150%; font-weight: 300').tailwind("drop-shadow", "font-bold")
            self.mgmtable = ui.row().classes("w-full")
            with self.mgmtable:
                ui.label("Table not yet available.")
            self.mgmtplot = ui.row().style("width: 100%")
            with self.mgmtplot:  # "size-full"):
                ui.label("Plot not yet available.")

    async def process_bam(self, bamfile, timestamp):
        MGMT_BED = f"{HVPATH}/bin/mgmt_hg38.bed"
        tempbamfile = tempfile.NamedTemporaryFile(dir=self.output, suffix=".bam")

        await run.cpu_bound(run_bedtools, bamfile, MGMT_BED, tempbamfile.name)

        if pysam.AlignmentFile(tempbamfile.name, "rb").count(until_eof=True) > 0:
            ui.notify("Running MGMT predictor - MGMT sites found.", type="positive", position="top")
            if not self.MGMTbamfile:
                #self.MGMTbamfile = tempfile.NamedTemporaryFile(dir=self.output, suffix=".bam")
                self.MGMTbamfile = os.path.join(self.output, "mgmt.bam")
                shutil.copy2(tempbamfile.name, self.MGMTbamfile)
            else:
                tempbamholder = tempfile.NamedTemporaryFile(dir=self.output, suffix=".bam")
                pysam.cat(
                    "-o", tempbamholder.name, self.MGMTbamfile, tempbamfile.name
                )
                shutil.copy2(tempbamholder.name, self.MGMTbamfile)
            tempmgmtdir = tempfile.TemporaryDirectory(dir=self.output)

            await run.cpu_bound(run_modkit, tempmgmtdir.name, self.MGMTbamfile, self.threads)
            ui.notify("MGMT predictor done.", type="positive", position="top")
            results = pd.read_csv(
                os.path.join(tempmgmtdir.name, "live_analysis_mgmt_status.csv")
            )
            self.counter += 1
            plot_out = os.path.join(self.output, f"{self.counter}_mgmt.png")
            await run.cpu_bound(run_methylartist, tempmgmtdir.name, plot_out)

            self.mgmtable.clear()
            with self.mgmtable:
                self.tabulate(results)
            results.to_csv(os.path.join(self.output, f"{self.counter}_mgmt.csv"))
            if os.path.exists(plot_out):
                self.mgmtplot.clear()
                with self.mgmtplot.classes("w-full"):
                    ui.image(plot_out).props("fit=scale-up")
            tempmgmtdir.cleanup()
            ui.notify("MGMT predictor complete.", type="positive", position="top")
        else:
            ui.notify("No new MGMT sites found.", type="warning", position="top")
        await asyncio.sleep(0.1)
        self.running = False

    def tabulate(self, results):
        ui.aggrid.from_pandas(
            results,
            theme="material",
            options={
                "defaultColDef": {
                    "sortable": True,
                    "resizable": True,
                },
                "columnDefs": [
                    {
                        "headerName": "Average",
                        "field": "average",
                        "filter": "agTextColumnFilter",
                        "floatingFilter": False,
                    },
                    {
                        "headerName": "Score",
                        "field": "pred",
                        "filter": "agNumberColumnFilter",
                        "floatingFilter": False,
                    },
                    {
                        "headerName": "Status",
                        "field": "status",
                        "filter": "agNumberColumnFilter",
                        "floatingFilter": False,
                    },
                ],
                "pagination": True,
                # "paginationAutoPageSize": True,
            },
        ).classes("w-full").style("height: 200px")
    def show_previous_data(self, watchfolder):
        for file in natsort.natsorted(os.listdir(watchfolder)):
            if file.endswith("_mgmt.csv"):
                results = pd.read_csv(os.path.join(watchfolder, file))
                plot_out = os.path.join(watchfolder, file.replace(".csv", ".png"))
                self.mgmtable.clear()
                with self.mgmtable:
                    self.tabulate(results)
                if os.path.exists(plot_out):
                    self.mgmtplot.clear()
                    with self.mgmtplot.classes("w-full"):
                        ui.image(plot_out).props('fit=scale-up')


def test_me(port: int, threads: int, watchfolder: str, output:str, reload: bool = False, browse: bool = False):
    my_connection = None
    with theme.frame("MGMT Data", my_connection):
        TestObject = MGMT_Object(threads, output, progress=True)
    if not browse:
        path = watchfolder
        searchdirectory = os.fsencode(path)
        for root, d_names, f_names in os.walk(searchdirectory):
            directory = os.fsdecode(root)
            for f in f_names:
                filename = os.fsdecode(f)
                if filename.endswith(".bam"):
                    TestObject.add_bam(os.path.join(directory, filename))
    else:
        TestObject.progress_trackers.visible=False
        TestObject.show_previous_data(output)
    ui.run(port=port,reload=reload)

@click.command()
@click.option(
    "--port",
    default=12345,
    help="Port for GUI",
)
@click.option(
    "--threads",
    default=4,
    help="Number of threads available."
)
@click.argument(
    "watchfolder",
    type=click.Path(
        exists=True, file_okay=False, dir_okay=True, resolve_path=True, path_type=Path
    ),
    required=False,
)
@click.argument(
    "output",
    type=click.Path(
        exists=True, file_okay=False, dir_okay=True, resolve_path=True, path_type=Path
    ),
    required=False,
)
@click.option(
    "--browse",
    is_flag=True,
    show_default=True,
    default=False,
    help="Browse Historic Data.",
)
def main(port, threads, watchfolder, output, browse):
    """
    Helper function to run the app.
    :param port: The port to serve the app on.
    :param reload: Should we reload the app on changes.
    :return:
    """
    if browse:
        # Handle the case when --browse is set
        click.echo("Browse mode is enabled. Only the output folder is required.")
        test_me(
            port=port,
            reload=False,
            threads=threads,
            #simtime=simtime,
            watchfolder=None,
            output=output,
            #sequencing_summary=sequencing_summary,
            #showerrors=showerrors,
            browse=browse,
            #exclude=exclude,
        )
        # Your logic for browse mode
    else:
        # Handle the case when --browse is not set
        click.echo(f"Watchfolder: {watchfolder}, Output: {output}")
        if watchfolder is None or output is None:
            click.echo("Watchfolder and output are required when --browse is not set.")
            sys.exit(1)
        test_me(
            port=port,
            reload=False,
            threads=threads,
            #simtime=simtime,
            watchfolder=watchfolder,
            output=output,
            #sequencing_summary=sequencing_summary,
            #showerrors=showerrors,
            browse=browse,
            #exclude=exclude,
        )
    test_me(port, threads, watchfolder, browse)

if __name__ in {"__main__", "__mp_main__"}:
    print("GUI launched by auto-reload function.")
    main()


