from nicegui import Tailwind, ui, app
from methnicegui import theme
from cnv_from_bam import iterate_bam_file
import gff3_parser
import tqdm

import os

os.environ["CI"] = "1"

import natsort
import numpy as np
import pandas as pd
from dna_features_viewer import BiopythonTranslator, GraphicFeature, GraphicRecord

from matplotlib import pyplot as plt


from methnicegui import resources


STRAND={"+":1, "-":-1}
import queue


class Fusion_Panel:
    def __init__(self):
        self.gene_bed = os.path.join(
                os.path.dirname(os.path.abspath(resources.__file__)), "unique_genes.bed"
            )
        self.gene_gff3 = os.path.join(
            os.path.dirname(os.path.abspath(resources.__file__)), "gencode.v45.annotation.gff3"
        )
        self.gene_table = gff3_parser.parse_gff3(self.gene_gff3, verbose = False, parse_attributes = True)
        #print(self.gene_table[self.gene_table['gene_name'].eq('BRAF')])
        #print(self.gene_table[self.gene_table['gene_name'].eq('KIAA1549')])
        #print(self.gene_table[self.gene_table['gene_name'].eq('BRAF')]['Type'].unique())
        #for index, row in self.gene_table[self.gene_table['gene_name'].eq('BRAF')].iterrows():
        #    if row['Type'] == 'gene': # and row['transcript_type'] == 'protein_coding':
        #        print(row['Start'])
        #        print(row['End'])
        #        print(row['Strand'])
        #        print(row)



    def setup_ui(self, fusion):
        with ui.card().style("width: 100%"):
            ui.label("Fusion Candidates").tailwind("drop-shadow", "font-bold")
            # self.mgmtcontent=ui.column()
            # with self.mgmtcontent:
            self.fusionplot = ui.row()
            with self.fusionplot.classes("w-full"):
                ui.label("Plot not yet available.")
                # ui.image("/tmp/run2/targetsbams/25_sorted.png").props("fit=scale-down")
            self.fusiontable = ui.row()
            with self.fusiontable:
                ui.label("Table not yet available.")

    def create_fusion_plot(self, title,reads):
        #ui.label(title).tailwind("drop-shadow", "font-bold")
        #print(reads)
        with ui.card().classes('no-shadow border-[1px]'):
            with ui.pyplot(figsize=(15, 2)):
                ax = plt.gca()
                plt.title(title)
                features = []
                tracker = []
                first_index = 0
                sequence_length = 0
                for index, row in self.gene_table[self.gene_table['gene_name'].eq(title)].iterrows():
                    if [row['Start'], row['End']] in tracker:
                        continue
                    else:
                        # tracker.append([row['Start'], row['End']])
                        if row['Type'] == 'gene':
                            features.append(
                                GraphicFeature(start=int(row['Start']), end=int(row['End']), strand=STRAND[row['Strand']],
                                               thickness=2, color="#ffd700",
                                               #label=row['gene_name']
                                               )
                            )
                            first_index = int(row['Start']) - 1000
                            sequence_length = int(row['End']) - int(row['Start']) + 2000
                        if row['Type'] == 'CDS':  # and row['transcript_type'] == 'protein_coding':
                            features.append(
                                GraphicFeature(start=int(row['Start']), end=int(row['End']), strand=STRAND[row['Strand']],
                                               color="#ffcccc", ))  # label=row['exon_id']))
                for index,row in reads.iterrows():
                    features.append(
                        GraphicFeature(start=int(row[5]), end=int(row[6]), strand=STRAND[row[9]],
                                       color="#ccffcc"))

                record = GraphicRecord(sequence_length=sequence_length, first_index=first_index, features=features)
                record.plot(ax=ax)

    def parse_bams(self, bampath):
        bamfiles = natsort.natsorted(os.listdir(bampath))
        # Check if all bamfiles have already been subset - if not subset them with the bed file and only keep reads with supplementary alignments
        for file in bamfiles:
            if file.endswith(".bam"):
                if file[0].isdigit():
                    #print(file)
                    if not os.path.exists(os.path.join(bampath, f"subset_{file}")):
                        #print("Subsetting")
                        subset_file = f"subset_{file}"
                        os.system(
                         f"samtools view -L {self.gene_bed} -d SA -o {os.path.join(bampath, subset_file)} {os.path.join(bampath, file)}"
                        )
                        #os.system(f"samtools index {os.path.join(bampath, file)}")

        # Now we merge the newly formed bamfiles:
        os.system(f"samtools cat -o {os.path.join(bampath, 'merged.bam')} {os.path.join(bampath ,'subset_*')}")
        os.system(
            f"bedtools intersect -a {self.gene_bed} -b {os.path.join(bampath, 'merged.bam')} -wa -wb > {os.path.join(bampath ,'mappings.txt')}"
        )
        #print(f"File being read is {os.path.join(bampath ,'mappings.txt')}")
        self.fusion_candidates = pd.read_csv(os.path.join(bampath ,'mappings.txt'), sep="\t", header=None)

        self.fusion_candidates = self.fusion_candidates[self.fusion_candidates[8].gt(50)]
        uniques = self.fusion_candidates[7].duplicated(keep=False)
        doubles = self.fusion_candidates[uniques]
        counts = doubles.groupby(7)[3].transform('nunique')
        self.fusiontable.clear()
        #print(doubles[counts > 1].sort_values(by=7))
        #print(doubles[counts > 1].sort_values(by=7).groupby([3,5]).aggregate("count"))
        #print(doubles[counts > 1].sort_values(by=7).groupby([3,6]).aggregate("count"))
        #result = doubles[counts > 1].sort_values(by=7).groupby(3).apply(self.count_shared_values).reset_index(name='shared_count')
        #print(result)
        #print(type(result))
        with self.fusiontable:
            #ui.table(result).classes("w-full")
            #ui.table.from_pandas(result)
            ui.table.from_pandas(doubles[counts > 1].sort_values(by=7).rename(columns={0:"chromBED",1:"BS",2:"BE",3:"Gene",4:"chrom",5:"mS",6:"mE",7:"readID",8:"mapQ",9:"strand"})).classes("w-full")
        self.fusionplot.clear()
        with self.fusionplot.classes("w-full"):
            for gene in doubles[counts > 1].sort_values(by=7)[3].unique():
                if len(doubles[counts > 1].sort_values(by=7)[doubles[counts > 1].sort_values(by=7)[3].eq(gene)]) > 1:
                    self.create_fusion_plot(gene, doubles[counts > 1].sort_values(by=7)[doubles[counts > 1].sort_values(by=7)[3].eq(gene)])
            #self.create_fusion_plot("Fusion Plot")
            #self.update_fusion_plot(self.fusion_candidates)
            pass


    def count_shared_values(self, group):
        return group[7].nunique()

    def fusion_predictor(self, bamfile):
        #bedtools
        #intersect - a
        #~ / GIT / niceGUI / cnsmeth / src / methnicegui / resources / unique_genes.bed - b
        #sorted.bam - wa - wb > mappings.txt
        os.system(
            f"bedtools intersect -a  sort --write-index -@{self.threads} -o {self.sortedbamfile} {tempbam.name} "
            f">/dev/null 2>&1"
        )


def index_page() -> None:
    # from check_connection import ConnectionDialog
    # initial_ip = "127.0.0.1"
    # my_connection = ConnectionDialog(initial_ip)
    my_connection = None
    with theme.frame("MethClass Interactive", my_connection):
        # my_connection.connect_to_minknow()
        fusion = Fusion_Panel()
        fusion.setup_ui(None)
        fusion.parse_bams("/Users/mattloose/datasets/nicegui_ds1305/donebams/")

        #CNV_PLOT.create_cnv_scatter("CNV Scatter")
        #CNV_PLOT.cnv_plotting("/Users/mattloose/datasets/ds1305_sort.hg38.h2m.bam")
        # my_object = MinknowHistograms(my_connection.positions[0])


def run_class(port: int, reload: bool):
    """
    Helper function to run the app.
    :param port: The port to serve the app on.
    :param reload: Should we reload the app on changes.
    :return:
    """
    index_page()
    ui.run(
        port=port, reload=reload, title="MethClass NiceGUI"
    )  # , native=True, fullscreen=False, window_size=(1200, 1000))


def main():  # , threads, simtime, watchfolder, output, sequencing_summary):
    from check_connection import ConnectionDialog

    """
    Entrypoint for when GUI is launched directly.
    :return: None
    """
    run_class(port=12398, reload=False)


# Entrypoint for when GUI is launched by the CLI.
# e.g.: python my_app/my_cli.py
if __name__ in {"__main__", "__mp_main__"}:
    """
    Entrypoint for when GUI is launched by the CLI
    :return: None
    """
    if __name__ == "__mp_main__":
        print("GUI launched by auto-reload")

    run_class(port=12398, reload=True)
