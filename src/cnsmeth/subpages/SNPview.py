# Python imports.
from __future__ import annotations
from nicegui import ui
import pandas as pd
import sys
import numpy as np


class SNPview:
    def __init__(self):
        super().__init__()
        self.snptable = None
        self.pathogenictable = None

    def renderme(self):
        self.placeholder = ui.card().tight().classes("overflow-x-auto")
        with self.placeholder:
            ui.label(
                "Candidate SNPs will be displayed here. SNPs are called based on available data at that time."
            )

    def process_annotations(self, record: dict) -> dict:
        """
        This function takes a dictionary record from a vcf file and explodes the record into multiple records based on the contents of the INFO field.
        We expect some unit of 16 entries in the INFO field. Where there are multiples of 16 entries, we split them into a new record entry for that specific mutation.
        """
        # print(record["INFO"])
        if "INFO" not in record.keys():
            return {}, {}
        annotations = record["INFO"]
        # This dictionary holds the information for a single record
        rec_dict = {}

        # This dictionary holds one or more records derived from the annotation field.
        ann_dict = {}

        for ann in annotations.split(";"):
            if "=" in ann:
                mykey = ann.split("=")[0]
                myvalue = ann.split("=")[1]
                if "|" in myvalue:
                    if mykey == "ANN":
                        if len(myvalue.split("|")) == 16:
                            try:
                                myvalues = myvalue.split("|")
                                count = 0
                                ann_dict[count] = dict()
                                ann_dict[count]["Allele"] = myvalues[0]
                                ann_dict[count]["Annotation"] = myvalues[1]
                                ann_dict[count]["Annotation_Impact"] = myvalues[2]
                                ann_dict[count]["Gene_Name"] = myvalues[3]
                                ann_dict[count]["Gene_ID"] = myvalues[4]
                                ann_dict[count]["Feature_Type"] = myvalues[5]
                                ann_dict[count]["Feature_ID"] = myvalues[6]
                                ann_dict[count]["Transcript_BioType"] = myvalues[7]
                                ann_dict[count]["Rank"] = myvalues[8]
                                ann_dict[count]["HGVS.c"] = myvalues[9]
                                ann_dict[count]["HGVS.p"] = myvalues[10]
                                ann_dict[count]["cDNA.pos / cDNA.length"] = myvalues[11]
                                ann_dict[count]["CDS.pos / CDS.length"] = myvalues[12]
                                ann_dict[count]["AA.pos / AA.length"] = myvalues[13]
                                ann_dict[count]["Distance"] = myvalues[14]
                                ann_dict[count]["ERRORS / WARNINGS / INFO"] = myvalues[
                                    15
                                ]
                            except Exception as e:
                                print(f"Error: {e}")
                                sys.exit()
                        elif len(myvalue.split("|")) > 16:
                            count = 0
                            for chunk in myvalue.split(","):
                                try:
                                    myvalues = chunk.split("|")
                                    ann_dict[count] = dict()
                                    ann_dict[count]["Allele"] = myvalues[0]
                                    ann_dict[count]["Annotation"] = myvalues[1]
                                    ann_dict[count]["Annotation_Impact"] = myvalues[2]
                                    ann_dict[count]["Gene_Name"] = myvalues[3]
                                    ann_dict[count]["Gene_ID"] = myvalues[4]
                                    ann_dict[count]["Feature_Type"] = myvalues[5]
                                    ann_dict[count]["Feature_ID"] = myvalues[6]
                                    ann_dict[count]["Transcript_BioType"] = myvalues[7]
                                    ann_dict[count]["Rank"] = myvalues[8]
                                    ann_dict[count]["HGVS.c"] = myvalues[9]
                                    ann_dict[count]["HGVS.p"] = myvalues[10]
                                    ann_dict[count]["cDNA.pos / cDNA.length"] = (
                                        myvalues[11]
                                    )
                                    ann_dict[count]["CDS.pos / CDS.length"] = myvalues[
                                        12
                                    ]
                                    ann_dict[count]["AA.pos / AA.length"] = myvalues[13]
                                    ann_dict[count]["Distance"] = myvalues[14]
                                    ann_dict[count]["ERRORS / WARNINGS / INFO"] = (
                                        myvalues[15]
                                    )
                                    # ann_dict[count]["unknown"] = "\n".join(myvalues[16:])
                                    count += 1
                                except Exception as e:
                                    print(f"Error: {e}")
                                    sys.exit()
                else:
                    # print(mykey, myvalue)
                    # print(len(myvalue.split("|")))
                    rec_dict[mykey] = myvalue
            # else:
            #    rec_dict[mykey] = myvalue

        return ann_dict, rec_dict

    def parse_vcf(self, vcf):
        with self.placeholder:
            header = "CHROM POS ID REF ALT QUAL FILTER INFO FORMAT GT".split()
            self.vcf = pd.read_csv(vcf, delimiter="\t", comment="#", names=header)
            if len(self.vcf) > 0:
                explodedvcf = []

                # print(self.vcf.to_dict('records')[:10])
                for record in self.vcf.to_dict("records"):
                    result, result2 = self.process_annotations(record)
                    if len(result) > 1:
                        # print("original record\n",record)
                        # print("Annotation Dictionary:\n",result)
                        # print("Record Dictionary:\n",result2)
                        for res in result:
                            # print(record)
                            dat = {**record, **result[res], **result2}
                            explodedvcf.append(dat)
                            # print(dat)
                            # sys.exit()

                self.vcf = pd.DataFrame.from_records(explodedvcf)
                if "INFO" in self.vcf.columns:
                    self.vcf = self.vcf.drop(columns=["INFO"]).drop_duplicates()
                else:
                    self.vcf = self.vcf.drop_duplicates()

                def set_unique_values(series):
                    set_series = set(series)
                    if len(set_series) > 0:
                        return "{}".format(", ".join(map(str, set_series)))
                    return None

                if len(self.vcf) > 0:
                    # print (self.vcf)
                    # Define shared columns
                    shared_columns = [
                        "CHROM",
                        "POS",
                        "ID",
                        "REF",
                        "ALT",
                        "QUAL",
                        "FILTER",
                        "FORMAT",
                        "GT",
                        "Allele",
                    ]

                    # Define columns to be aggregated
                    non_shared_columns = [
                        col for col in self.vcf.columns if col not in shared_columns
                    ]

                    self.vcf = self.vcf.replace({np.nan: None})
                    result = (
                        self.vcf.groupby(shared_columns)[non_shared_columns]
                        .agg(set_unique_values)
                        .reset_index()
                    )

                    self.vcf = result

                    self.placeholder.clear()
                    with self.placeholder:
                        self.snptable = (
                            ui.table.from_pandas(self.vcf, pagination=25)
                            .props("dense")
                            .style("height: 900px")
                            .style("font-size: 100%; font-weight: 300")
                        )
                        for col in self.snptable.columns:
                            col["sortable"] = True

                        def toggle(column: dict, visible: bool) -> None:
                            column["classes"] = "" if visible else "hidden"
                            column["headerClasses"] = "" if visible else "hidden"
                            self.snptable.update()

                        def set_pathogenic(value: bool) -> None:
                            self.snptable.filter = "pathogenic" if value else None

                        with self.snptable.add_slot("top-left"):

                            def toggle_fs() -> None:
                                self.snptable.toggle_fullscreen()
                                button.props(
                                    "icon=fullscreen_exit"
                                    if self.snptable.is_fullscreen
                                    else "icon=fullscreen"
                                )

                            button = ui.button(
                                "Toggle fullscreen",
                                icon="fullscreen",
                                on_click=toggle_fs,
                            ).props("flat")

                            with ui.button(icon="menu"):
                                with ui.menu(), ui.column().classes("gap-0 p-2"):
                                    for column in self.snptable.columns:
                                        ui.switch(
                                            column["label"],
                                            value=True,
                                            on_change=lambda e, column=column: toggle(
                                                column, e.value
                                            ),
                                        )

                            ui.switch(
                                "Show potentially pathogenic SNPs only",
                                value=False,
                                on_change=lambda e: set_pathogenic(e.value),
                            )

                        with self.snptable.add_slot("top-right"):
                            with ui.input(placeholder="Search").props(
                                "type=search"
                            ).bind_value(self.snptable, "filter").add_slot("append"):
                                ui.icon("search")

                        # with self.snptable.add_slot('top-left'):
                        #    def toggle() -> None:
                        #        self.snptable.toggle_fullscreen()
                        #        button.props('icon=fullscreen_exit' if self.snptable.is_fullscreen else 'icon=fullscreen')

                        #    button = ui.button('Toggle fullscreen', icon='fullscreen', on_click=toggle).props('flat')

                    # else:
                    #    self.snptable.update_rows(self.vcf.to_dict(orient='records'))


def index_page() -> None:
    # from check_connection import ConnectionDialog
    # initial_ip = "127.0.0.1"
    # my_connection = ConnectionDialog(initial_ip)
    # my_connection = None
    # with theme.frame("MethClass Interactive"):
    # my_connection.connect_to_minknow()
    mySNPview = SNPview()
    mySNPview.renderme()
    mySNPview.parse_vcf("test_intraop0049_b/clair3/snpsift_output.vcf")

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
