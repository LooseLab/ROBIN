"""
File included for development purposes.

Needs to be refactored.

"""

# Python imports.
from __future__ import annotations
from nicegui import ui, app
from pathlib import Path

import theme
import numpy as np
import os
import natsort


def add_child(data_structure, strand_id, new_child):
    """
    Adds a new child to the children list of the specified strand (+ or -) in the given data structure.

    Args:
        data_structure (list): The list containing dictionaries for each strand.
        strand_id (str): The strand identifier, either '+' or '-'.
        new_child (dict): The new child item to be added to the children list.

    Returns:
        bool: True if the item was added successfully, False if the strand_id was not found.
    """
    for item in data_structure:
        if item["id"] == strand_id:
            item["children"].append(new_child)
            return True
    return False


@ui.page("/", response_timeout=30)
def index_page() -> None:
    # from check_connection import ConnectionDialog
    # initial_ip = "127.0.0.1"
    # my_connection = ConnectionDialog(initial_ip)
    # my_connection = None
    output = "/Users/mattloose/GIT/niceGUI/cnsmeth/robin_output_test/_ds1305_Intraop0051_c_ULK_2_highFRA"
    CNVResults = np.load(
        os.path.join(output, "ruptures.npy"), allow_pickle="TRUE"
    ).item()
    with theme.frame(
        "<strong><font color='#000000'>R</font></strong>apid nanop<strong><font color='#000000'>O</font></strong>re <strong><font color='#000000'>B</font></strong>rain intraoperat<strong><font color='#000000'>I</font></strong>ve classificatio<strong><font color='#000000'>N</font></strong>",
        smalltitle="<strong><font color='#000000'>R.O.B.I.N</font></strong>",
    ):
        # my_connection.connect_to_minknow()
        with ui.card().classes("w-full"):

            ui.label("Target Information")

            trees = []
            for key in natsort.natsorted(CNVResults):
                tree_dict = {}
                # print(key)
                tree_dict["id"] = key
                tree_dict["description"] = f"Targets On Chromosome {key}"
                tree_dict["children"] = [
                    {"id": "+", "description": "Forward Strand", "children": []},
                    {"id": "-", "description": "Reverse Strand", "children": []},
                ]
                # print(CNVResults[key])
                for line in CNVResults[key]["bed_data"].split("\n"):
                    # print(line)
                    # print(tree_dict['children'])
                    # Split each line into its components
                    chrom, start, end, _, _, strand = line.split("\t")
                    if strand == "+":
                        add_child(
                            tree_dict["children"],
                            strand,
                            {"id": f"{chrom}:{start}-{end}"},
                        )
                    else:
                        add_child(
                            tree_dict["children"],
                            strand,
                            {"id": f"{chrom}:{start}-{end}"},
                        )

                # print(dir(tree))
                # print(tree.id)
                trees.append(tree_dict)
            with ui.card().classes("w-full border-[1px]"):
                tree = ui.tree(
                    [
                        {
                            "id": "chromosome",
                            "description": "Targets On Each Chromosome",
                            "children": trees,
                        }
                    ],
                    label_key="id",
                    on_select=lambda e: ui.notify(e.value),
                )

                tree.add_slot(
                    "default-header",
                    """
                    <span :props="props"><strong>{{ props.node.id }}</strong></span>
                """,
                )
                tree.add_slot(
                    "default-body",
                    """
                    <span :props="props">Description: "{{ props.node.description }}"</span>
                """,
                )
                # ui.label(f"{CNVResults}")

                with ui.row():
                    ui.button("+ all", on_click=tree.expand)
                    ui.button("- all", on_click=tree.collapse)
                    ui.button("open + strand", on_click=lambda: tree.expand(["+"]))
                    ui.button("close + strand", on_click=lambda: tree.collapse(["+"]))

        # my_object = MinknowHistograms(my_connection.positions[0])


def run_class(port: int, reload: bool):
    """
    Helper function to run the app.
    :param port: The port to serve the app on.
    :param reload: Should we reload the app on changes.
    :return:
    """
    ui.add_css(
        """
        .shadows-into light-regular {
            font-family: "Shadows Into Light", cursive;
            font-weight: 400;
            font-style: normal;
        }
    """
    )
    app.add_static_files("/fonts", str(Path(__file__).parent / "fonts"))
    # app.on_startup(mainpage.index_page)
    # index_page()
    ui.run(
        port=port,
        reload=reload,
        title="MethClass NiceGUI",
        storage_secret="slartibartfast",
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
