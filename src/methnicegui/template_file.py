# Python imports.
from __future__ import annotations
from nicegui import Tailwind, ui, app
from nicegui.events import ValueChangeEventArguments
from methnicegui import images

import threading

import time

from random import random

import theme

from minknow_api.manager import Manager
import minknow_api.manager_pb2 as manager_pb2
from minknow_api.statistics_pb2 import (
    ReadLengthHistogramSplit,
    ReadEndReason,
    DataSelection,
)

import time
import numpy as np
import math
import os




def index_page() -> None:
    #from check_connection import ConnectionDialog
    #initial_ip = "127.0.0.1"
    #my_connection = ConnectionDialog(initial_ip)
    my_connection = None
    with theme.frame('MethClass Interactive',my_connection):
        #my_connection.connect_to_minknow()
        ui.label(f"Hello")
        #my_object = MinknowHistograms(my_connection.positions[0])




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


def main(): # , threads, simtime, watchfolder, output, sequencing_summary):
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