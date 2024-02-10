from nicegui import ui


class MGMT_Panel:
    def setup_ui(self, mgmt):
        with ui.card().style("width: 100%"):
            ui.label("MGMT Methylation").tailwind("drop-shadow", "font-bold")
            self.mgmtplot = ui.row()
            with self.mgmtplot.style("width: 100%"):  # "size-full"):
                ui.label("Plot not yet available.")
            self.mgmtable = ui.row()
            with self.mgmtable:
                ui.label("Table not yet available.")
