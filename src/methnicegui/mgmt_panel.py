from nicegui import Tailwind, ui, app


class MGMT_Panel:
    def setup_ui(self, mgmt):
        with ui.card().style("width: 100%"):
            ui.label("MGMT Methylation").tailwind("drop-shadow", "font-bold")
            # self.mgmtcontent=ui.column()
            # with self.mgmtcontent:
            self.mgmtplot = ui.row()
            with self.mgmtplot.classes("w-full"):
                ui.label("Plot not yet available.")
                # ui.image("/tmp/run2/targetsbams/25_sorted.png").props("fit=scale-down")
            self.mgmtable = ui.row()
            with self.mgmtable:
                ui.label("Table not yet available.")
