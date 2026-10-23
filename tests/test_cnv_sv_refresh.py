import tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
from nicegui import ui
from robin.gui.components import cnv
from robin.gui.plotting_preferences import PlottingPreferencesConfig

class CNVSVRefreshTests(unittest.TestCase):
    def test_partial_files_keep_cache_eligible_for_retry(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            for name in ['CNV.npy','CNV3.npy']:(root/name).write_bytes(b'incomplete')
            payload=cnv._cnv_load_binary_payload(root,cnv_dict_npy_changed=False,cnv_npy_changed=True,cnv2_npy_changed=False,cnv3_npy_changed=True,data_array_reload=False,xy_pkl_changed=False)
            self.assertNotIn('cnv',payload);self.assertNotIn('cnv3',payload)
            self.assertTrue(payload['retry_cnv']);self.assertTrue(payload['retry_cnv3'])

    def test_highlight_uses_dom_and_retains_scatter_data(self):
        sample=Path('/tmp/cnv-sv-refresh-test');key=str(sample)
        launcher=SimpleNamespace(_cnv_state={},batman_mode=False,plotting_preferences=PlottingPreferencesConfig())
        def timer(*args,**kwargs):return SimpleNamespace(cancel=lambda:None,deactivate=lambda:None)
        with patch.object(cnv,'client_timer',timer),patch.object(cnv,'register_theme_sync_callback',return_value=lambda:None),patch.object(cnv,'_cnv_echart_push_update') as push,patch.object(ui,'run_javascript') as js,patch('robin.gui.plotting_preferences.resolve_cnv_summary_normalized',return_value=False):
            cnv.add_cnv_section(launcher,sample)
            launcher._cnv_state[key].update(cnv={'chr7':[2.0]*100,'chr9':[2.0]*150},cnv3={'chr7':[0.0]*100,'chr9':[0.0]*150},cnv_dict={'bin_width':1000000},xy='XX')
            callback=ui.context.client._sv_cnv_highlighters[key]
            for chrom in ['chr7','chr9','chr7']:
                callback(dict(id=1,c1=chrom,p1=51000000,c2=chrom,p2=56000000))
            scripts='\n'.join(str(call.args[0]) for call in js.call_args_list)
            self.assertIn('robin-sv-band',scripts)
            self.assertIn("transition:'opacity 1500ms ease-out'",scripts)
            self.assertNotIn('chart.setOption',scripts)
            self.assertNotIn('requestAnimationFrame(frame)',scripts)
            charts={call.args[0].id:call.args[0] for call in push.call_args_list}
            self.assertEqual(len(charts),2)
            for chart in charts.values():
                self.assertTrue(any(s.get('type')=='scatter' and s.get('data') for s in chart.options['series']))
                self.assertFalse(any(s.get('id')=='robin-sv-highlight' for s in chart.options['series']))
