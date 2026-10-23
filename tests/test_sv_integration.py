import unittest
from pathlib import Path
import tomllib
from robin import cli, workflow_ray
from robin.state_tracker import DEFAULT_STAGES
from robin.gui.display_config import is_section_visible

class SVIntegrationTests(unittest.TestCase):
    def test_native_workflow_registration(self):
        self.assertIn('sv',cli.VALID_JOB_TYPES)
        self.assertEqual(cli.QUEUE_MAPPING['sv'],'sv')
        self.assertEqual(cli._validate_workflow_steps(['sv','sv:sv']),['sv','sv:sv'])
        self.assertEqual(workflow_ray.QUEUE_TO_TYPES['sv'],{'sv'})
        self.assertIn('sv',workflow_ray.TRIGGERS['preprocessing'])
        self.assertIn('sv',workflow_ray.NEEDS_WORK_DIR)
        self.assertIsNotNone(workflow_ray._sv_handler)
    def test_viewer_and_checkpoint_stage(self):
        self.assertTrue(is_section_visible('sv',workflow_steps=['sv']))
        self.assertEqual(len({s.bit for s in DEFAULT_STAGES.values()}),len(DEFAULT_STAGES))
    def test_entrypoint_and_repeat_resource(self):
        root=Path(__file__).resolve().parents[1]
        config=tomllib.loads((root/'pyproject.toml').read_text())
        self.assertEqual(config['project']['scripts']['robin-sv-replay'],'robin.sv_replay.replay:main')
        import hashlib
        resource=root/'src/robin/resources/sv_grch38_tandem_repeats.bed'
        self.assertEqual(hashlib.sha256(resource.read_bytes()).hexdigest(),'729cbb38bff84534070033db1ae004b6368f6856eef12dd7ce93e8ee55f37407')
