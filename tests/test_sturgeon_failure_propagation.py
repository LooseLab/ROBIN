"""Prediction failure must not be reported as successful analysis."""
import unittest
from unittest.mock import patch

from robin.analysis.sturgeon_analysis import sturgeon_bam_background_work


class SturgeonFailureTests(unittest.TestCase):
    def test_input_failure_propagates_to_workflow(self):
        with patch('robin.analysis.sturgeon_analysis.load_modkit_data',
                   side_effect=RuntimeError('prediction input failed')):
            with self.assertRaisesRegex(RuntimeError, 'prediction input failed'):
                sturgeon_bam_background_work('missing.parquet','/tmp','missing.bed',0)


if __name__=='__main__':unittest.main()
