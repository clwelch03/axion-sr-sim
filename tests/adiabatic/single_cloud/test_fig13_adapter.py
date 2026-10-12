"""Check the relocated Python-to-Julia boundary without running a mode solve."""
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from adiabatic.single_cloud import fig13_points


class FigureAdapterChecks(unittest.TestCase):
    def test_batch_uses_existing_model_driver_and_project(self):
        references = {"n=0": {0.42: .4088 + 1e-7j},
                      "n=1": {0.42: .415 + 5e-8j}}

        def run_julia(command, check):
            self.assertTrue(check)
            self.assertEqual(command[0], "test-julia")
            project = Path(command[1].removeprefix("--project="))
            self.assertTrue((project / "Project.toml").is_file())
            self.assertEqual(Path(command[2]), project /
                             "src/adiabatic/single_cloud/drivers/fig13_points.jl")
            self.assertTrue(Path(command[2]).is_file())
            seeds = np.loadtxt(command[3], delimiter=",")
            np.testing.assert_array_equal(seeds[:, 0], [0, 1])
            np.testing.assert_allclose(seeds[:, 1], .42)
            np.savetxt(command[4],
                       [[0, .42, .4088, 1e-7, 1e-12],
                        [1, .42, .415, 5e-8, 1e-12]],
                       delimiter=",", comments="",
                       header="radial_overtone,alpha,omega_real,omega_imag,residual_norm")

        with tempfile.TemporaryDirectory() as directory:
            with patch.object(fig13_points.subprocess, "run", side_effect=run_julia):
                points = fig13_points.growth_rate_points(
                    [.42], references, .99, Path(directory) / "points.csv",
                    julia="test-julia")
        np.testing.assert_array_equal(points["n=0"], [.4088 + 1e-7j])
        np.testing.assert_array_equal(points["n=1"], [.415 + 5e-8j])

    def test_unconverged_output_is_rejected(self):
        def run_julia(command, check):
            Path(command[4]).write_text(
                "radial_overtone,alpha,omega_real,omega_imag,residual_norm\n"
                "0,.42,.4088,1e-7,1e-4\n1,.42,.415,5e-8,1e-12\n")

        references = {"n=0": {.42: .4088 + 1e-7j}, "n=1": {.42: .415 + 5e-8j}}
        with tempfile.TemporaryDirectory() as directory:
            with patch.object(fig13_points.subprocess, "run", side_effect=run_julia):
                with self.assertRaisesRegex(ValueError, "unconverged"):
                    fig13_points.growth_rate_points(
                        [.42], references, .99, Path(directory) / "points.csv")


if __name__ == "__main__":
    unittest.main()
