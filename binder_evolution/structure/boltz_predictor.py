"""
Boltz-based structure predictor for protein–peptide complexes.

Boltz (https://github.com/jwohlwend/boltz) is a biomolecular structure
predictor trained on the PDB.  This module wraps its CLI so that
ProdigyBindingPredictor can request a fresh complex structure for every
peptide sequence variant produced by the evolutionary optimizer.

Install:
    pip install boltz

Typical usage::

    from binder_evolution.structure.boltz_predictor import BoltzStructurePredictor
    from binder_evolution.binding.prodigy_predictor import ProdigyBindingPredictor

    boltz = BoltzStructurePredictor(
        receptor_sequence="MPRFM...",   # full receptor AA sequence
        output_dir=Path("boltz_out"),
        use_msa_server=True,            # auto-generate MSAs via mmseqs2
    )

    predictor = ProdigyBindingPredictor(
        structure_fn=boltz,             # callable interface
        receptor_chains="A",
        peptide_chains="B",
    )
"""

from __future__ import annotations

import hashlib
import logging
import subprocess
import tempfile
import textwrap
import time
from pathlib import Path
from typing import Optional

log = logging.getLogger(__name__)


class BoltzStructurePredictor:
    """
    Callable structure predictor that runs Boltz to produce a protein–peptide
    complex for any given peptide sequence.

    This class satisfies the ``structure_fn(sequence, target_pdb) -> Path``
    interface expected by ``ProdigyBindingPredictor``.  The *target_pdb*
    argument is accepted but ignored — the receptor is provided as a plain
    amino-acid sequence in ``receptor_sequence``.

    Predictions are cached on disk under ``output_dir``: if a CIF file
    already exists for a sequence hash it is returned immediately without
    re-running Boltz.

    Args:
        receptor_sequence:
            Single-letter amino-acid sequence of the receptor protein.
        output_dir:
            Root directory for Boltz outputs.  Created if it does not exist.
        receptor_chain:
            Chain ID assigned to the receptor in the Boltz YAML (default "A").
        peptide_chain:
            Chain ID assigned to the peptide in the Boltz YAML (default "B").
        recycling_steps:
            Number of recycling steps passed to Boltz (default 3).
        diffusion_samples:
            Number of diffusion samples; only the top-ranked model is used
            (default 1).
        sampling_steps:
            Number of diffusion sampling steps (default 200).
        use_msa_server:
            If True, pass ``--use_msa_server`` to Boltz so that MSAs are
            generated automatically via mmseqs2.  Requires internet access
            (default False — no MSA, faster but potentially less accurate).
        use_potentials:
            If True, pass ``--use_potentials`` for inference-time potentials
            (default False).
        boltz_cmd:
            Name or path of the Boltz executable (default ``"boltz"``).
        timeout:
            Maximum seconds to wait for a single Boltz call before raising
            ``TimeoutError`` (default 1800 — 30 minutes).
    """

    def __init__(
        self,
        receptor_sequence: str,
        output_dir: Path,
        receptor_chain: str = "A",
        peptide_chain: str = "B",
        recycling_steps: int = 3,
        diffusion_samples: int = 1,
        sampling_steps: int = 200,
        use_msa_server: bool = False,
        use_potentials: bool = False,
        boltz_cmd: str = "boltz",
        timeout: int = 1800,
    ) -> None:
        self.receptor_sequence = receptor_sequence.strip().upper()
        self.output_dir = Path(output_dir)
        self.receptor_chain = receptor_chain
        self.peptide_chain = peptide_chain
        self.recycling_steps = recycling_steps
        self.diffusion_samples = diffusion_samples
        self.sampling_steps = sampling_steps
        self.use_msa_server = use_msa_server
        self.use_potentials = use_potentials
        self.boltz_cmd = boltz_cmd
        self.timeout = timeout

        self.output_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Callable interface
    # ------------------------------------------------------------------

    def __call__(
        self, sequence: str, target_pdb: Optional[Path] = None
    ) -> Optional[Path]:
        """
        Generate (or retrieve from cache) the complex CIF for *sequence*.

        Args:
            sequence:   Peptide amino-acid sequence (single-letter code).
            target_pdb: Ignored; receptor is taken from ``receptor_sequence``.

        Returns:
            Path to the predicted complex ``.cif`` file, or ``None`` on failure.
        """
        sequence = sequence.strip().upper()
        job_name = self._job_name(sequence)
        cif_path = self._expected_cif(job_name)

        if cif_path.exists():
            log.debug("Cache hit for %s → %s", sequence, cif_path)
            return cif_path

        yaml_path = self._write_yaml(sequence, job_name)
        success = self._run_boltz(yaml_path, job_name)
        if not success:
            return None

        if not cif_path.exists():
            log.error("Boltz finished but expected CIF not found: %s", cif_path)
            return None

        return cif_path

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _job_name(self, sequence: str) -> str:
        """Stable, filesystem-safe identifier for *sequence*."""
        digest = hashlib.sha1(sequence.encode()).hexdigest()[:10]
        return f"complex_{digest}"

    def _expected_cif(self, job_name: str) -> Path:
        """
        Boltz writes the top-ranked model to:
            <output_dir>/predictions/<job_name>/<job_name>_model_0.cif
        """
        return self.output_dir / "predictions" / job_name / f"{job_name}_model_0.cif"

    def _write_yaml(self, sequence: str, job_name: str) -> Path:
        """Write the Boltz input YAML and return its path."""
        yaml_dir = self.output_dir / "inputs"
        yaml_dir.mkdir(parents=True, exist_ok=True)
        yaml_path = yaml_dir / f"{job_name}.yaml"

        content = textwrap.dedent(f"""\
            version: 1
            sequences:
              - protein:
                  id: {self.receptor_chain}
                  sequence: {self.receptor_sequence}
              - protein:
                  id: {self.peptide_chain}
                  sequence: {sequence}
        """)
        yaml_path.write_text(content)
        log.debug("Wrote Boltz YAML: %s", yaml_path)
        return yaml_path

    def _run_boltz(self, yaml_path: Path, job_name: str) -> bool:
        """
        Invoke Boltz via subprocess.

        Returns True on success, False on failure.
        """
        cmd = [
            self.boltz_cmd,
            "predict",
            str(yaml_path),
            "--out_dir", str(self.output_dir),
            "--recycling_steps", str(self.recycling_steps),
            "--diffusion_samples", str(self.diffusion_samples),
            "--sampling_steps", str(self.sampling_steps),
        ]
        if self.use_msa_server:
            cmd.append("--use_msa_server")
        if self.use_potentials:
            cmd.append("--use_potentials")

        log.info("Running Boltz for job %s …", job_name)
        log.debug("Command: %s", " ".join(cmd))

        t0 = time.monotonic()
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=self.timeout,
            )
        except subprocess.TimeoutExpired:
            log.error("Boltz timed out after %d s for job %s", self.timeout, job_name)
            return False
        except FileNotFoundError:
            raise RuntimeError(
                f"Boltz executable '{self.boltz_cmd}' not found. "
                "Install it with:  pip install boltz"
            )

        elapsed = time.monotonic() - t0
        if result.returncode != 0:
            log.error(
                "Boltz failed (exit %d, %.1f s) for job %s.\nstderr:\n%s",
                result.returncode,
                elapsed,
                job_name,
                result.stderr[-2000:],  # last 2 KB
            )
            return False

        log.info("Boltz finished in %.1f s for job %s", elapsed, job_name)
        return True
