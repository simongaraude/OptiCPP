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
        diffusion_samples=5,            # predict 5 structures, use the best
        min_iptm=0.6,                   # discard poorly modelled complexes
    )

    predictor = ProdigyBindingPredictor(
        structure_fn=boltz,             # callable interface
        receptor_chains="A",
        peptide_chains="B",
    )
"""

from __future__ import annotations

import hashlib
import json
import logging
import subprocess
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

    When ``diffusion_samples > 1``, Boltz produces multiple candidate
    structures.  This class reads the per-model confidence JSON that Boltz
    writes alongside each CIF and returns the model with the highest ipTM
    (interface pTM) score — the standard quality metric for predicted
    protein–protein complexes.  Models whose ipTM falls below ``min_iptm``
    are rejected entirely (returns ``None``), signalling to the downstream
    PRODIGY predictor that the variant should receive a score of 0.

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
            Number of diffusion samples to generate per sequence.  The model
            with the highest ipTM is selected (default 1).
        sampling_steps:
            Number of diffusion sampling steps (default 200).
        min_iptm:
            Minimum acceptable ipTM for a predicted complex.  Structures below
            this threshold are treated as prediction failures and ``None`` is
            returned (default 0.5).  Set to 0.0 to disable the filter.
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
        min_iptm: float = 0.5,
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
        self.min_iptm = min_iptm
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
        Generate (or retrieve from cache) the best-ranked complex CIF for
        *sequence*.

        Args:
            sequence:   Peptide amino-acid sequence (single-letter code).
            target_pdb: Ignored; receptor is taken from ``receptor_sequence``.

        Returns:
            Path to the best-ranked ``.cif`` file whose ipTM ≥ ``min_iptm``,
            or ``None`` if Boltz fails or all models fall below the threshold.
        """
        sequence = sequence.strip().upper()
        job_name = self._job_name(sequence)

        # Only re-run Boltz if model_0 (the minimum output) is absent
        if not self._cif_path(job_name, 0).exists():
            yaml_path = self._write_yaml(sequence, job_name)
            if not self._run_boltz(yaml_path, job_name):
                return None

        return self._best_model(job_name)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _job_name(self, sequence: str) -> str:
        """Stable, filesystem-safe identifier for *sequence*."""
        digest = hashlib.sha256(sequence.encode()).hexdigest()[:12]
        return f"complex_{digest}"

    def _cif_path(self, job_name: str, model_idx: int) -> Path:
        """
        Boltz writes each model to:
            <output_dir>/predictions/<job_name>/<job_name>_model_<i>.cif
        """
        return (
            self.output_dir
            / "predictions"
            / job_name
            / f"{job_name}_model_{model_idx}.cif"
        )

    def _confidence_path(self, job_name: str, model_idx: int) -> Path:
        """
        Boltz writes per-model confidence to:
            <output_dir>/predictions/<job_name>/confidence_<job_name>_model_<i>.json
        """
        return (
            self.output_dir
            / "predictions"
            / job_name
            / f"confidence_{job_name}_model_{model_idx}.json"
        )

    def _read_iptm(self, job_name: str, model_idx: int) -> Optional[float]:
        """
        Return the ipTM score for *model_idx*, or None if the file is missing
        or malformed.
        """
        conf_path = self._confidence_path(job_name, model_idx)
        if not conf_path.exists():
            log.warning("Confidence file not found: %s", conf_path)
            return None
        try:
            data = json.loads(conf_path.read_text())
            return float(data["iptm"])
        except Exception:
            log.warning("Could not read ipTM from %s", conf_path)
            return None

    def _best_model(self, job_name: str) -> Optional[Path]:
        """
        Select the model with the highest ipTM from all generated samples.

        Returns the CIF path of the winner, or None if no model meets
        ``min_iptm``.
        """
        best_cif: Optional[Path] = None
        best_iptm: float = -1.0

        for i in range(self.diffusion_samples):
            cif = self._cif_path(job_name, i)
            if not cif.exists():
                continue  # Boltz may produce fewer models than requested

            iptm = self._read_iptm(job_name, i)
            if iptm is None:
                continue

            log.debug("Model %d: ipTM=%.3f (%s)", i, iptm, cif.name)
            if iptm > best_iptm:
                best_iptm = iptm
                best_cif = cif

        if best_cif is None:
            log.error("No valid Boltz models found for job %s", job_name)
            return None

        if best_iptm < self.min_iptm:
            log.warning(
                "Best ipTM %.3f < threshold %.3f for job %s — rejecting complex",
                best_iptm,
                self.min_iptm,
                job_name,
            )
            return None

        log.info(
            "Selected model %s (ipTM=%.3f) for job %s",
            best_cif.name,
            best_iptm,
            job_name,
        )
        return best_cif

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
                result.stderr[-2000:],
            )
            return False

        log.info("Boltz finished in %.1f s for job %s", elapsed, job_name)
        return True
