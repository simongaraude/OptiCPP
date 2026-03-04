"""
PRODIGY-based binding affinity predictor.

PRODIGY (PROtein binDIng enerGY) predicts binding affinity (ΔG, kcal/mol)
and dissociation constant (Kd, M) from a protein–peptide complex PDB structure,
using intermolecular contact counts and non-interface surface residue fractions.

Reference:
    Vangone & Bonvin (2015). eLife, 4:e07454.
    Xue et al. (2016). Bioinformatics, 32(12):i534–i542.

GitHub: https://github.com/haddocking/prodigy
Install: pip install prodigy-prot
"""

from __future__ import annotations

import math
import logging
import tempfile
from pathlib import Path
from typing import Callable, Optional

log = logging.getLogger(__name__)


class ProdigyBindingPredictor:
    """
    Binding affinity predictor that uses PRODIGY to compute Kd from a
    protein–peptide complex PDB structure.

    Because PRODIGY requires a 3D structure, complex generation is delegated
    to a user-supplied callable (``structure_fn``). A reference Kd is
    established from the seed complex on the first call and is used to
    normalise all variant scores to [0, 1].

    Score normalisation (log10 scale):
        score = 1.0 − (log10(Kd_variant) − log10(Kd_reference)) / score_window
        clamped to [0, 1].

    This means:
        - variant Kd == reference Kd  →  score 1.0
        - variant Kd 10× weaker       →  score 0.5  (score_window=2)
        - variant Kd 100× weaker      →  score 0.0  (score_window=2)
        - variant Kd stronger         →  score capped at 1.0

    Args:
        structure_fn:
            Callable(sequence: str, target_pdb: Optional[Path]) -> Path.
            Must return the path to a PDB file of the fully assembled
            peptide–target complex ready for PRODIGY.
        target_pdb:
            Optional path to the target protein PDB, forwarded to
            ``structure_fn`` as its second argument.
        peptide_chains:
            Chain ID(s) of the peptide in the complex PDB.
            Single string (e.g. ``"B"``) or list (e.g. ``["B"]``).
        receptor_chains:
            Chain ID(s) of the receptor in the complex PDB.
            Single string (e.g. ``"A"``) or list (e.g. ``["A"]``).
        temperature:
            Temperature in °C used by PRODIGY for Kd calculation (default 25.0).
        distance_cutoff:
            Contact distance cutoff in Å (default 5.5).
        acc_threshold:
            Minimum fractional SASA for a residue to be counted as
            surface-exposed (default 0.05).
        score_window:
            Number of log10 Kd units over which the score decreases from
            1.0 to 0.0 (default 2.0, i.e. two orders of magnitude).
        reference_kd:
            If provided, skip the automatic reference computation and use
            this Kd (M) value directly. Useful when the seed complex
            structure is not available up-front.
    """

    def __init__(
        self,
        structure_fn: Callable[[str, Optional[Path]], Path],
        target_pdb: Optional[Path] = None,
        peptide_chains: str | list[str] = "B",
        receptor_chains: str | list[str] = "A",
        temperature: float = 25.0,
        distance_cutoff: float = 5.5,
        acc_threshold: float = 0.05,
        score_window: float = 2.0,
        reference_kd: Optional[float] = None,
    ) -> None:
        self.structure_fn = structure_fn
        self.target_pdb = target_pdb
        self.peptide_chains = (
            [peptide_chains] if isinstance(peptide_chains, str) else list(peptide_chains)
        )
        self.receptor_chains = (
            [receptor_chains] if isinstance(receptor_chains, str) else list(receptor_chains)
        )
        self.temperature = temperature
        self.distance_cutoff = distance_cutoff
        self.acc_threshold = acc_threshold
        self.score_window = score_window
        self.reference_kd: Optional[float] = reference_kd

        # Cache (sequence → kd) to avoid re-predicting the same variant
        self._kd_cache: dict[str, float] = {}

    # ------------------------------------------------------------------
    # Public interface (matches EvolutionController expectations)
    # ------------------------------------------------------------------

    def predict_binding(self, sequence: str) -> float:
        """
        Return a binding score in [0, 1] for *sequence*.

        On the very first call the result is stored as the reference Kd
        (score == 1.0). Subsequent variants are scored relative to it.

        Args:
            sequence: Amino-acid sequence of the peptide (single-letter code).

        Returns:
            Binding score in [0, 1]; 1.0 means at least as good as the
            reference, 0.0 means ``score_window`` log-units weaker.

        Raises:
            RuntimeError: If PRODIGY is not installed or structure
                          generation fails.
        """
        kd = self._get_kd(sequence)
        if kd is None:
            log.warning("PRODIGY prediction failed for %s – returning score 0.0", sequence)
            return 0.0

        # First successful prediction establishes the reference
        if self.reference_kd is None:
            self.reference_kd = kd
            log.info("Reference Kd set to %.3e M (sequence: %s)", kd, sequence)

        return self._kd_to_score(kd)

    def predict_kd(self, sequence: str) -> Optional[float]:
        """
        Return the raw predicted Kd (M) for *sequence*, or None on failure.
        """
        return self._get_kd(sequence)

    def predict_dg(self, sequence: str) -> Optional[float]:
        """
        Return the predicted ΔG (kcal/mol) for *sequence*, or None on failure.
        """
        complex_pdb = self._build_complex(sequence)
        if complex_pdb is None:
            return None
        result = self._run_prodigy(complex_pdb)
        if result is None:
            return None
        _kd, dg = result
        return dg

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _get_kd(self, sequence: str) -> Optional[float]:
        """Return Kd (M) from cache or by running PRODIGY."""
        if sequence in self._kd_cache:
            return self._kd_cache[sequence]

        complex_pdb = self._build_complex(sequence)
        if complex_pdb is None:
            return None

        result = self._run_prodigy(complex_pdb)
        if result is None:
            return None

        kd, _dg = result
        self._kd_cache[sequence] = kd
        return kd

    def _build_complex(self, sequence: str) -> Optional[Path]:
        """Call structure_fn to produce a complex PDB for *sequence*."""
        try:
            pdb_path = self.structure_fn(sequence, self.target_pdb)
            if pdb_path is None or not Path(pdb_path).exists():
                log.error("structure_fn returned no valid PDB for sequence %s", sequence)
                return None
            return Path(pdb_path)
        except Exception:
            log.exception("structure_fn raised an exception for sequence %s", sequence)
            return None

    def _run_prodigy(self, complex_pdb: Path) -> Optional[tuple[float, float]]:
        """
        Run PRODIGY on *complex_pdb* and return (Kd in M, ΔG in kcal/mol).

        Returns None if PRODIGY is not installed or prediction fails.
        """
        try:
            from Bio.PDB import PDBParser  # type: ignore
            from prodigy_prot.modules.prodigy import Prodigy  # type: ignore
        except ImportError as exc:
            raise RuntimeError(
                "PRODIGY is not installed. Install it with:\n"
                "    pip install prodigy-prot\n"
                f"(Original error: {exc})"
            ) from exc

        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("complex", str(complex_pdb))
            model = structure[0]

            selection = self.receptor_chains + self.peptide_chains
            prodigy = Prodigy(model, selection=selection, temp=self.temperature)
            prodigy.predict(
                distance_cutoff=self.distance_cutoff,
                acc_threshold=self.acc_threshold,
            )

            kd: float = prodigy.kd_val   # Kd in M
            dg: float = prodigy.ba_val   # ΔG in kcal/mol

            log.debug(
                "PRODIGY: %s → ΔG=%.2f kcal/mol, Kd=%.3e M",
                complex_pdb.name,
                dg,
                kd,
            )
            return kd, dg

        except Exception:
            log.exception("PRODIGY failed on %s", complex_pdb)
            return None

    def _kd_to_score(self, kd: float) -> float:
        """
        Convert a Kd (M) to a normalised binding score in [0, 1].

        Uses a log10-linear scale relative to ``self.reference_kd``.
        A variant ``score_window`` log-units weaker than the reference
        receives score 0.0; equal or stronger receives 1.0.
        """
        if self.reference_kd is None:
            raise RuntimeError("reference_kd not set – call predict_binding with the seed first.")

        log_kd = math.log10(kd)
        log_ref = math.log10(self.reference_kd)
        score = 1.0 - (log_kd - log_ref) / self.score_window
        return max(0.0, min(1.0, score))


# ---------------------------------------------------------------------------
# Convenience: PRODIGY on a pre-built complex (no structure prediction needed)
# ---------------------------------------------------------------------------

class StaticProdigyBindingPredictor(ProdigyBindingPredictor):
    """
    PRODIGY predictor for cases where a single fixed complex PDB is available.

    This is useful for validating the parent peptide's Kd without running
    structure prediction. All sequence variants receive the same fixed score
    and should be used together with the ``StructuralFilter`` (which generates
    per-variant structures) for a complete evaluation.

    Args:
        complex_pdb: Path to the protein–peptide complex PDB.
        peptide_chains: Chain ID(s) of the peptide.
        receptor_chains: Chain ID(s) of the receptor.
        temperature: Temperature in °C (default 25.0).
    """

    def __init__(
        self,
        complex_pdb: Path,
        peptide_chains: str | list[str] = "B",
        receptor_chains: str | list[str] = "A",
        temperature: float = 25.0,
        **kwargs,
    ) -> None:
        def _fixed_structure(sequence: str, target_pdb: Optional[Path]) -> Path:
            return complex_pdb

        super().__init__(
            structure_fn=_fixed_structure,
            peptide_chains=peptide_chains,
            receptor_chains=receptor_chains,
            temperature=temperature,
            **kwargs,
        )
        self._fixed_complex_pdb = complex_pdb

    def evaluate_reference(self) -> tuple[float, float]:
        """
        Run PRODIGY on the fixed complex and return (Kd in M, ΔG in kcal/mol).

        Sets ``self.reference_kd`` as a side-effect.
        """
        result = self._run_prodigy(self._fixed_complex_pdb)
        if result is None:
            raise RuntimeError(f"PRODIGY failed on {self._fixed_complex_pdb}")
        kd, dg = result
        if self.reference_kd is None:
            self.reference_kd = kd
        return kd, dg
