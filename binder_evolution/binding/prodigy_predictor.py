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

import json
import math
import logging
from pathlib import Path
from typing import Callable, Optional

log = logging.getLogger(__name__)

# Type alias for a cached PRODIGY result
_ProdigyResult = tuple[float, float]  # (kd_M, dg_kcal_per_mol)


class ProdigyBindingPredictor:
    """
    Binding affinity predictor that uses PRODIGY to compute Kd from a
    protein–peptide complex PDB/CIF structure.

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
            Must return the path to a PDB or CIF file of the fully assembled
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

        # In-memory cache: sequence → (kd, dg).  Both values are stored so
        # predict_kd and predict_dg never run PRODIGY twice for the same sequence.
        self._result_cache: dict[str, _ProdigyResult] = {}

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
            RuntimeError: If PRODIGY is not installed.
        """
        result = self._get_result(sequence)
        if result is None:
            log.warning("PRODIGY prediction failed for %s – returning score 0.0", sequence)
            return 0.0

        kd, _ = result
        if self.reference_kd is None:
            self.reference_kd = kd
            log.info("Reference Kd set to %.3e M (sequence: %s)", kd, sequence)

        return self._kd_to_score(kd)

    def predict_kd(self, sequence: str) -> Optional[float]:
        """Return the raw predicted Kd (M) for *sequence*, or None on failure."""
        result = self._get_result(sequence)
        return result[0] if result is not None else None

    def predict_dg(self, sequence: str) -> Optional[float]:
        """Return the predicted ΔG (kcal/mol) for *sequence*, or None on failure."""
        result = self._get_result(sequence)
        return result[1] if result is not None else None

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _get_result(self, sequence: str) -> Optional[_ProdigyResult]:
        """
        Return (kd, dg) from in-memory cache, disk cache, or PRODIGY.

        Lookup order:
          1. In-memory cache  (fastest)
          2. Disk cache (.prodigy.json next to the CIF)  (survives restarts)
          3. Run PRODIGY and populate both caches
        """
        if sequence in self._result_cache:
            return self._result_cache[sequence]

        complex_path = self._build_complex(sequence)
        if complex_path is None:
            return None

        # Try disk cache before running PRODIGY
        cached = self._load_disk_cache(complex_path)
        if cached is not None:
            self._result_cache[sequence] = cached
            return cached

        result = self._run_prodigy(complex_path)
        if result is None:
            return None

        self._result_cache[sequence] = result
        self._save_disk_cache(complex_path, result)
        return result

    def _build_complex(self, sequence: str) -> Optional[Path]:
        """Call structure_fn to produce a complex PDB/CIF for *sequence*."""
        try:
            path = self.structure_fn(sequence, self.target_pdb)
            if path is None or not Path(path).exists():
                log.error("structure_fn returned no valid structure file for sequence %s", sequence)
                return None
            return Path(path)
        except Exception:
            log.exception("structure_fn raised an exception for sequence %s", sequence)
            return None

    def _run_prodigy(self, complex_path: Path) -> Optional[_ProdigyResult]:
        """
        Run PRODIGY on *complex_path* and return (Kd in M, ΔG in kcal/mol).

        Accepts both PDB (``.pdb``) and mmCIF (``.cif``) files.  Boltz
        produces ``.cif`` output, so mmCIF support is required when using
        ``BoltzStructurePredictor``.

        Returns None if PRODIGY is not installed or prediction fails.
        """
        try:
            from Bio.PDB import MMCIFParser, PDBParser  # type: ignore
            from prodigy_prot.modules.prodigy import Prodigy  # type: ignore
        except ImportError as exc:
            raise RuntimeError(
                "PRODIGY is not installed. Install it with:\n"
                "    pip install prodigy-prot\n"
                f"(Original error: {exc})"
            ) from exc

        try:
            suffix = complex_path.suffix.lower()
            if suffix in {".cif", ".mmcif"}:
                parser = MMCIFParser(QUIET=True)
            else:
                parser = PDBParser(QUIET=True)

            structure = parser.get_structure("complex", str(complex_path))
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
                complex_path.name,
                dg,
                kd,
            )
            return kd, dg

        except Exception:
            log.exception("PRODIGY failed on %s", complex_path)
            return None

    # ------------------------------------------------------------------
    # Disk cache helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _disk_cache_path(complex_path: Path) -> Path:
        return complex_path.with_suffix(complex_path.suffix + ".prodigy.json")

    def _load_disk_cache(self, complex_path: Path) -> Optional[_ProdigyResult]:
        cache_file = self._disk_cache_path(complex_path)
        if not cache_file.exists():
            return None
        try:
            data = json.loads(cache_file.read_text())
            result: _ProdigyResult = (float(data["kd"]), float(data["dg"]))
            log.debug("Disk cache hit: %s", cache_file)
            return result
        except Exception:
            log.warning("Corrupt PRODIGY disk cache, ignoring: %s", cache_file)
            return None

    def _save_disk_cache(self, complex_path: Path, result: _ProdigyResult) -> None:
        cache_file = self._disk_cache_path(complex_path)
        try:
            cache_file.write_text(json.dumps({"kd": result[0], "dg": result[1]}))
        except Exception:
            log.warning("Could not write PRODIGY disk cache: %s", cache_file)

    # ------------------------------------------------------------------
    # Score conversion
    # ------------------------------------------------------------------

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
    PRODIGY predictor for cases where a single fixed complex PDB/CIF is available.

    Useful for validating the parent peptide's Kd without running structure
    prediction.  Combine with ``BoltzStructurePredictor`` for a full pipeline.

    Args:
        complex_pdb: Path to the protein–peptide complex PDB or CIF.
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

    def evaluate_reference(self) -> _ProdigyResult:
        """
        Run PRODIGY on the fixed complex and return (Kd in M, ΔG in kcal/mol).

        Checks the disk cache first; sets ``self.reference_kd`` as a side-effect.
        """
        cached = self._load_disk_cache(self._fixed_complex_pdb)
        if cached is not None:
            result = cached
        else:
            result = self._run_prodigy(self._fixed_complex_pdb)
            if result is None:
                raise RuntimeError(f"PRODIGY failed on {self._fixed_complex_pdb}")
            self._save_disk_cache(self._fixed_complex_pdb, result)

        kd, _ = result
        if self.reference_kd is None:
            self.reference_kd = kd
        return result
