"""
Example: use PRODIGY as the binding affinity predictor in OptiCPP.

Two scenarios are shown:

  1. StaticProdigyBindingPredictor  – you already have a complex PDB (e.g.
     from the PDB, HADDOCK, or a previous docking run).  PRODIGY is run
     once to establish the reference Kd.  Combine this with StructuralFilter
     (AlphaFold2/Boltz) to generate per-variant structures and assess them.

  2. ProdigyBindingPredictor with a custom structure_fn  – you supply a
     callable that builds a complex PDB for any peptide sequence.  PRODIGY
     is called for every evaluated variant.

Requirements:
    pip install prodigy-prot biopython
"""

from pathlib import Path
from typing import Optional

# ---------------------------------------------------------------------------
# Scenario 1: fixed complex PDB (e.g. 1YCR – MDM2 / p53 peptide from RCSB)
# ---------------------------------------------------------------------------

def scenario_static(complex_pdb: Path):
    """Predict Kd from a known complex PDB, then run evolution."""
    from binder_evolution.binding.prodigy_predictor import StaticProdigyBindingPredictor

    predictor = StaticProdigyBindingPredictor(
        complex_pdb=complex_pdb,
        receptor_chains="A",   # MDM2 chain
        peptide_chains="B",    # p53 peptide chain
        temperature=25.0,
    )

    # Evaluate the fixed complex once to set the reference Kd
    kd, dg = predictor.evaluate_reference()
    print(f"Reference Kd : {kd:.3e} M")
    print(f"Reference ΔG : {dg:.2f} kcal/mol")

    # The predictor is now ready to be passed to EvolutionController.
    # All variants will receive score 1.0 unless their predicted structures
    # are evaluated separately (e.g. via StructuralFilter → PRODIGY pipeline).
    return predictor


# ---------------------------------------------------------------------------
# Scenario 2: structure_fn supplies per-variant complex PDBs
# ---------------------------------------------------------------------------

def make_structure_fn(target_pdb: Path):
    """
    Factory that returns a structure_fn compatible with ProdigyBindingPredictor.

    Replace the body of `_predict_complex` with your actual structure
    prediction tool (HADDOCK, AlphaFold-Multimer, Boltz, RoseTTAFold2, …).

    Args:
        target_pdb: PDB of the receptor (fixed across all variants).

    Returns:
        Callable(sequence, target_pdb) -> Path  pointing to the complex PDB.
    """

    def _predict_complex(sequence: str, target_pdb: Optional[Path]) -> Path:
        # --- replace this block with a real structure predictor ---
        # Example stub: write a minimal placeholder PDB so the code runs.
        import tempfile, textwrap
        tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
        tmp.write(textwrap.dedent(f"""\
            REMARK  Stub complex for sequence {sequence}
            END
        """).encode())
        tmp.close()
        return Path(tmp.name)
        # ----------------------------------------------------------

    return _predict_complex


def scenario_dynamic(target_pdb: Path):
    """Predict Kd for each variant via a user-supplied structure predictor."""
    from binder_evolution.binding.prodigy_predictor import ProdigyBindingPredictor

    structure_fn = make_structure_fn(target_pdb)

    predictor = ProdigyBindingPredictor(
        structure_fn=structure_fn,
        target_pdb=target_pdb,
        receptor_chains="A",
        peptide_chains="B",
        temperature=25.0,
        distance_cutoff=5.5,
        score_window=2.0,   # score drops 0→1 over 2 log10 units of Kd
    )

    # Plug into EvolutionController
    # from binder_evolution.config import EvolutionConfig
    # from binder_evolution.ml.cpp_scorer import create_cpp_scorer
    # from binder_evolution.evolution.controller import EvolutionController
    #
    # config = EvolutionConfig(binding_weight=0.6, cpp_weight=0.4)
    # cpp_scorer = create_cpp_scorer(Path("cpp_predictor_production.pkl"))
    # hotspots = {2, 6, 9}
    #
    # controller = EvolutionController(config, cpp_scorer, predictor, hotspots)
    # best, pop, history = controller.run("SQETFSDLWKLLPEN", "MDM2")

    return predictor


# ---------------------------------------------------------------------------
# Quick smoke test (no real structures needed)
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import math

    from binder_evolution.binding.prodigy_predictor import ProdigyBindingPredictor

    # Fake structure_fn that returns a hardcoded Kd via monkey-patching
    # (demonstrates score normalisation without real PDB files)
    def fake_structure_fn(sequence: str, target_pdb) -> Path:
        return Path("/dev/null")   # will be intercepted below

    predictor = ProdigyBindingPredictor(
        structure_fn=fake_structure_fn,
        receptor_chains="A",
        peptide_chains="B",
    )

    # Override _run_prodigy to return synthetic Kd values
    kd_map = {
        "SQETFSDLWKLLPEN": 1e-7,   # reference: 100 nM
        "RQRTFRDLWKRLPKN": 1e-7,   # same Kd → score 1.0
        "AAAAAAAAAAAAAAAA": 1e-5,   # 100× weaker → score 0.0
        "RQRTFADLWKALPKN": 5e-7,   # ~5× weaker → score ~0.65
    }

    def fake_run_prodigy(pdb):
        # Retrieve which sequence we are scoring (last call's sequence)
        return None  # intentionally return None; we'll call predict_kd below

    sequences = list(kd_map.keys())
    seed = sequences[0]

    # Manually inject cached Kd values to bypass structure/PRODIGY calls
    predictor._kd_cache = kd_map.copy()
    predictor.reference_kd = kd_map[seed]

    print(f"{'Sequence':<25} {'Kd (M)':>12} {'Score':>8}")
    print("-" * 48)
    for seq, kd in kd_map.items():
        score = predictor._kd_to_score(kd)
        print(f"{seq:<25} {kd:>12.2e} {score:>8.3f}")
