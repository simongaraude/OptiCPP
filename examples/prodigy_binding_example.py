"""
Example: use PRODIGY + Boltz as the binding affinity predictor in OptiCPP.

Two scenarios are shown:

  1. StaticProdigyBindingPredictor  – you already have a complex PDB/CIF (e.g.
     from the RCSB PDB or a previous docking run).  PRODIGY is run once on the
     fixed complex.

  2. BoltzStructurePredictor + ProdigyBindingPredictor  – Boltz generates a
     fresh complex CIF for every peptide variant; PRODIGY scores each one.
     This is the recommended pipeline for evolutionary optimization.

Requirements:
    pip install boltz prodigy-prot biopython
"""

from pathlib import Path

# ---------------------------------------------------------------------------
# Scenario 1: fixed complex (PDB or CIF)
# ---------------------------------------------------------------------------

def scenario_static(complex_path: Path):
    """Predict Kd from a known complex PDB/CIF, then run evolution."""
    from binder_evolution.binding.prodigy_predictor import StaticProdigyBindingPredictor

    predictor = StaticProdigyBindingPredictor(
        complex_pdb=complex_path,
        receptor_chains="A",   # receptor chain ID in the file
        peptide_chains="B",    # peptide chain ID in the file
        temperature=25.0,
    )

    kd, dg = predictor.evaluate_reference()
    print(f"Reference Kd : {kd:.3e} M")
    print(f"Reference ΔG : {dg:.2f} kcal/mol")
    return predictor


# ---------------------------------------------------------------------------
# Scenario 2: Boltz → PRODIGY pipeline (recommended for evolution)
# ---------------------------------------------------------------------------

def scenario_boltz(receptor_sequence: str, seed_peptide: str):
    """
    Generate per-variant complex structures with Boltz, score with PRODIGY.

    Args:
        receptor_sequence: Full amino-acid sequence of the target protein.
        seed_peptide:      Seed CPP sequence.
    """
    from binder_evolution.binding.prodigy_predictor import ProdigyBindingPredictor
    from binder_evolution.structure.boltz_predictor import BoltzStructurePredictor

    boltz = BoltzStructurePredictor(
        receptor_sequence=receptor_sequence,
        output_dir=Path("boltz_out"),
        receptor_chain="A",
        peptide_chain="B",
        recycling_steps=3,
        diffusion_samples=1,
        use_msa_server=True,   # auto MSA via mmseqs2 (requires internet)
    )

    predictor = ProdigyBindingPredictor(
        structure_fn=boltz,        # BoltzStructurePredictor is callable
        receptor_chains="A",
        peptide_chains="B",
        temperature=25.0,
        distance_cutoff=5.5,
    )

    seed_dg = predictor.predict_binding(seed_peptide)
    print(f"Seed peptide : {seed_peptide}")
    print(f"Seed ΔG      : {seed_dg:.2f} kcal/mol")

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
    # best, pop, history = controller.run(seed_peptide, "MDM2")

    return predictor


# ---------------------------------------------------------------------------
# Quick smoke test (no real Boltz/PRODIGY installation needed)
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    from binder_evolution.binding.prodigy_predictor import ProdigyBindingPredictor

    def fake_structure_fn(sequence, target_pdb=None):
        return Path("/dev/null")   # bypassed below

    predictor = ProdigyBindingPredictor(
        structure_fn=fake_structure_fn,
        receptor_chains="A",
        peptide_chains="B",
    )

    dg_map = {
        "SQETFSDLWKLLPEN": -9.5,   # strong binder
        "RQRTFRDLWKRLPKN": -9.5,   # same ΔG
        "RQRTFADLWKALPKN": -8.2,   # weaker binder
        "AAAAAAAAAAAAAAAA": -4.0,  # very weak binder
    }

    # Inject fake results directly into the in-memory cache
    predictor._result_cache = {
        seq: (0.0, dg) for seq, dg in dg_map.items()
    }

    print(f"{'Sequence':<25} {'ΔG (kcal/mol)':>14}")
    print("-" * 42)
    for seq in dg_map:
        dg = predictor.predict_binding(seq)
        print(f"{seq:<25} {dg:>14.2f}")
