# API Documentation

This document provides detailed API documentation for the CPP Binder Evolution framework.

## Core Modules

### binder_evolution.ml.cpp_scorer

Machine learning-based CPP prediction module.

#### CPPScorer

Main class for CPP probability prediction using the production ML model.

```python
class CPPScorer:
    """
    Production ML-based CPP scorer.
    
    Uses ensemble of Random Forest + Gradient Boosting + SVM
    trained on 1,075 experimentally validated CPPs.
    
    Attributes:
        predictor: MLCPPPredictor instance
        
    Performance:
        Accuracy: 89%
        AUC: 0.95
        Precision: 87%
        Recall: 91%
    """
    
    def __init__(self, model_path: Path):
        """
        Initialize CPP scorer.
        
        Args:
            model_path: Path to cpp_predictor_production.pkl
            
        Raises:
            FileNotFoundError: If model file not found
            ImportError: If ml_cpp_predictor module not available
        """
    
    def score(self, sequence: str) -> float:
        """
        Predict CPP probability.
        
        Args:
            sequence: Amino acid sequence (single-letter code)
            
        Returns:
            CPP probability (0-1), higher indicates stronger CPP properties
            
        Example:
            >>> scorer = CPPScorer(Path("cpp_predictor_production.pkl"))
            >>> prob = scorer.score("RQIKIWFQNRRMKWKK")
            >>> print(f"CPP probability: {prob:.3f}")
            CPP probability: 0.999
        """
    
    def predict(self, sequence: str) -> CPPPrediction:
        """
        Get full CPP prediction with details.
        
        Args:
            sequence: Amino acid sequence
            
        Returns:
            CPPPrediction object containing:
                - cpp_probability: 0-1 score
                - confidence: prediction confidence
                - cpp_type: 'arginine-rich', 'amphipathic', 'hydrophobic', 'mixed'
                - feature_importances: dict of feature contributions
                
        Example:
            >>> pred = scorer.predict("RQIKIWFQNRRMKWKK")
            >>> print(f"Type: {pred.cpp_type}")
            Type: mixed
        """
    
    def is_cpp(self, sequence: str, threshold: float = 0.5) -> bool:
        """
        Binary CPP classification.
        
        Args:
            sequence: Amino acid sequence
            threshold: CPP probability threshold (default 0.5)
            
        Returns:
            True if predicted as CPP, False otherwise
        """
```

#### Factory Function

```python
def create_cpp_scorer(model_path: Path = None, use_fallback: bool = False) -> CPPScorer:
    """
    Factory function to create CPP scorer.
    
    Args:
        model_path: Path to production model (default: cpp_predictor_production.pkl)
        use_fallback: If True, use heuristic fallback (65% accurate)
        
    Returns:
        CPPScorer or CPPScorerFallback instance
        
    Example:
        >>> scorer = create_cpp_scorer(Path("model.pkl"))
        >>> # or for testing without model:
        >>> scorer = create_cpp_scorer(use_fallback=True)
    """
```

### binder_evolution.evolution.controller

Main evolution controller for CPP optimization.

#### EvolutionController

```python
class EvolutionController:
    """
    Main evolution controller for CPP optimization.
    
    Evolves a known binding peptide to add CPP properties
    while maintaining binding affinity.
    
    Attributes:
        config: EvolutionConfig
        cpp_scorer: CPPScorer
        binding_predictor: User-provided binding predictor
        binding_hotspots: Set of protected positions
        structural_filter: Optional StructuralFilter
        generator: VariantGenerator
        history: List of fitness evaluations
        best_fitness_history: Fitness over generations
        best_peptide_history: Best peptides over generations
    """
    
    def __init__(
        self,
        config: EvolutionConfig,
        cpp_scorer: CPPScorer,
        binding_predictor,
        binding_hotspots: Set[int] = None,
        structural_filter = None
    ):
        """
        Initialize evolution controller.
        
        Args:
            config: Evolution configuration
            cpp_scorer: CPP probability scorer (production ML model)
            binding_predictor: Binding affinity predictor with predict_binding(seq) method
            binding_hotspots: Set of binding-critical positions to protect
            structural_filter: Optional StructuralFilter for pose validation
            
        Raises:
            RuntimeError: If cpp_scorer or binding_predictor is None
        """
    
    def run(
        self,
        seed_sequence: str,
        target_protein: str = None
    ) -> Tuple[str, List[str], List[Dict]]:
        """
        Run evolutionary optimization.
        
        Args:
            seed_sequence: Known binding peptide to optimize
            target_protein: Target protein identifier (optional)
            
        Returns:
            Tuple of:
                - best_peptide: Optimized peptide sequence
                - final_population: Last generation population
                - history: Complete evolution history
                
        Example:
            >>> controller = EvolutionController(
            ...     config, cpp_scorer, binding_predictor, {2, 6, 9}
            ... )
            >>> best, pop, history = controller.run("SQETFSDLWKLLPEN", "MDM2")
            >>> print(f"Optimized: {best}")
            Optimized: RQRTFRDLWKRLPKN
        """
    
    def export_results(self, output_path: Path):
        """
        Export optimization results to JSON file.
        
        Args:
            output_path: Path to output file
            
        Creates:
            JSON file with:
                - config: Configuration used
                - best_peptide: Final optimized peptide
                - best_fitness: Final fitness value
                - generations: Number of generations
                - history: Evaluation history (last 100 records)
        """
```

### binder_evolution.evolution.generator

Variant generation with hotspot protection.

#### VariantGenerator

```python
class VariantGenerator:
    """
    Generate peptide variants with smart mutations.
    
    Features:
        - Protects binding hotspots
        - Biases towards CPP-friendly residues
        - Implements crossover
        - Conservative mutations for hotspots if needed
        
    Attributes:
        binding_hotspots: Set of protected positions
    """
    
    def __init__(self, seed: int, binding_hotspots: Set[int] = None):
        """
        Initialize variant generator.
        
        Args:
            seed: Random seed for reproducibility
            binding_hotspots: Set of position indices that are binding-critical
        """
    
    def mutate(
        self,
        sequence: str,
        mutation_rate: float = 0.15,
        cpp_bias: float = 0.7,
        allow_hotspot_mutations: bool = False
    ) -> str:
        """
        Mutate sequence while protecting binding hotspots.
        
        Args:
            sequence: Original sequence
            mutation_rate: Probability of mutating each position
            cpp_bias: Probability of using CPP-friendly amino acid
            allow_hotspot_mutations: If True, allow conservative mutations in hotspots
            
        Returns:
            Mutated sequence
            
        Example:
            >>> gen = VariantGenerator(seed=42, binding_hotspots={2, 6, 9})
            >>> variant = gen.mutate("SQETFSDLWKLLPEN")
            >>> # Positions 2, 6, 9 protected, others may be mutated
        """
    
    def crossover(
        self,
        parent1: str,
        parent2: str,
        method: str = 'single_point'
    ) -> Tuple[str, str]:
        """
        Perform crossover between two parent sequences.
        
        Args:
            parent1: First parent sequence
            parent2: Second parent sequence
            method: 'single_point', 'two_point', or 'uniform'
            
        Returns:
            Tuple of (child1, child2)
            
        Raises:
            ValueError: If parents have different lengths or invalid method
            
        Example:
            >>> child1, child2 = gen.crossover("RQIKI...", "GRKKR...")
        """
    
    def generate_variants(
        self,
        parent: str,
        n_variants: int,
        mutation_rate: float = 0.15
    ) -> List[str]:
        """
        Generate multiple variants from a parent.
        
        Args:
            parent: Parent sequence
            n_variants: Number of variants to generate
            mutation_rate: Mutation rate
            
        Returns:
            List of variant sequences
        """
```

#### Utility Functions

```python
def calculate_diversity(population: List[str]) -> float:
    """
    Calculate population diversity (average pairwise dissimilarity).
    
    Args:
        population: List of sequences
        
    Returns:
        Diversity score (0-1), higher = more diverse
        
    Example:
        >>> pop = ["AAA", "AAB", "AAC"]
        >>> div = calculate_diversity(pop)
        >>> print(f"Diversity: {div:.2f}")
        Diversity: 0.67
    """

def maintain_diversity(
    population: List[str],
    seed_sequence: str,
    generator: VariantGenerator,
    diversity_threshold: float = 0.15,
    injection_rate: float = 0.1
) -> List[str]:
    """
    Maintain population diversity by injecting random variants.
    
    Args:
        population: Current population
        seed_sequence: Original sequence to mutate from
        generator: VariantGenerator instance
        diversity_threshold: Minimum acceptable diversity
        injection_rate: Fraction of population to replace
        
    Returns:
        Population with diversity maintained
    """
```

### binder_evolution.structure.structural_filter

Structural validation for binding pose preservation.

#### StructuralFilter

```python
class StructuralFilter:
    """
    Filter variants based on structural validation.
    
    This implements Step 3.5 of the workflow.
    
    Attributes:
        predictor: StructurePredictor
        analyzer: BindingPoseAnalyzer
        min_contact_overlap: Minimum fraction of parent contacts
        max_rmsd: Maximum RMSD from parent (Angstroms)
        require_pocket_match: Require binding in same pocket
        parent_interface: Cached parent binding interface
        parent_structure: Cached parent structure
    """
    
    def __init__(
        self,
        structure_predictor: StructurePredictor,
        pose_analyzer: BindingPoseAnalyzer,
        min_contact_overlap: float = 0.7,
        max_rmsd: float = 3.0,
        require_pocket_match: bool = True
    ):
        """
        Initialize structural filter.
        
        Args:
            structure_predictor: Structure prediction interface
            pose_analyzer: Binding pose analyzer
            min_contact_overlap: Minimum fraction of parent contacts to keep (default 0.7)
            max_rmsd: Maximum RMSD from parent in Angstroms (default 3.0)
            require_pocket_match: Require binding in same pocket (default True)
        """
    
    def set_parent_complex(
        self,
        parent_sequence: str,
        target_sequence: str = None,
        target_pdb: Path = None
    ):
        """
        Set parent complex as reference.
        
        Args:
            parent_sequence: Parent peptide sequence
            target_sequence: Target protein sequence (if predicting from sequence)
            target_pdb: Target protein structure file (if using known structure)
            
        Note:
            This must be called before validate_variant() or filter_population()
        """
    
    def validate_variant(self, variant_sequence: str) -> Tuple[bool, PoseComparison]:
        """
        Validate variant binding pose against parent.
        
        Args:
            variant_sequence: Variant peptide sequence
            
        Returns:
            Tuple of (is_valid, comparison_results)
            
        Raises:
            RuntimeError: If parent complex not set
            
        Example:
            >>> filter.set_parent_complex("SQETFSDLWKLLPEN", target_pdb=Path("mdm2.pdb"))
            >>> is_valid, comparison = filter.validate_variant("RQETFSDLWKLLPEN")
            >>> if is_valid:
            ...     print(f"Valid pose, overlap: {comparison.contact_overlap:.2f}")
        """
    
    def filter_population(
        self,
        variants: List[str]
    ) -> Tuple[List[str], List[PoseComparison]]:
        """
        Filter population based on structural validation.
        
        Args:
            variants: List of variant sequences
            
        Returns:
            Tuple of (valid_variants, comparison_results)
            
        Example:
            >>> variants = ["RQETF...", "GRKKR...", ...]
            >>> valid, comparisons = filter.filter_population(variants)
            >>> print(f"Valid: {len(valid)}/{len(variants)}")
            Valid: 75/100
        """
```

#### PoseComparison

```python
@dataclass
class PoseComparison:
    """
    Results of comparing two binding poses.
    
    Attributes:
        contact_overlap: Fraction of parent contacts maintained (0-1)
        rmsd: Root mean square deviation in Angstroms
        pocket_match: True if binds in same pocket
        pose_similarity: Overall similarity score (0-1)
        critical_contacts_kept: Set of hotspot positions still binding
        new_contacts: List of new interactions formed
        lost_contacts: List of parent interactions lost
        energy_change: Interaction energy change from parent (kcal/mol)
    """
```

### binder_evolution.structure.cpp_structures

CPPsite2.0 PDB structure database integration.

#### CPPStructureDatabase

```python
class CPPStructureDatabase:
    """
    Database of CPPsite2.0 experimental structures.
    
    Manages 58 PDB structures of experimentally validated CPPs.
    
    Attributes:
        pdb_dir: Directory containing PDB files
        structures: Dict mapping PDB ID to CPPStructure
        sequence_to_pdb: Dict mapping sequence to PDB ID
        parser: PDBParser instance
    """
    
    def __init__(self, pdb_dir: Path):
        """
        Initialize structure database.
        
        Args:
            pdb_dir: Directory containing CPPsite2.0 PDB files
            
        Example:
            >>> db = CPPStructureDatabase(Path("pdb_structures"))
            >>> print(f"Loaded {len(db.structures)} structures")
            Loaded 58 structures
        """
    
    def has_structure(self, sequence: str) -> bool:
        """
        Check if exact structure exists for sequence.
        
        Args:
            sequence: Amino acid sequence
            
        Returns:
            True if experimental structure available
        """
    
    def get_structure(self, sequence: str) -> Optional[CPPStructure]:
        """
        Get structure for exact sequence match.
        
        Args:
            sequence: Amino acid sequence
            
        Returns:
            CPPStructure if available, None otherwise
            
        Example:
            >>> struct = db.get_structure("RQIKIWFQNRRMKWKK")
            >>> if struct:
            ...     print(f"PDB ID: {struct.pdb_id}")
        """
    
    def find_similar_structure(
        self,
        sequence: str,
        min_similarity: float = 0.8
    ) -> Optional[CPPStructure]:
        """
        Find most similar structure by sequence similarity.
        
        Args:
            sequence: Target sequence
            min_similarity: Minimum sequence identity threshold
            
        Returns:
            Most similar CPPStructure if above threshold, None otherwise
        """
    
    def get_statistics(self) -> dict:
        """
        Get database statistics.
        
        Returns:
            Dict with:
                - total_structures: Number of structures
                - unique_sequences: Number of unique sequences
                - min_length: Minimum peptide length
                - max_length: Maximum peptide length
                - mean_length: Average peptide length
                - total_residues: Total number of residues
        """
```

### binder_evolution.config

Configuration management.

#### EvolutionConfig

```python
@dataclass
class EvolutionConfig:
    """
    Configuration for evolutionary optimization.
    
    All parameters are validated and have sensible defaults.
    """
    
    # Random seed
    seed: int = 42
    
    # Evolution parameters
    max_generations: int = 500
    population_size: int = 100
    parents_kept: int = 20
    
    # Mutation parameters
    mutation_rate: float = 0.15
    crossover_rate: float = 0.7
    tournament_size: int = 3
    
    # Fitness weights
    binding_weight: float = 0.6
    cpp_weight: float = 0.4
    
    # Constraints
    min_binding_score: float = 0.7
    binding_tolerance: float = 0.2
    
    # Diversity maintenance
    diversity_threshold: float = 0.85
    diversity_injection_rate: float = 0.1
    
    # Convergence detection
    convergence_window: int = 50
    convergence_threshold: float = 0.01
    
    # Paths
    cpp_model_path: Path = Path("cpp_predictor_production.pkl")
    pdb_structures_dir: Path = Path("pdb_structures")
    
    # GPU
    use_gpu: bool = True
```

## Usage Examples

### Basic CPP Optimization

```python
from pathlib import Path
from binder_evolution.config import EvolutionConfig
from binder_evolution.ml.cpp_scorer import create_cpp_scorer
from binder_evolution.evolution.controller import EvolutionController

# Configure
config = EvolutionConfig(
    seed=42,
    max_generations=200,
    population_size=50
)

# Initialize
cpp_scorer = create_cpp_scorer(Path("cpp_predictor_production.pkl"))
binding_predictor = MyBindingPredictor()
hotspots = {2, 6, 9}

# Run
controller = EvolutionController(config, cpp_scorer, binding_predictor, hotspots)
best, pop, history = controller.run("SQETFSDLWKLLPEN")
```

### With Structural Validation

```python
from binder_evolution.structure.structural_filter import (
    StructurePredictor, BindingPoseAnalyzer, StructuralFilter
)

# Setup structure prediction
predictor = StructurePredictor("alphafold2")
analyzer = BindingPoseAnalyzer(contact_distance=4.5)
structural_filter = StructuralFilter(predictor, analyzer)

# Run with validation
controller = EvolutionController(
    config, cpp_scorer, binding_predictor, hotspots,
    structural_filter=structural_filter
)
best, pop, history = controller.run("SQETFSDLWKLLPEN")
```

### Batch Processing

```python
# Optimize multiple peptides
peptides = ["SQETFSDLWKLLPEN", "LTFEHYWAQLTS", "EQMTLRLLQE"]
results = {}

for peptide in peptides:
    best, _, _ = controller.run(peptide)
    results[peptide] = best

# Export all results
import json
with open("batch_results.json", "w") as f:
    json.dump(results, f, indent=2)
```

## Error Handling

All methods raise appropriate exceptions:

- `ValueError`: Invalid parameters
- `FileNotFoundError`: Missing required files
- `RuntimeError`: Invalid state (e.g., parent complex not set)
- `TypeError`: Type mismatches

Example:

```python
try:
    controller.run(seed_sequence)
except ValueError as e:
    print(f"Invalid parameter: {e}")
except FileNotFoundError as e:
    print(f"Missing file: {e}")
except RuntimeError as e:
    print(f"Runtime error: {e}")
```

## Performance Considerations

- Use GPU for structure prediction when available
- Cache expensive computations
- Batch predictions when possible
- Monitor memory usage for large populations
- Consider parallel evaluation for independent variants

## Version Compatibility

This API is stable as of v1.0.0. Breaking changes will increment major version.

For implementation details, see source code documentation.
