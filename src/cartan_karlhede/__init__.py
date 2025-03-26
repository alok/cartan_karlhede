"""Cartan-Karlhede algorithm for determining the equivalence of spacetimes in general relativity."""

from cartan_karlhede.algorithm import compare_metrics
from cartan_karlhede.metric import Metric, LorentzFrame, NullFrame
from cartan_karlhede.examples import compare_schwarzschild_examples, compare_3d_examples
from cartan_karlhede.simple_examples import run_simple_examples
from cartan_karlhede.cli import main as cli_main


def main() -> None:
    """Run the Cartan-Karlhede algorithm examples.

    In this simplified version, we only run basic examples to avoid complex
    calculations with frames that might be difficult to compute correctly
    in this initial implementation.
    """
    # Run simple examples to verify the implementation
    run_simple_examples()
