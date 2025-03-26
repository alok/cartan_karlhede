"""Command-line interface for the Cartan-Karlhede algorithm."""

import argparse
import time
from typing import List

from cartan_karlhede.examples import (
    compare_schwarzschild_examples,
    compare_3d_examples,
)


def parse_args(args: List[str] = None) -> argparse.Namespace:
    """Parse command-line arguments.

    Args:
        args: Command-line arguments (defaults to sys.argv)

    Returns:
        Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="Cartan-Karlhede algorithm for spacetime equivalence"
    )

    parser.add_argument(
        "--example",
        "-e",
        choices=["schwarzschild", "3d", "all"],
        default="all",
        help="Example to run (default: all)",
    )

    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Enable verbose output"
    )

    return parser.parse_args(args)


def main(args: List[str] = None) -> None:
    """Run the specified examples.

    Args:
        args: Command-line arguments
    """
    parsed_args = parse_args(args)

    if parsed_args.verbose:
        print("Running Cartan-Karlhede algorithm with verbose output")

    if parsed_args.example in ["schwarzschild", "all"]:
        if parsed_args.verbose:
            print("\nRunning Schwarzschild examples...")
            start_time = time.time()

        compare_schwarzschild_examples()

        if parsed_args.verbose:
            elapsed = time.time() - start_time
            print(f"Completed in {elapsed:.2f} seconds")

    if parsed_args.example in ["3d", "all"]:
        if parsed_args.verbose:
            print("\nRunning 3D examples...")
            start_time = time.time()

        compare_3d_examples()

        if parsed_args.verbose:
            elapsed = time.time() - start_time
            print(f"Completed in {elapsed:.2f} seconds")


if __name__ == "__main__":
    main()
