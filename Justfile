# Justfile for the Cartan-Karlhede algorithm

# List available commands
default:
    @just --list

# Install the package in development mode
install:
    uv pip install -e ".[dev]"

# Install documentation dependencies
install-docs:
    uv pip install sphinx sphinx_rtd_theme

# Run the examples
run-examples:
    python -m cartan_karlhede --verbose

# Run the 3D examples
run-3d:
    python -m cartan_karlhede --example 3d --verbose

# Run the Schwarzschild examples
run-schwarzschild:
    python -m cartan_karlhede --example schwarzschild --verbose

# Run tests
test:
    python -m pytest

# Run type checking
typecheck:
    python -m mypy src

# Format the code
format:
    python -m black src
    python -m isort src

# Run linting checks
lint:
    python -m flake8 src
    python -m black --check src
    python -m isort --check src
    python -m mypy src

# Build the documentation
build-docs:
    cd docs && make html

# Clean the documentation build directory
clean-docs:
    cd docs && make clean

# Open the documentation in a browser (macOS)
open-docs:
    open docs/_build/html/index.html

# Full documentation workflow: install dependencies, build, and open
docs: install-docs build-docs open-docs 