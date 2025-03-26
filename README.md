# Cartan-Karlhede Algorithm

Implementation of the Cartan-Karlhede algorithm for determining the equivalence of spacetimes in general relativity.

## Description

This project implements the Cartan-Karlhede algorithm, which provides a systematic procedure to determine whether two spacetime metrics represent the same physical geometry, despite being expressed in different coordinate systems. The algorithm works by computing and comparing the functional relations between curvature tensors and their covariant derivatives in appropriately chosen frames.

Based on the paper "On the Equivalence of Spacetimes, the Cartan-Karlhede Algorithm" by Thiago M. Mergulhão and Carlos Batista.

## Features

- Representation of metrics and frame systems
- Computation of curvature tensors (Riemann, Ricci, Weyl)
- Computation of curvature derivatives
- Framework for comparing metrics using the Cartan-Karlhede algorithm
- Support for both Lorentzian frames and null frames (4D)
- Example implementations from the paper

## Current Status

The current implementation provides the fundamental building blocks of the Cartan-Karlhede algorithm:

- ✅ Metric representation and tensor calculations
- ✅ Curvature tensor calculations (Riemann, Ricci, Weyl)
- ✅ Basic frame support (Lorentzian and null frames)
- ✅ Basic functionality for comparing metrics
- ✅ Simple examples for testing

The following features are still being developed:

- ❌ Robust implementation of orthonormal frames for complex metrics
- ❌ Correct computation of functional relations between tensor components
- ❌ Complete implementation of isotropy group calculations
- ❌ Advanced examples from the paper (3D and Schwarzschild)

## Installation

This project uses `uv` for package management. To install:

```bash
# Clone the repository
git clone https://github.com/yourusername/cartan-karlhede.git
cd cartan-karlhede

# Create a virtual environment and install dependencies
uv venv .venv
source .venv/bin/activate
uv pip install -e .
```

## Usage

```python
from cartan_karlhede.metric import Metric
from sympy import Matrix, Symbol, diag

# Create a simple metric
t = Symbol('t')
x = Symbol('x')
y = Symbol('y')
z = Symbol('z')
g = diag(-1, 1, 1, 1)
metric = Metric(g, [t, x, y, z], "Minkowski 4D")

# Compute curvature
riemann = metric.riemann_tensor()
ricci = metric.ricci_tensor()
ricci_scalar = metric.ricci_scalar()
```

For more examples, see the `simple_examples.py` module.

## Implementation Notes

This implementation focuses on the core conceptual aspects of the Cartan-Karlhede algorithm. Several simplifications have been made:

1. Computing the minimum value of t_i (τ_i) requires testing all possible frames, which is computationally intensive. This implementation approximates τ_i as t_i.

2. Finding the isotropy group at each step requires sophisticated group theory algorithms. This implementation provides a framework where these groups could be calculated, but doesn't fully implement the computation.

3. The covariant derivatives are computed in a simplified way that works for coordinate-aligned frames but would need to be expanded for general frames.

4. The frame validation and orthogonalization need refinement for complex metrics.

## Future Improvements

1. Improving the frame construction for complex metrics
2. Implementing robust calculation of functional relations between tensor components
3. Adding complete isotropy group calculations
4. Creating a more sophisticated method for computing covariant derivatives
5. Expanding the examples from the paper with corrected implementations
6. Adding visualization tools for curvature components and functional relations

## License

MIT

## References

1. Mergulhão, T. M., & Batista, C. (2023). On the Equivalence of Spacetimes, the Cartan-Karlhede Algorithm.

2. Karlhede, A. (1980). A Review of the Geometrical Equivalence of Metrics in General Relativity. General Relativity and Gravitation, 12, 693-707.

3. Cartan, É. (1951). Leçons sur la geometrie des espaces de Riemann. Gauthier-Villars Paris.
