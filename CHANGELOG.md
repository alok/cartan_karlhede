# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2023-07-31

### Added
- Initial implementation of the Cartan-Karlhede algorithm core components
- Metric class for representing spacetime metrics
- Tensor calculations (Riemann, Ricci, Weyl)
- Frame classes (Lorentzian and null frames)
- Curvature component and derivative computations
- Algorithm framework for comparing metrics
- Simple examples with Minkowski spacetime

### Known Issues
- Frame construction for complex metrics needs improvement
- Functional relation computation between tensor components is incomplete
- Isotropy group calculations are simplified
- Examples from the paper (3D and Schwarzschild) need refinement 