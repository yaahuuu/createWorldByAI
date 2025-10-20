"""
Quantum-inspired, educational particle simulation framework.
This package provides a simplified 2D continuous-space simulation whose
microscopic entities are 19 canonical particle types (12 fermions + 7 bosons).
It does NOT implement full quantum field theory; instead, it enforces a set of
coarse conservation rules and probabilistic interactions to explore whether
macro-like phenomena can emerge from micro rules.
"""

__all__ = [
    "ParticleType",
    "Particle",
    "Universe",
]

from .particles import ParticleType, Particle  # noqa: E402
from .world import Universe  # noqa: E402
