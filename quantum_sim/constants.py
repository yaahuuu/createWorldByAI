from __future__ import annotations

# Physical constants (arbitrary units scaled for stability and simplicity)
C = 299792458.0  # speed of light, not directly used but kept for reference
K_COULOMB = 0.05  # scaled Coulomb constant
G_NEWTON = 1e-4   # scaled gravitational constant
DT = 0.1          # time step (arbitrary units)
SOFTENING = 0.5   # softening length to avoid singularities in 1/r^2 forces
MAX_SPEED = 5.0   # cap speeds to keep numerical stability
WORLD_SIZE = 200.0  # size of the square world (periodic boundaries)

# Strong force parameters (toy model)
STRONG_ATTRACTION = 0.4  # attraction when color-neutralization possible
STRONG_RANGE = 5.0       # interaction range for strong force
BIND_DISTANCE = 2.0      # threshold to consider particles bonded
BOND_K = 0.8             # spring constant for bonds
BOND_DAMPING = 0.2       # damping in bonds

# Quantum-like stochastic jitter (toy Heisenberg-style diffusion)
JITTER_STD = 0.05

# Macro detection parameters
CLUSTER_NEIGHBOR_RADIUS = 3.0
MACRO_CLUSTER_THRESHOLD = 30  # number of particles that counts as macro cluster
COARSE_GRID_SIZE = 20  # coarse-grain grid dimension per axis for macro fields
