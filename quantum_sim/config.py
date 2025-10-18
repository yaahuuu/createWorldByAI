from __future__ import annotations

from typing import Dict

from .particles import ParticleType

# Default initial counts for each particle type.
# These are intentionally modest to keep the simulation responsive.
DEFAULT_COUNTS: Dict[ParticleType, int] = {
    # quarks — enough to form a few baryons
    ParticleType.UP: 60,
    ParticleType.DOWN: 60,
    ParticleType.CHARM: 4,
    ParticleType.STRANGE: 6,
    ParticleType.TOP: 0,       # top is heavy and short-lived; omit by default
    ParticleType.BOTTOM: 2,

    # leptons
    ParticleType.ELECTRON: 40,
    ParticleType.MUON: 6,
    ParticleType.TAU: 2,
    ParticleType.ELECTRON_NEUTRINO: 40,
    ParticleType.MUON_NEUTRINO: 20,
    ParticleType.TAU_NEUTRINO: 10,

    # bosons (as free excitations - optional)
    ParticleType.PHOTON: 30,
    ParticleType.GLUON: 10,
    ParticleType.W_PLUS: 0,
    ParticleType.W_MINUS: 0,
    ParticleType.Z_BOSON: 0,
    ParticleType.HIGGS: 0,
    ParticleType.GRAVITON: 0,
}
