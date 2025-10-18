from __future__ import annotations

import math
import random
from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, Optional, Tuple

Color = Optional[str]  # 'r','g','b' for quarks, None otherwise


class ParticleType(Enum):
    # Fermions (12)
    UP = "up_quark"
    DOWN = "down_quark"
    CHARM = "charm_quark"
    STRANGE = "strange_quark"
    TOP = "top_quark"
    BOTTOM = "bottom_quark"

    ELECTRON = "electron"
    MUON = "muon"
    TAU = "tau"
    ELECTRON_NEUTRINO = "electron_neutrino"
    MUON_NEUTRINO = "muon_neutrino"
    TAU_NEUTRINO = "tau_neutrino"

    # Bosons (7) — treated as stable mediators in this toy model
    PHOTON = "photon"
    GLUON = "gluon"
    W_PLUS = "W_plus"
    W_MINUS = "W_minus"
    Z_BOSON = "Z_boson"
    HIGGS = "higgs"
    GRAVITON = "graviton"


# Masses (MeV/c^2, approximate), but we use them as arbitrary units
# Simplified and rounded for stability
PARTICLE_MASS: Dict[ParticleType, float] = {
    # quarks (current masses; toy values)
    ParticleType.UP: 2.3,
    ParticleType.DOWN: 4.8,
    ParticleType.CHARM: 1275.0,
    ParticleType.STRANGE: 95.0,
    ParticleType.TOP: 173_000.0,
    ParticleType.BOTTOM: 4180.0,
    # leptons
    ParticleType.ELECTRON: 0.511,
    ParticleType.MUON: 105.7,
    ParticleType.TAU: 1776.86,
    ParticleType.ELECTRON_NEUTRINO: 0.0,
    ParticleType.MUON_NEUTRINO: 0.0,
    ParticleType.TAU_NEUTRINO: 0.0,
    # bosons
    ParticleType.PHOTON: 0.0,
    ParticleType.GLUON: 0.0,
    ParticleType.W_PLUS: 80_379.0,
    ParticleType.W_MINUS: 80_379.0,
    ParticleType.Z_BOSON: 91_187.6,
    ParticleType.HIGGS: 125_100.0,
    ParticleType.GRAVITON: 0.0,
}

# Electric charge in units of e
PARTICLE_CHARGE: Dict[ParticleType, float] = {
    ParticleType.UP: +2.0 / 3.0,
    ParticleType.CHARM: +2.0 / 3.0,
    ParticleType.TOP: +2.0 / 3.0,
    ParticleType.DOWN: -1.0 / 3.0,
    ParticleType.STRANGE: -1.0 / 3.0,
    ParticleType.BOTTOM: -1.0 / 3.0,

    ParticleType.ELECTRON: -1.0,
    ParticleType.MUON: -1.0,
    ParticleType.TAU: -1.0,
    ParticleType.ELECTRON_NEUTRINO: 0.0,
    ParticleType.MUON_NEUTRINO: 0.0,
    ParticleType.TAU_NEUTRINO: 0.0,

    ParticleType.PHOTON: 0.0,
    ParticleType.GLUON: 0.0,
    ParticleType.W_PLUS: +1.0,
    ParticleType.W_MINUS: -1.0,
    ParticleType.Z_BOSON: 0.0,
    ParticleType.HIGGS: 0.0,
    ParticleType.GRAVITON: 0.0,
}

# Spin (in units of ħ), simplified
PARTICLE_SPIN: Dict[ParticleType, float] = {
    ParticleType.UP: 0.5,
    ParticleType.DOWN: 0.5,
    ParticleType.CHARM: 0.5,
    ParticleType.STRANGE: 0.5,
    ParticleType.TOP: 0.5,
    ParticleType.BOTTOM: 0.5,

    ParticleType.ELECTRON: 0.5,
    ParticleType.MUON: 0.5,
    ParticleType.TAU: 0.5,
    ParticleType.ELECTRON_NEUTRINO: 0.5,
    ParticleType.MUON_NEUTRINO: 0.5,
    ParticleType.TAU_NEUTRINO: 0.5,

    ParticleType.PHOTON: 1.0,
    ParticleType.GLUON: 1.0,
    ParticleType.W_PLUS: 1.0,
    ParticleType.W_MINUS: 1.0,
    ParticleType.Z_BOSON: 1.0,
    ParticleType.HIGGS: 0.0,
    ParticleType.GRAVITON: 2.0,  # hypothetical
}

# Baryon and lepton numbers
PARTICLE_BNUMBER: Dict[ParticleType, float] = {
    ParticleType.UP: 1.0 / 3.0,
    ParticleType.DOWN: 1.0 / 3.0,
    ParticleType.CHARM: 1.0 / 3.0,
    ParticleType.STRANGE: 1.0 / 3.0,
    ParticleType.TOP: 1.0 / 3.0,
    ParticleType.BOTTOM: 1.0 / 3.0,

    ParticleType.ELECTRON: 0.0,
    ParticleType.MUON: 0.0,
    ParticleType.TAU: 0.0,
    ParticleType.ELECTRON_NEUTRINO: 0.0,
    ParticleType.MUON_NEUTRINO: 0.0,
    ParticleType.TAU_NEUTRINO: 0.0,

    ParticleType.PHOTON: 0.0,
    ParticleType.GLUON: 0.0,
    ParticleType.W_PLUS: 0.0,
    ParticleType.W_MINUS: 0.0,
    ParticleType.Z_BOSON: 0.0,
    ParticleType.HIGGS: 0.0,
    ParticleType.GRAVITON: 0.0,
}

PARTICLE_LNUMBER: Dict[ParticleType, float] = {
    ParticleType.UP: 0.0,
    ParticleType.DOWN: 0.0,
    ParticleType.CHARM: 0.0,
    ParticleType.STRANGE: 0.0,
    ParticleType.TOP: 0.0,
    ParticleType.BOTTOM: 0.0,

    ParticleType.ELECTRON: 1.0,
    ParticleType.MUON: 1.0,
    ParticleType.TAU: 1.0,
    ParticleType.ELECTRON_NEUTRINO: 1.0,
    ParticleType.MUON_NEUTRINO: 1.0,
    ParticleType.TAU_NEUTRINO: 1.0,

    ParticleType.PHOTON: 0.0,
    ParticleType.GLUON: 0.0,
    ParticleType.W_PLUS: 0.0,
    ParticleType.W_MINUS: 0.0,
    ParticleType.Z_BOSON: 0.0,
    ParticleType.HIGGS: 0.0,
    ParticleType.GRAVITON: 0.0,
}


QUARK_TYPES = {
    ParticleType.UP,
    ParticleType.DOWN,
    ParticleType.CHARM,
    ParticleType.STRANGE,
    ParticleType.TOP,
    ParticleType.BOTTOM,
}

LEPTON_TYPES = {
    ParticleType.ELECTRON,
    ParticleType.MUON,
    ParticleType.TAU,
    ParticleType.ELECTRON_NEUTRINO,
    ParticleType.MUON_NEUTRINO,
    ParticleType.TAU_NEUTRINO,
}

BOSON_TYPES = {
    ParticleType.PHOTON,
    ParticleType.GLUON,
    ParticleType.W_PLUS,
    ParticleType.W_MINUS,
    ParticleType.Z_BOSON,
    ParticleType.HIGGS,
    ParticleType.GRAVITON,
}


@dataclass
class Particle:
    ptype: ParticleType
    pos: Tuple[float, float]
    vel: Tuple[float, float]
    mass: float = field(init=False)
    charge: float = field(init=False)
    spin: float = field(init=False)
    bnum: float = field(init=False)
    lnum: float = field(init=False)
    color: Color = field(default=None)
    id: int = field(default_factory=lambda: random.randrange(10**9))

    def __post_init__(self) -> None:
        self.mass = PARTICLE_MASS[self.ptype]
        self.charge = PARTICLE_CHARGE[self.ptype]
        self.spin = PARTICLE_SPIN[self.ptype]
        self.bnum = PARTICLE_BNUMBER[self.ptype]
        self.lnum = PARTICLE_LNUMBER[self.ptype]
        if self.ptype in QUARK_TYPES and self.color is None:
            self.color = random.choice(["r", "g", "b"])  # naive color assignment
        if self.ptype not in QUARK_TYPES:
            self.color = None

    @property
    def momentum(self) -> Tuple[float, float]:
        return (self.mass * self.vel[0], self.mass * self.vel[1])

    def kinetic_energy(self) -> float:
        return 0.5 * self.mass * (self.vel[0] ** 2 + self.vel[1] ** 2)

    def copy(self) -> "Particle":
        p = Particle(self.ptype, self.pos, self.vel, color=self.color)
        p.id = self.id
        return p

    @staticmethod
    def random_particle(ptype: ParticleType, world_size: float) -> "Particle":
        x = random.uniform(0.0, world_size)
        y = random.uniform(0.0, world_size)
        vx = random.uniform(-1.0, 1.0)
        vy = random.uniform(-1.0, 1.0)
        return Particle(ptype, (x, y), (vx, vy))


def is_color_neutral(colors: Tuple[Color, ...]) -> bool:
    # Color neutral in this toy model: contains r,g,b exactly once each
    filtered = [c for c in colors if c in ("r", "g", "b")]
    return sorted(filtered) == ["b", "g", "r"]


def distance_periodic(a: Tuple[float, float], b: Tuple[float, float], size: float) -> float:
    dx = abs(a[0] - b[0])
    dy = abs(a[1] - b[1])
    if dx > size / 2:
        dx = size - dx
    if dy > size / 2:
        dy = size - dy
    return math.hypot(dx, dy)
