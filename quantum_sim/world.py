from __future__ import annotations

import math
import random
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Set, Tuple

from .constants import (
    BIND_DISTANCE,
    BOND_DAMPING,
    BOND_K,
    DT,
    G_NEWTON,
    JITTER_STD,
    K_COULOMB,
    MAX_SPEED,
    SOFTENING,
    STRONG_ATTRACTION,
    STRONG_RANGE,
    WORLD_SIZE,
)
from .particles import (
    BOSON_TYPES,
    Particle,
    ParticleType,
    QUARK_TYPES,
    distance_periodic,
)


@dataclass(frozen=True)
class Bond:
    a: int  # particle id
    b: int  # particle id
    rest: float = BIND_DISTANCE


class Universe:
    def __init__(
        self,
        size: float = WORLD_SIZE,
        seed: Optional[int] = None,
        counts: Optional[Dict[ParticleType, int]] = None,
    ) -> None:
        self.size = size
        self.rng = random.Random(seed)
        self.particles: List[Particle] = []
        self.bonds: Set[Bond] = set()
        if counts:
            self.spawn(counts)

    # Initialization
    def spawn(self, counts: Dict[ParticleType, int]) -> None:
        for ptype, n in counts.items():
            for _ in range(max(0, n)):
                p = Particle.random_particle(ptype, self.size)
                # use our seeded rng for determinism
                x = self.rng.uniform(0.0, self.size)
                y = self.rng.uniform(0.0, self.size)
                vx = self.rng.uniform(-1.0, 1.0)
                vy = self.rng.uniform(-1.0, 1.0)
                p.pos = (x, y)
                p.vel = (vx, vy)
                self.particles.append(p)

    # Core loop
    def step(self, steps: int = 1) -> None:
        for _ in range(steps):
            self._detect_new_baryon_bonds()
            acc = self._compute_accelerations()
            self._integrate(acc)
            self._apply_periodic_boundaries()
            self._maybe_break_bonds()

    # Forces & interactions
    def _compute_accelerations(self) -> List[Tuple[float, float]]:
        n = len(self.particles)
        ax = [0.0] * n
        ay = [0.0] * n
        # Pair interactions
        for i in range(n):
            pi = self.particles[i]
            for j in range(i + 1, n):
                pj = self.particles[j]
                dx, dy, r, inv_r3 = self._pair_geometry(pi.pos, pj.pos)
                if r < 1e-6:
                    continue
                # Electromagnetic (Coulomb)
                qprod = pi.charge * pj.charge
                if qprod != 0.0:
                    f = K_COULOMB * qprod * inv_r3
                    ax[i] += f * dx / pi.mass
                    ay[i] += f * dy / pi.mass
                    ax[j] -= f * dx / pj.mass
                    ay[j] -= f * dy / pj.mass
                # Gravity (Newtonian)
                g = G_NEWTON * pi.mass * pj.mass * inv_r3
                ax[i] -= g * dx / pi.mass
                ay[i] -= g * dy / pi.mass
                ax[j] += g * dx / pj.mass
                ay[j] += g * dy / pj.mass
                # Strong force (toy): short-range attraction among quarks
                if (pi.ptype in QUARK_TYPES) and (pj.ptype in QUARK_TYPES):
                    if r < STRONG_RANGE:
                        # slight attraction trying to reduce separation
                        s = STRONG_ATTRACTION * (1.0 - r / STRONG_RANGE)
                        ax[i] -= s * dx / (r * pi.mass)
                        ay[i] -= s * dy / (r * pi.mass)
                        ax[j] += s * dx / (r * pj.mass)
                        ay[j] += s * dy / (r * pj.mass)
        # Bonded springs
        id_to_index = {p.id: k for k, p in enumerate(self.particles)}
        for bond in self.bonds:
            if bond.a not in id_to_index or bond.b not in id_to_index:
                continue
            ia = id_to_index[bond.a]
            ib = id_to_index[bond.b]
            pa = self.particles[ia]
            pb = self.particles[ib]
            dx, dy, r, _ = self._pair_geometry(pa.pos, pb.pos)
            if r < 1e-6:
                continue
            # spring force: F = -k (r - r0) along direction
            extension = r - bond.rest
            fx = -BOND_K * extension * (dx / r)
            fy = -BOND_K * extension * (dy / r)
            # damping proportional to relative velocity
            rvx = pa.vel[0] - pb.vel[0]
            rvy = pa.vel[1] - pb.vel[1]
            fx -= BOND_DAMPING * rvx
            fy -= BOND_DAMPING * rvy
            ax[ia] += fx / pa.mass
            ay[ia] += fy / pa.mass
            ax[ib] -= fx / pb.mass
            ay[ib] -= fy / pb.mass

        # stochastic jitter
        for i in range(n):
            ax[i] += self.rng.gauss(0.0, JITTER_STD)
            ay[i] += self.rng.gauss(0.0, JITTER_STD)
        return list(zip(ax, ay))

    def _pair_geometry(self, a: Tuple[float, float], b: Tuple[float, float]) -> Tuple[float, float, float, float]:
        # minimal image under periodic boundary
        dx = b[0] - a[0]
        dy = b[1] - a[1]
        if dx > self.size / 2:
            dx -= self.size
        elif dx < -self.size / 2:
            dx += self.size
        if dy > self.size / 2:
            dy -= self.size
        elif dy < -self.size / 2:
            dy += self.size
        r2 = dx * dx + dy * dy + SOFTENING * SOFTENING
        r = math.sqrt(r2)
        inv_r3 = 1.0 / (r2 * r)
        return dx, dy, r, inv_r3

    def _integrate(self, acc: List[Tuple[float, float]]) -> None:
        for i, p in enumerate(self.particles):
            vx = p.vel[0] + acc[i][0] * DT
            vy = p.vel[1] + acc[i][1] * DT
            # speed cap for stability
            speed = math.hypot(vx, vy)
            if speed > MAX_SPEED:
                scale = MAX_SPEED / max(speed, 1e-9)
                vx *= scale
                vy *= scale
            x = (p.pos[0] + vx * DT) % self.size
            y = (p.pos[1] + vy * DT) % self.size
            p.vel = (vx, vy)
            p.pos = (x, y)

    def _apply_periodic_boundaries(self) -> None:
        # positions are wrapped in _integrate; nothing else to do here
        return

    def _detect_new_baryon_bonds(self) -> None:
        # naive O(N^3): search for close triads r,g,b quarks and bond them if not yet bonded
        quark_indices = [i for i, p in enumerate(self.particles) if p.ptype in QUARK_TYPES]
        n = len(quark_indices)
        if n < 3:
            return
        # map particle id to index for bond existence checks
        id_to_index = {p.id: i for i, p in enumerate(self.particles)}
        existing_pairs = {(min(b.a, b.b), max(b.a, b.b)) for b in self.bonds}
        for a_idx in range(n):
            ia = quark_indices[a_idx]
            pa = self.particles[ia]
            for b_idx in range(a_idx + 1, n):
                ib = quark_indices[b_idx]
                pb = self.particles[ib]
                if pa.color == pb.color:
                    continue
                if distance_periodic(pa.pos, pb.pos, self.size) > STRONG_RANGE:
                    continue
                for c_idx in range(b_idx + 1, n):
                    ic = quark_indices[c_idx]
                    pc = self.particles[ic]
                    colors = (pa.color, pb.color, pc.color)
                    if sorted(colors) != ["b", "g", "r"]:
                        continue
                    if distance_periodic(pa.pos, pc.pos, self.size) > STRONG_RANGE:
                        continue
                    if distance_periodic(pb.pos, pc.pos, self.size) > STRONG_RANGE:
                        continue
                    # create bonds among the triad if not already present
                    pairs = [
                        (min(pa.id, pb.id), max(pa.id, pb.id)),
                        (min(pa.id, pc.id), max(pa.id, pc.id)),
                        (min(pb.id, pc.id), max(pb.id, pc.id)),
                    ]
                    new_bonds = []
                    for a, b in pairs:
                        if (a, b) not in existing_pairs:
                            new_bonds.append(Bond(a, b))
                    if new_bonds:
                        self.bonds.update(new_bonds)
                        existing_pairs.update((min(b.a, b.b), max(b.a, b.b)) for b in new_bonds)
                    # only one triad per outer loops to limit O(N^3) cost
                    break

    def _maybe_break_bonds(self) -> None:
        # Break bonds that are stretched too far
        to_remove: List[Bond] = []
        id_to_index = {p.id: i for i, p in enumerate(self.particles)}
        for b in self.bonds:
            if b.a not in id_to_index or b.b not in id_to_index:
                to_remove.append(b)
                continue
            ia = id_to_index[b.a]
            ib = id_to_index[b.b]
            pa = self.particles[ia]
            pb = self.particles[ib]
            d = distance_periodic(pa.pos, pb.pos, self.size)
            if d > 3.0 * STRONG_RANGE:
                to_remove.append(b)
        for b in to_remove:
            self.bonds.discard(b)

    # Utilities
    def total_energy(self) -> float:
        return sum(p.kinetic_energy() for p in self.particles)

    def counts_by_type(self) -> Dict[ParticleType, int]:
        counts: Dict[ParticleType, int] = {}
        for p in self.particles:
            counts[p.ptype] = counts.get(p.ptype, 0) + 1
        return counts

    def snapshot(self, limit: Optional[int] = None) -> List[Tuple[int, str, float, float, float, float, Optional[str]]]:
        items = []
        for i, p in enumerate(self.particles):
            if limit is not None and i >= limit:
                break
            items.append((
                p.id, p.ptype.value, p.pos[0], p.pos[1], p.vel[0], p.vel[1], p.color
            ))
        return items
