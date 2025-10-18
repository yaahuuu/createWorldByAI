from __future__ import annotations

import math
from typing import Dict, List, Set, Tuple

from .constants import (
    CLUSTER_NEIGHBOR_RADIUS,
    COARSE_GRID_SIZE,
)
from .world import Universe
from .particles import distance_periodic


class UnionFind:
    def __init__(self, n: int) -> None:
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: int, b: int) -> None:
        ra = self.find(a)
        rb = self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            self.parent[ra] = rb
        elif self.rank[ra] > self.rank[rb]:
            self.parent[rb] = ra
        else:
            self.parent[rb] = ra
            self.rank[ra] += 1


def compute_clusters(u: Universe) -> List[List[int]]:
    n = len(u.particles)
    uf = UnionFind(n)
    # neighbor linking
    for i in range(n):
        pi = u.particles[i]
        for j in range(i + 1, n):
            pj = u.particles[j]
            d = distance_periodic(pi.pos, pj.pos, u.size)
            if d <= CLUSTER_NEIGHBOR_RADIUS:
                uf.union(i, j)
    # collect clusters
    roots: Dict[int, List[int]] = {}
    for i in range(n):
        r = uf.find(i)
        roots.setdefault(r, []).append(i)
    return list(roots.values())


def coarse_grain_density(u: Universe, grid: int = COARSE_GRID_SIZE) -> List[List[int]]:
    cell = u.size / grid
    field = [[0 for _ in range(grid)] for _ in range(grid)]
    for p in u.particles:
        ix = int(p.pos[0] // cell) % grid
        iy = int(p.pos[1] // cell) % grid
        field[iy][ix] += 1
    return field


def macro_metrics(u: Universe) -> Dict[str, float]:
    clusters = compute_clusters(u)
    largest = max((len(c) for c in clusters), default=0)
    macro_clusters = sum(1 for c in clusters if len(c) >= 30)
    # temperature proxy: average kinetic energy per particle
    n = len(u.particles)
    avg_ke = 0.0
    if n:
        avg_ke = sum(p.kinetic_energy() for p in u.particles) / n
    return {
        "num_particles": float(n),
        "num_clusters": float(len(clusters)),
        "largest_cluster": float(largest),
        "macro_clusters": float(macro_clusters),
        "avg_ke": avg_ke,
    }
