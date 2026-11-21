"""Fixed-volume compartment graphs and a grouped turnover-chain sampler."""

from collections.abc import Mapping
import json
from math import isfinite
from operator import index
import re

from . import _native
from ._common import model_bytes
from .builder import mass_action
from .rules import _crn_error


class SpatialModel:
    """Build an ordinary reaction network from named well-mixed compartments.

    The graph uses directed per-molecule transport hazards. Compartment
    volumes are fixed. Rate templates can interpolate local species names
    and the fixed volume.
    """

    def __init__(self, name, volumes):
        from . import Model

        self.model = Model(name)
        self._volumes = dict(volumes)
        if not self._volumes or any(
            not isinstance(region, str)
            or not re.fullmatch(r"[A-Za-z0-9_]+", region)
            or not isfinite(volume)
            or volume <= 0
            for region, volume in self._volumes.items()
        ):
            raise ValueError(
                "compartments require identifiers and positive finite volumes"
            )
        self._species = {}

    def add_species(self, name, initial=0, *, reservoirs=()):
        """Add a species in every compartment, with optional fixed counts."""
        from collections.abc import Mapping

        fixed = set(reservoirs)
        if name in self._species or fixed - self._volumes.keys():
            raise ValueError("duplicate species or unknown reservoir")
        if isinstance(initial, Mapping) and initial.keys() - self._volumes.keys():
            raise ValueError("initial counts name an unknown compartment")
        regions = {}
        for region in self._volumes:
            count = index(
                initial.get(region, 0) if isinstance(initial, Mapping) else initial
            )
            label = f"{name}_{region}"
            self.model.add_species(label, count)
            regions[region] = (label, count, region in fixed)
        self._species[name] = regions
        return self

    def label(self, species, region):
        """Return the species label used by the lowered model."""
        return self._species[species][region][0]

    def local(
        self, name, rate, reactants, products, *, combinatorial=False, regions=None
    ):
        """Instantiate a local reaction template in the selected regions."""
        reactants, products = tuple(reactants), tuple(products)
        for region in self._volumes if regions is None else regions:
            volume = self._volumes[region]
            inputs = [self.label(s, region) for s in reactants]
            outputs = [self.label(s, region) for s in products]
            if any(self._species[s][region][2] for s in (*reactants, *products)):
                raise ValueError("local reactions may not change reservoir species")
            if isinstance(rate, str):
                symbols = {s: item[region][0] for s, item in self._species.items()}
                expr = rate.format_map({**symbols, "volume": volume})
                value = mass_action(expr) if combinatorial else expr
            else:
                value = rate
            self.model.add_reaction(value, inputs, outputs, name=f"{name}_{region}")
        return self

    def transport(self, species, origin, target, rate):
        """Add a directed, per-copy first-order hop."""
        if origin == target or not isfinite(rate) or rate < 0:
            raise ValueError(
                "transport needs distinct regions and a finite nonnegative rate"
            )
        a, copies, a_fixed = self._species[species][origin]
        b, _, b_fixed = self._species[species][target]
        if not (a_fixed and b_fixed):
            self.model.add_reaction(
                float(rate) * copies if a_fixed else float(rate),
                [] if a_fixed else [a],
                [] if b_fixed else [b],
                name=f"move_{a}_to_{target}",
            )
        return self

    def diffusion(self, species, first, second, diffusivity, spacing):
        """Add both directions on a uniform lattice at rate D / spacing²."""
        if not isfinite(spacing) or spacing <= 0:
            raise ValueError("lattice spacing must be positive and finite")
        if self._volumes[first] != self._volumes[second]:
            raise ValueError("unequal volumes require explicit transport rates")
        rate = diffusivity / spacing**2
        return self.transport(species, first, second, rate).transport(
            species, second, first, rate
        )

    def profile(self, values, species):
        """Select region columns from a trajectory/report array."""
        return values.sel(variable=[self.label(species, i) for i in self._volumes])


def spatial_turnover(
    initial,
    *,
    end,
    target,
    feedback,
    mutant_weight,
    death,
    hop,
    seed=0,
    mutation=0.0,
    fixed_left=False,
    labels=False,
    max_events=3_000_000_000,
):
    """Sample the endpoint of a two-genotype nearest-neighbor turnover chain.

    Each copy has birth rate max(0, death + feedback * (target - W - weight*M)).
    Death and migration have constant per-copy rates. New mutations are
    marked at birth, and labels have no effect on subsequent hazards.
    """
    if not isinstance(fixed_left, bool) or not isinstance(labels, bool):
        raise TypeError("fixed_left and labels are booleans")
    if isinstance(initial, Mapping):
        raise TypeError("initial must be an ordered sequence of count pairs")
    counts = [(index(a), index(b)) for a, b in initial]
    spec = {
        "target": float(target),
        "feedback": float(feedback),
        "mutant_weight": float(mutant_weight),
        "death": float(death),
        "hop": float(hop),
        "mutation": float(mutation),
        "fixed_left": fixed_left,
    }
    try:
        return json.loads(
            _native.turnover_json(
                model_bytes(spec),
                counts,
                float(end),
                index(seed),
                index(max_events),
                labels,
            )
        )
    except _native.NativeError as error:
        raise _crn_error(error) from None
