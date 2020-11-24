from chemreact.models2 import Constant, Domain, Variable
import numpy as np
import re


class Species(Variable):
    def __init__(self, order, stoichiometry, is_reagent, **kwargs):
        super().__init__(**kwargs)
        self.order = order
        self.stoichiometry = stoichiometry
        self.is_reagent = is_reagent


class Reaction(Domain):
    """A single reaction, doesn't know about chemistry its in"""

    REVERSIBLE = r"<=>"
    IRREVERSIBLE = r"=>"
    EQUILIBRIUM = r"="

    def __init__(self, s, k0=1, idx=None, **kwargs):
        s = s.replace(" ", "")
        species = self._read(s)

        super().__init__(
            name=s,
            variables=species,
            constants=[Constant("k0", k0, index=0)],
            index=idx,
            **kwargs,
        )
        self._stoichiometry = dict().setdefault(0)
        self._orders = dict().setdefault(0)

    def __repr__(self):
        return f"Reaction {self._name}"

    def _read(self, s, reac_type=IRREVERSIBLE):
        def get_coeff(s):
            match = re.search(r"\b[\d\.]+", s)
            if match:
                return float(match.group())
            else:
                return 1

        def get_species(s):
            return re.findall(r"[A-z]\w*", s)[0]

        def read_side(side, is_reagents, group):
            for i in re.findall(r"[\.\w]+", side):
                name = get_species(i)
                coeff = get_coeff(i)
                coeff *= -1 if is_reagents else 1
                order = abs(coeff) if is_reagents else 0
                species = Species(
                    name=name,
                    stoichiometry=coeff,
                    order=order,
                    is_reagent=is_reagents,
                    initial_value=0,
                )
            if name not in group:
                group[name] = species
            else:
                group[name].stoichiometry += coeff
                group[name].order += order
            return group

        left, right = re.split(reac_type, s)
        species = read_side(left, True, dict())
        species = read_side(right, False, species)

        return species_group.values()

    def ode(self, t, y):
        return np.dot(self.k0 * np.prod(y ** self.orders), self.stoichiometry)

    @property
    def species(self):
        return self._vars

    @property
    def stoichiometry(self):
        st = []
        for s in self.species:
            st.append(s.stoichiometry)
        return st

    @property
    def orders(self):
        o = []
        for s in self.species:
            o.append(s.order)
        return o


class Chemistry(Domain):
    """Creates system of reactions
    solves ODE in a matrix form
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._stoichiometry = np.array([], ndmin=2)
        self._rate_constants = np.array([])
        self._orders = np.array([], ndmin=2, dtype=int)
        self._reactions = []

    def ode(self, t, y):
        return np.dot(
            np.prod(y ** self._orders, 1) * self._rate_constants, self._stoichiometry
        )

    def reaction(self, s, k=1):
        # Adds reaction
        r = Reaction(s, k)
        for sp in r.species:
            if sp not in self._vars:
                self.add_variable(sp.name, sp)
        self._reactions.append(r)
        self._stoichiometry, self._orders, self._rate_constants = self._upd()
        return r

    def initial_concentrations(self, **kwargs):
        for sp, ic in kwargs.items():
            if sp in self.__dict__:
                getattr(self, sp).initial_value = ic

    def _upd(self):
        # prepares to start simulation
        s = np.zeros((len(self._reactions), len(self._vars)))
        o = np.zeros((len(self._reactions), len(self._vars)), dtype=int)
        rk = np.zeros(len(self._reactions))
        for i, r in enumerate(self._reactions):
            s[i] = r.stoichiometry
            o[i] = r.orders
            rk[i] = r.k0.value
        self._upd_idx()
        self.y0  # updates indexe
        return s, o, rk
