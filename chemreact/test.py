import re
import numpy as np
from scipy.integrate import solve_ivp


def ode(t, y):
    dcdt = -y[0]
    return dcdt


y0 = [1]
t_span = [0, 1]
sol = solve_ivp(ode, t_span, y0)
print(sol)


class Chemistry:
    def __init__(self, *args, **kwargs):
        """
        Creating a new modelling domain
        """
        self.reactions = dict()
        if args:
            self.add_reaction(self, *args, **kwargs)

    def add_reaction(self, string: str, k=1, kf=1, kr=1):
        """
        Adding a new chemical reaction to the system of reactions 
        in the current modelling domain.
        """

        self.reactions[len(self.reactions)] = Reaction(string, k, kf, kr)

        pass

    pass


class Reaction:
    def __init__(self, string: str, k=1, kf=1, kr=1):
        self.reversible = None

        string = string.replace(" ", "")

        if re.search(r"=>", string):
            # Irreversible
            self.reversible = False
            reagents, products = re.split(r"=>", string)
            self.forward = self._elementary_reaction(reagents, products, k)
        elif re.search(r"<=>", string):
            # Reversible
            self.reversible = True
            reagents, products = re.split(r"<=>", string)
            self.forward = self._elementary_reaction(reagents, products, kf)
            self.reverse = self._elementary_reaction(products, reagents, kr)
        else:
            raise ValueError(f"No reactions found in the input: {string}")

    def _elementary_reaction(self, reagents, products, k):
        # Adds an elementary reaction to the reaction group

        def get_coeff(s):
            match = re.search(r"\b[\d\.]+", s)
            if match:
                return float(match.group())
            else:
                return 1

        def get_species(s):
            return re.findall(r"[A-z]\w*", s)[0]

        for i in re.findall(r"[\.\w]+", reagents):
            self.species[get_species(i)] = -get_coeff(i)

        for i in re.findall(r"[\.\w]+", products):
            product = get_species(i)
            coeff = get_coeff(i)
            self.species.setdefault(product, 0)
            self.species[product] += coeff

