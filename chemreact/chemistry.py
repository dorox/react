from chemreact.models2 import Domain, Variable
from typing import List, Any, Union
import numpy as np
import re

array_like = Union[list, np.array, np.ndarray]
number = Union[int, float]

class Reaction:
    REVERSIBLE = r"<=>"
    IRREVERSIBLE = r"=>"
    EQUILIBRIUM = r"="

    def __init__(self, s:str , k:float) -> None:
        s = s.replace(" ", "")
        self._name = s
        self._reagents:list = []
        self._products:list = []
        self._k = k
        (left, right) = re.split(self.IRREVERSIBLE, s)
        self._read_side(left, reagents=True)
        self._read_side(right, reagents=False)

    def __repr__(self)->str:
        return f'{self._name} - {self.__class__}'

    def _read_side(self, s:str, reagents:bool)->None:
        
        def get_coeff(s:str) -> float:
            match = re.search(r"\b[\d\.]+", s)
            if match:
                return float(match.group())
            else:
                return 1
        
        def get_species(s:str):
            return re.findall(r"[A-z]\w*", s)[0]

        for i in re.findall(r"[\.\w]+", s):
            species = get_species(i)
            coeff = get_coeff(i)
            v = Variable(name=species, initial_value=0)
            if reagents:
                v.is_reagent = True
                v.coeff = -coeff
                v.order = coeff
                self._reagents.append(v)
            else:
                v.is_reagent = False
                v.coeff = coeff
                v.order = 0
                self._products.append(v)


class Chemistry(Domain):
    def __init__(self) -> None:
        super().__init__()
        self._stoichiometry : np.ndarray = np.array([], ndmin=2)
        self._rate_constants : np.ndarray = np.array([])
        self._orders : np.ndarray = np.array([], ndmin=2)
        self._reactions : List[Reaction] = []

    def ode(self, t:float, y:array_like) -> array_like:
        dcdt = np.dot(
            np.prod(y**self._orders, 1) * self._rate_constants, 
            self._stoichiometry
            )
        return dcdt

    def reaction(self, s:str, k:float = 1):
        r = Reaction(s, k)
        self._reactions.append(r)
        self._upd()

    def _upd(self):
        st = []
        orders = []
        rk = []
        for r in self._reactions:
            rk.append(r._k)
            for v in r._reagents + r._products:
                st.append(v.coeff)
                orders.append(v.order)
                self.add_variable(v.name, v=v)
        self._stoichiometry = np.array(st, ndmin=2)
        self._orders = np.array(orders, ndmin = 2)
        self._rate_constants = np.array(rk)
    