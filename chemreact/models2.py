from typing import Any, Callable, List, Tuple, Union
import numpy as np
from scipy.integrate import solve_ivp


class Variable:
    def __init__(
        self,
        name: str,
        initial_value: float,
        index: int = None,
        var_type: str = None,
        domain=None,
    ):
        self.initial_value = initial_value
        self.name = name
        self._idx = index
        self.var_type = var_type
        self._domain = domain

    def __get__(self, obj, klass):
        if isinstance(obj, Domain) and obj._is_running:
            return obj._y[self._idx]
        else:
            return self

    def __set__(self, obj, val):
        self.initial_value = float(val)

    def __repr__(self) -> str:
        return f"Variable {self.name}, iv={self.initial_value}, idx={self._idx}"


class Constant:
    def __init__(self, name: str, value: float, index: int = None):
        self.name = name
        self.value = value
        self._idx = index

    def __get__(self, obj, klass):
        if isinstance(obj, Domain) and obj._is_running:
            return self.value
        else:
            return self

    def __set__(self, obj, val):
        if isinstance(obj, Domain):
            self.initial_value = float(val)
        else:
            NotImplemented

    def __repr__(self) -> str:
        return f"Constant {self.name}={self}"


class Domain:
    name: str = ""
    _is_running: bool = False
    _y: list = []

    def __new__(cls, name: str):
        """creates new subclass for descriptor attachment"""
        cls = type(f"{cls.__name__}_{name}", (cls,), {})
        cls.name = name
        return super().__new__(cls)

    def __init__(
        self,
        name: str,
        variables: List[Variable] = [],
        constants: List[Constant] = [],
        subdomains: list = [],
        index: int = None,
    ) -> None:
        s = self.__class__
        s._vars: List[Variable] = []
        s._const: List[Constant] = []
        s._subdomains: list = []
        s._idx: int = index
        for v in variables:
            self.add_variable(v.name, v)
        for c in constants:
            self.add_constant(c.name, c.value, c)
        for d in subdomains:
            self.add_subdomain(d.name, d)

    @classmethod
    def add_variable(cls, n: Union[str, Variable], v: Variable = None, **kwargs):
        if type(n) == str and v == None:
            v = Variable(name=n, **kwargs)
        elif type(n) == str and v == Variable:
            v.name = n
        elif isinstance(n, Variable):
            v = n
        if n in cls.__dict__:
            raise ValueError(f"{v} is already in {cls}")
        cls._vars.append(v)
        type.__setattr__(cls, v.name, v)

    @classmethod
    def add_constant(cls, name, value: float, c: Constant = None, **kwargs):
        if c == None:
            c = Constant(name=name, value=value, **kwargs)
        if c in cls.__dict__:
            raise ValueError
        cls._const.append(c)
        type.__setattr__(cls, c.name, c)

    @classmethod
    def add_subdomain(cls, name: str, d):
        d.name = name
        cls._subdomains.append(d)
        type.__setattr__(cls, d.name, d)

    def _ode(self, t: float, y: list) -> Callable:
        Domain._y = y
        return self.ode(t, y)

    def _upd_idx(self, idx0: int = 0):
        for idx, v in enumerate(self._vars):
            v._idx = idx + idx0
        for d in self._subdomains:
            idx0 = d._upd_idx(idx + idx0 + 1)
        return idx0

    @property
    def y0(self) -> list:
        y0 = []
        for v in self._vars:
            y0.append(v.initial_value)
        for d in self._subdomains:
            y0.extend(d.y0)
        return y0

    def run(self, *args, **kwargs):
        t_span = (0, 100)
        Domain._is_running = True
        sol = solve_ivp(self._ode, t_span, self.y0)
        Domain._is_running = False
        return sol

    def __repr__(self):
        return f"Domain {self.name}"

