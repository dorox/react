from typing import Any, Callable, List, Tuple
import numpy as np
from scipy.integrate import solve_ivp, OdeSolution

class Variable:
    def __init__(self, name:str=None, initial_value:float=None, index:int=None, var_type:str=None):
        self.initial_value = initial_value
        self.name = name
        self._idx = index
        self.var_type = var_type

    def __repr__(self) -> str:
        return f'Variable {self.name}, iv={self.initial_value}, idx={self._idx}'
    
class Constant:
    def __init__(self, name:str=None, value:float=None, index:int=None):
        self.name = name
        self.value = value
        self._idx = index
    def __repr__(self) -> str:
        return f'Constant {self.name}={self}'

class Domain:
    '''Simulation domain:
    repersents a 1D object '''
    _y : list = []
    _is_running : bool = False
    def __init__(self) -> None:
        self._vars : List[Variable]= []
        self._const : List[Constant] = []
        self.subdomains : List[Domain]= []

    def __getattribute__(self, attr:str) -> Any:
        if not Domain._is_running:
            return object.__getattribute__(self, attr)
        d = object.__getattribute__(self, '__dict__')
        if attr in d and type(d[attr]) is Variable:
            return Domain._y[d[attr]._idx]
        else:
            return object.__getattribute__(self, attr)

    def add_variable(self, name, v:Variable=None, **kwargs):
        if v==None:
            v = Variable(name = name, **kwargs)
        self._vars.append(v)
        object.__setattr__(self, v.name, v)

    def add_constant(self, name:str, value:float, c:Constant=None, **kwargs):
        if c==None:
            c = Constant(name=name, value=value, **kwargs)
        self._const.append(c)
        object.__setattr__(self, c.name, c)

    def _ode(self, t:float, y:list) -> Callable:
        Domain._y = y
        return self.ode(t,y)

    def _get_y0(self, y0:list=[], idx0:int=0) -> Tuple[list, int]:
        for idx, v in enumerate(self._vars):
            v._idx = idx+idx0
            y0.append(v.initial_value)
        return y0, idx+1

    @property
    def y0(self) -> list :
        y0, idx0 = self._get_y0(y0=[],idx0=0)
        if len(self.subdomains) != 0:
            for d in self.subdomains:
                y0, idx0 = d._get_y0(y0, idx0)
        return y0

    def run(self) -> OdeSolution:
        t_span = (0, 100)
        Domain._is_running = True
        sol = solve_ivp(self._ode, t_span, self.y0)
        Domain._is_running = False
        return sol