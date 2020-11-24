from scipy.integrate import solve_ivp


class Variable:
    def __init__(
        self,
        name,
        initial_value=0,
    ):
        self.initial_value = initial_value
        self.name = name
        self._idx = 0
        self._ode = lambda: 0

    def __get__(self, obj, cls):
        if isinstance(obj, Domain) and obj._is_running:
            return obj._y[self._idx]
        else:
            return self

    @property
    def dt(self):
        return self._ode()

    def __set__(self, obj, val):
        self.initial_value = float(val)

    def __repr__(self):
        return f"Variable {self.name}, iv={self.initial_value}, idx={self._idx}"


class Constant:
    def __init__(self, name, value=0):
        self.name = name
        self.value = value
        self._idx = 0

    def __get__(self, obj, cls):
        if isinstance(obj, Domain) and obj._is_running:
            return self.value
        else:
            return self

    def __set__(self, obj, val):
        if isinstance(obj, Domain):
            self.initial_value = float(val)
        else:
            NotImplemented

    def __repr__(self):
        return f"Constant {self.name}={self}"


class Domain:
    _is_running = False
    _y = []

    def __new__(cls, name):
        """creates a new subclass for descriptor attachment"""
        cls = type(f"{cls.__name__}_{name}", (Domain,), {})
        cls.name = name
        return super().__new__(cls)

    def __init__(
        self,
        name,
        variables=[],
        constants=[],
        subdomains=[],
        index=None,
    ):
        self._vars = []
        self._const = []
        self._subdomains = []
        self._idx = index
        for v in variables:
            self.add_variable(v.name, v)
        for c in constants:
            self.add_constant(c.name, c.value, c)
        for d in subdomains:
            self.add_subdomain(d.name, d)

    def __setattr__(self, name, value):
        if isinstance(value, Variable):
            self._vars.append(value)
        elif isinstance(value, Constant):
            self._const.append(value)
        else:
            object.__setattr__(self, name, value)
            return
        type.__setattr__(self.__class__, name, value)

    def add_variable(self, v):
        if not isinstance(v, Variable):
            raise TypeError(f"{v} is not a Variable")
        elif v in self._vars:
            raise ValueError(f"{v} is already in {self}")
        elif v.name in self.__dir__():
            raise ValueError(f"Name {v.name} is occupied")
        self.__setattr__(v.name, v)

    def new_variables(self, s):
        if type(s) == str:
            self.add_variable(Variable(s))
        else:
            for v in s:
                self.add_variable(Variable(v))

    def add_constant(self, c):
        if not isinstance(c, Constant):
            raise TypeError(f"{c} is not a Constant")
        elif c in self._const:
            raise ValueError(f"{c} is already in {self}")
        elif c.name in self.__dir__():
            raise ValueError(f"Name {c.name} is occupied")
        self.__setattr__(c.name, c)

    def new_constants(self, s):
        if type(s) == str:
            self.add_constant(Constant(s))
        else:
            for c in s:
                self.add_constant(Constant(c))

    def add_subdomain(self, d):
        self._subdomains.append(d)
        type.__setattr__(self.__class__, d.name, d)

    def _ode(self, t, y):
        Domain._y = y
        return self.ode(t, y)

    def ode(self, t, y):
        # placeholder for user defined ODE function
        pass

    def _upd_idx(self, idx0=0):
        for idx, v in enumerate(self._vars):
            v._idx = idx + idx0
        for d in self._subdomains:
            idx0 = d._upd_idx(idx + idx0 + 1)
        return idx0

    @property
    def y0(self):
        self._upd_idx()
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
