import re  # Regular expression module for reading reaction strings
from collections import OrderedDict
from re import escape
from time import time
import matplotlib.pyplot as p
import numpy as np
from scipy.integrate import solve_ivp, OdeSolution
from scipy.sparse import bsr_matrix
from scipy.optimize import minimize
from . import tools
import warnings


class Domain:
    """
    Basic functions of a modelling domain
    """

    def __init__(self):
        self.data = OrderedDict()  # keeps imported data
        self.solution = OrderedDict()  # keeps simulation results
        self.variables = (
            OrderedDict()
        )  # keeps variables(species) names and inlet functions
        self.initial_values = OrderedDict()  # Initial conditions as const or f(t)
        self.time_start = 0  # Start of simulation
        self.time_stop = 100  # End of simulation
        self.time_eval = None  # Time steps for simulation
        self._solver_params = None
        self.events = []

    def _ode(self, t, y):
        """
        Empty ode function, 
        y has to be in the order of self.variables
        """
        print("Empty ode function, y has to be in the order of self.variables")
        pass

    def import_data(self, file_name, plot=True):
        """
        Iporting .csv file with the format:

        time| variable1     | variable2     | ....
        ----------------------------------------
        0   | init.cond.1   | init.cond2.   |   ...
        ... |   ....        |   .....       |   ...

        """
        data = tools.get_data(file_name)
        self.data = data
        if self.data:
            self.time_start = self.data["t"][0]
            self.time_stop = self.data["t"][-1]
            self.time_eval = self.data["t"]
            for k in data:
                if k != "t":
                    self.initial_values[k] = data[k][0]

        if plot:
            for k in data:
                if k != "t":
                    p.plot(data["t"], data[k], "o", label=k)
            p.legend()
            p.show()

    def _sort_vars(self):
        """
        Sort internal variables to simplify indexing during integration
        """
        sorted_list = sorted(self.variables.keys())
        if type(self) is Chemistry:
            # sort stoichiometry
            for i, k in enumerate(self.variables.keys()):
                self.variables[k] = self.stoichiometry[:, i]
            # sorted stoichiometry
            self.stoichiometry = np.stack([self.variables[k] for k in sorted_list], 1)

            # sort orders
            for i, k in enumerate(self.variables.keys()):
                self.variables[k] = self.orders[:, i]
            # sorted orders
            self.orders = np.stack([self.variables[k] for k in sorted_list], 1)

        # sort initial values
        self.initial_values = OrderedDict(
            {k: self.initial_values[k] for k in sorted_list}
        )

        # sort species
        self.variables = OrderedDict({k: 0 for k in sorted_list})

    def solver_params(self, **kwargs):
        self._solver_params = kwargs
        return

    def run(self, plot=True, output=False):
        """ 
        Running simulation of a modelling domain
        """
        # TODO: pre-initialise solver separately
        # TODO: optimise overhead

        self.solution = OrderedDict.fromkeys(self.variables, [])
        self.solution["t"] = []
        self.dense = None
        t = self.time_start
        events = self.events.copy()

        iv = self.initial_values
        sol_initial_values = np.zeros(len(self.initial_values))
        for i, k in enumerate(iv):
            if callable(iv[k]):
                sol_initial_values[i] = iv[k](t)
            else:
                sol_initial_values[i] = iv[k]

        sol = None
        t0 = time()
        while True:
            params = {
                "t_eval": self.time_eval,
                "method": "BDF",
                "events": events,
                "dense_output": False,
            }
            if self._solver_params:
                params.update(self._solver_params)
            sol = solve_ivp(
                self._ode, [t, self.time_stop], sol_initial_values, **params
            )

            for i, k in enumerate(self.solution):
                if not k == "t":
                    self.solution[k] = np.append(self.solution[k], sol.y[i])
                    sol_initial_values[i] = sol.y[i][-1]
            self.solution["t"] = np.append(self.solution["t"], sol.t)

            # Dense output for solution interpolation
            if not self.dense:
                self.dense = sol.sol
            else:
                self.dense.interpolants.extend(sol.sol.interpolants)
                self.dense.t_max = sol.sol.t_max
                self.dense.n_segments += sol.sol.n_segments
                self.dense.ts = np.append(self.dense.ts, sol.sol.ts)
                self.dense.ts_sorted = np.append(
                    self.dense.ts_sorted, sol.sol.ts_sorted
                )

            # If aborted by an event: continue
            if sol.status == 1:
                for i, e in enumerate(sol.t_events):
                    if e:
                        events.pop(i)
                        t = e[0]
            # If aborted by t_stop: end
            elif sol.status == 0:
                break

        if type(self) is PFR:
            self.solution["t"][1:] += self.delay
        t1 = time()

        if plot:
            self.plot(all=True)
            print(f"run time: {t1-t0:.3f}s")

        if output:
            return sol

    def _obj_fun(self, parameters):
        # TODO: sensitivity matrix output
        self.rate_constants = parameters
        self.run(plot=False)
        obj = 0
        for key in self.data:
            if key != "t":
                obj += np.sum((self.solution[key] - self.data[key]) ** 2)
        return obj

    def fit(self, plot=True):
        """
        Fitting the model into imported data
        """
        # TODO: Change to Gauss-Newton method with sensitivity estimation
        # TODO: add checkboxes to choose components to optimise for
        res = minimize(
            self._obj_fun,
            self.rate_constants,
            # bounds = tuple((0,None) for i in self.rate_constants),
            method="Nelder-Mead",
        )
        print(res)
        # self.rate_constants = res
        self.run(plot)
        return res

    def plot(self, *args, all=False):
        """
        Plotting selected variables as a function of time
        """
        # TODO check if solution exists before plotting
        #      add checkboxes for components to plot
        if all:
            to_plot = self.variables
        else:
            to_plot = args

        fig = p.figure()
        ax = fig.add_subplot(111)
        ax.set_ylabel("concentration")
        ax.set_xlabel("time")
        for variable in to_plot:
            l = ax.plot(self.solution["t"], self.solution[variable], label=variable)
            if callable(self.initial_values[variable]):
                t = np.linspace(self.time_start, self.time_stop, 1000)
                l = ax.plot(
                    t,
                    [self.initial_values[variable](t) for t in t],
                    label="inlet " + variable,
                    # drawstyle = 'steps-post'
                )
            if variable in self.data:
                ax.plot(
                    self.data["t"],
                    self.data[variable],
                    "o",
                    label="exp." + variable,
                    color=l[-1].get_color(),
                )
        # ax.set_xlim(self.time_start,self.time_stop)
        ax.legend()
        p.show()
        # todo: return axis object
        return ax


class Chemistry(Domain):
    def __init__(self):
        super().__init__()
        self.stoichiometry = np.array([], ndmin=2)
        self.rate_constants = []
        self.orders = np.array([], ndmin=2)
        self.parameters = self.rate_constants

    def _rate(self, c):
        """
        Returns array of species generation rates.

        Warning: if rate fractional orders are implemented, np.power can't raise negative number to fractional exponent,
            see here: https://stackoverflow.com/questions/34898917/how-to-raise-arrays-with-negative-values-to-fractional-power-in-python
        """
        # ---Warning--- flooring c to 0
        # c[c < 0] = 0
        r = c ** self.orders
        r = np.prod(r, 1) * self.rate_constants

        return r

    def _ode(self, t, c):
        """
        ODE to be solved:
        only species generation term - reaction rate
        """
        dcdt = np.dot(self._rate(c), self.stoichiometry)
        return dcdt

    # Jacobian method: curenntly doesn't give more performance.
    # def _jac(self, t, c):
    #     Jo = self.orders.copy()
    #     Jo[Jo!=0]-=1
    #     J = np.dot(self.stoichiometry.transpose(), self.orders+Jo*c)
    #     return J

    def initial_concentrations(self, **kwargs):
        """
        Set the initial concentrations of species at t=0s
        """
        for species, conc in kwargs.items():
            if species in self.variables:
                self.initial_values.update({species: conc})

    def reaction(self, string, k=1, k1=1, k2=1):
        """
            Adding new chemical reaction to the system of reactions 
            in the current modelling domain.
            
            Parameters
            ----------
            string : str
                The chemical reaction in string form
            k : float
                rate constant in case of irreversible reaction
            k1 : float
                forward rate constant for reversible reaction
            k2 : float
                backward rate constant for reversible reaction
            
            Examples
            --------
            chem.reaction('A+B=>C')
            chem.reaction('A+B<=>C')
            chem.reaction('A+B=>C', k = 1e-6)
            chem.reaction('A+B<=>C', k1 = 3.5e-12, k2 = 4.5e-11)
        """
        for i in self.variables:
            self.variables[i] = 0

        string = string.replace(" ", "")

        irreversible = r"=>"
        reversible = r"<=>"

        if re.search(reversible, string):
            [reagents, products] = re.split(reversible, string)
            self._new_reaction(reagents, products, k1)
            self._new_reaction(products, reagents, k2)
        elif re.search(irreversible, string):
            [reagents, products] = re.split(irreversible, string)
            self._new_reaction(reagents, products, k)
        else:
            print("reactions not found")

        for s in self.variables:
            if s in self.initial_values:
                self.initial_values.update({s: self.initial_values[s]})
            else:
                self.initial_values.update({s: 0})
        # sorting variables for easier integration into other domains
        self._sort_vars()
        return

    def _new_reaction(self, reagents, products, k):
        """
        This could be optimised, maybe?
        dim0: all species of a single reaction
        dim1: all reactions of a single species
        """
        # TODO: rewrite from self.species to new dictionaries

        self.rate_constants.append(k)

        def get_coeff(s):
            match = re.search(r"\b[\d\.]+", s)
            if match:
                return float(match.group())
            else:
                return 1

        def get_species(s):
            return re.findall(r"[A-z]\w*", s)[0]

        # ---Creating stoichiometry matrix
        # -Using species dict as temporary storage for coefficients

        for i in re.findall(r"[\.\w]+", reagents):
            species = get_species(i)
            coeff = get_coeff(i)
            if species in self.variables:
                self.variables[species] += -coeff
            else:
                self.variables[species] = -coeff

        for i in re.findall(r"[\.\w]+", products):
            species = get_species(i)
            coeff = get_coeff(i)
            if species in self.variables:
                self.variables[species] += coeff
            else:
                self.variables.update({species: coeff})

        new_s = [self.variables.get(i) for i in self.variables]
        s = self.stoichiometry

        if s.size < 1:
            self.stoichiometry = np.array(new_s, ndmin=2)
        else:
            c = np.zeros((s.shape[0], len(new_s)))
            c[: s.shape[0], : s.shape[1]] = s
            c = np.vstack((c, new_s))
            self.stoichiometry = c

        # ---Creating species rate order matrix

        for i in self.variables:
            self.variables[i] = 0

        for i in re.findall(r"[\.\w]+", reagents):
            self.variables[get_species(i)] += abs(get_coeff(i))

        new_o = [self.variables.get(i) for i in self.variables]
        o = self.orders

        if o.size < 1:
            self.orders = np.array(new_o, ndmin=2)
        else:
            c = np.zeros((o.shape[0], len(new_o)))
            c[: o.shape[0], : o.shape[1]] = o
            c = np.vstack((c, new_o))
            self.orders = c

        if np.any(self.orders % 1 > 0):
            warnings.warn_explicit(
                message="Fractional rate orders are not supported and will be rounded to nearest integer",
                filename="models.py",
                lineno=406,
                category=Warning,
            )

        self.orders = self.orders.astype(int)

        # Cleaning species: temp storage for coeffitients
        for i in self.variables:
            self.variables[i] = 0

        return


class CSTR(Domain):
    def __init__(self, q=1, V=1):
        super().__init__()
        self.source = None
        self.destination = None
        self.q = q
        self.V = V
        self._chemistry = None
        self._chemistry_ind = False
        self._jac = None

    def _get_c_in(self, t):
        # used to get inlet variables values
        iv = self.initial_values
        c_in = np.zeros(len(iv))
        for i, k in enumerate(iv):
            if callable(iv[k]):
                c_in[i] = iv[k](t)
            else:
                c_in[i] = iv[k]
        return c_in

    def _ode(self, t, c):
        dmdt = self.q / self.V * (self._get_c_in(t) - c)
        if self._chemistry:
            dmdt[self._chemistry_ind] += self._chemistry._ode(t, c[self._chemistry_ind])
        return dmdt

    def inlet(self, **kwargs):
        for k, v in kwargs.items():
            if type(v) == tuple:
                self.initial_values.update({k: v[0]})
                self.events = v[1]
            elif callable(v):
                self.initial_values.update({k: v})
            else:
                self.initial_values.update({k: v})
            self.variables.update({k: 0})
        self._iv = [v for v in self.initial_values.values()]
        self._update_vars()
        self._sort_vars()

    @property
    def chemistry(self):
        return self._chemistry

    @chemistry.setter
    def chemistry(self, domain):
        """
        Adds chemistry domain into the reactor
        """
        if type(domain) != Chemistry:
            return

        self._chemistry = domain
        self._update_vars()
        self._sort_vars()
        return

    def _update_vars(self):
        if self._chemistry:
            ch_vars = set(self._chemistry.variables)
            self_vars = set(self.variables)
            diff = ch_vars - self_vars
            new_vars = dict().fromkeys(diff, 0)
            self.variables.update(new_vars)
            self.initial_values.update(new_vars)

            # Indexes of vars present in chemistry domain:
            self._chemistry_ind = [
                (v in self._chemistry.variables) for v in self.variables
            ]
        return


class PFR(Domain):
    def __init__(self, q=1, V=1):
        super().__init__()
        self.source = None
        self.destination = None
        self.q = q
        self.V = V
        self.delay = V / q
        self._store = []
        self._chemistry = None
        self._chemistry_ind = False
        self._jac = None

    def _get_c_in(self, t):
        # used to get inlet variables values
        iv = self.initial_values
        c_in = np.zeros(len(iv))
        for i, k in enumerate(iv):
            if callable(iv[k]):
                c_in[i] = iv[k](t)
            else:
                c_in[i] = iv[k]
        return c_in

    def _ode(self, t, c):
        c_in = self._get_c_in(t)
        dmdt = self.q / 0.0001 * (c_in - c)
        if self._chemistry:
            self._chemistry.initial_concentrations(
                **{
                    k: c_in[self._chemistry_ind][i]
                    for i, k in enumerate(self._chemistry.variables.keys())
                }
            )
            self._chemistry.time_stop = self.delay
            sol = self._chemistry.run(plot=False, output=True)
            dmdt[self._chemistry_ind] += (
                self.q / 0.0001 * ([v[-1] for v in sol.y] - c_in[self._chemistry_ind])
            )
        return dmdt

    def inlet(self, **kwargs):
        for k, v in kwargs.items():
            if type(v) == tuple:
                self.initial_values.update({k: v[0]})
                self.events = v[1]
            elif callable(v):
                self.initial_values.update({k: v})
            else:
                self.initial_values.update({k: v})
            self.variables.update({k: 0})
        self._iv = [v for v in self.initial_values.values()]
        self._update_vars()
        self._sort_vars()

    @property
    def chemistry(self):
        return self._chemistry

    @chemistry.setter
    def chemistry(self, domain):
        """
        Adds chemistry domain into the reactor
        """
        if type(domain) != Chemistry:
            return

        self._chemistry = domain
        self._update_vars()
        self._sort_vars()
        return

    def _update_vars(self):
        if self._chemistry:
            ch_vars = set(self._chemistry.variables)
            self_vars = set(self.variables)
            diff = ch_vars - self_vars
            new_vars = dict().fromkeys(diff, 0)
            self.variables.update(new_vars)
            self.initial_values.update(new_vars)

            # Indexes of vars present in chemistry domain:
            self._chemistry_ind = [
                (v in self._chemistry.variables) for v in self.variables
            ]
        return

