import re  # Regular expression module for reading reaction strings
from collections import OrderedDict
import matplotlib.pyplot as p
import numpy as np
from scipy.optimize import minimize
from .domains import Domain0d, Domain1d

class Chemistry(Domain0d):
    def __init__(self):
        super().__init__()
        self.rate_constants = []
        self.parameters = self.rate_constants

    def _rate(self, c):
        '''
        Returns array of species generation rates
        '''
        r = np.zeros(len(self.rate_constants))
        for i,reaction in enumerate(self.orders):
            r[i]=self.rate_constants[i]

            # np.prod is slower than a loop:
            # r[i]=np.prod(np.power(c,reaction))

            for j,coeff in enumerate(reaction):
                r[i]*=c[j]**(coeff)
        return r

    def _ode(self, t, c):
        '''
        ODE to be solved:
        only species generation term - reaction rate
        '''
        dcdt = np.dot(self._rate(c), self.stoichiometry)
        return dcdt

    # Jacobian method: curenntly doesn't give more performance.
    # def _jac(self, t, c):
    #     Jo = self.orders.copy()
    #     Jo[Jo!=0]-=1
    #     J = np.dot(self.stoichiometry.transpose(), self.orders+Jo*c)
    #     return J
    
    def initial_concentrations(self, **kwargs):
        '''
        Set the initial concentrations of species at t=0s
        '''
        for species, conc in kwargs.items():
            if species in self.variables:
                self.initial_values.update({species:conc})

    def reaction(self, string, k=1, k1=1, k2=1):
        '''
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
        '''
        for i in self.variables:
            self.variables[i] = 0

        string = string.replace(' ', '')

        irreversible = r'=>'
        reversible = r'<=>'
               
        if re.search(reversible, string):
            [reagents, products] = re.split(reversible, string)
            self._new_reaction(reagents, products, k1)
            self._new_reaction(products, reagents, k2)
        elif re.search(irreversible, string):
            [reagents, products] = re.split(irreversible, string)
            self._new_reaction(reagents, products, k)
        else:
            print('reactions not found')
        
        for s in self.variables:
            if s in self.initial_values:
                self.initial_values.update({s:self.initial_values[s]})
            else:
                self.initial_values.update({s:0})
        #sorting variables for easier integration into other domains
        self._sort_vars()
        return

    def _new_reaction(self, reagents, products, k):
        '''
        This could be optimised, maybe?
        dim0: all species of a single reaction
        dim1: all reactions of a single species
        '''
        #TODO: rewrite from self.species to new dictionaries

        self.rate_constants.append(k)

        def get_coeff(s): 
            match = re.search(r'\b[\d\.]+', s)
            if match:
                return float(match.group())
            else:
                return 1

        def get_species(s): return re.findall(r'[A-z]\w*', s)[0]

        #---Creating stoichiometry matrix
        #-Using species dict as temporary storage for coefficients
        self.variables.update({get_species(i):-get_coeff(i) for i in re.findall(r'[\.\w]+', reagents)})

        for i in re.findall(r'[\.\w]+', products):
            species = get_species(i)
            coeff = get_coeff(i)
            if species in self.variables:
                self.variables[species] += coeff
            else:
                self.variables.update({species:coeff})

        new_s = [self.variables.get(i) for i in self.variables]
        s = self.stoichiometry

        if s.size<1:
            self.stoichiometry = np.array(new_s, ndmin = 2)
        else:
            c = np.zeros((s.shape[0], len(new_s)))
            c[:s.shape[0], :s.shape[1]] = s
            c = np.vstack((c, new_s))
            self.stoichiometry = c
    
        #---Creating species rate order matrix

        for i in self.variables:
            self.variables[i] = 0
        
        self.variables.update({get_species(i):abs(get_coeff(i)) for i in re.findall(r'[\.\w]+', reagents)})
        new_o = [self.variables.get(i) for i in self.variables]
        o = self.orders

        if o.size<1:
            self.orders = np.array(new_o, ndmin = 2)
        else:
            c = np.zeros((o.shape[0], len(new_o)))
            c[:o.shape[0], :o.shape[1]] = o
            c = np.vstack((c, new_o))
            self.orders = c

        #Cleaning species: temp storage for coeffitients
        for i in self.variables:
            self.variables[i] = 0

        return

class CSTR(Domain0d):
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
        #used to get inlet variables values
        iv = self.initial_values
        c_in = np.zeros(len(iv))
        for i, k in enumerate(iv):
            if callable(iv[k]):
                c_in[i] = iv[k](t)
            else:
                c_in[i] = iv[k]
        return c_in

    def _ode(self, t, c):
        dmdt = self.q/self.V*(self._get_c_in(t)-c)
        if self._chemistry:
            dmdt[self._chemistry_ind] += self._chemistry._ode(t,c[self._chemistry_ind])
        return dmdt
    
    def inlet(self, **kwargs):
        for k, v in kwargs.items():
            if type(v)==tuple:
                self.initial_values.update({k:v[0]})
                self.events = v[1]
            elif callable(v):
                self.initial_values.update({k:v})
            else:
                self.initial_values.update({k:v})
            self.variables.update({k:0})
        self._iv = [v for v in self.initial_values.values()]
        self._update_vars()
        self._sort_vars()
    
    @property
    def chemistry(self):
        return self._chemistry

    @chemistry.setter
    def chemistry(self, domain):
        '''
        Adds chemistry domain into the reactor
        '''
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
            diff = ch_vars-self_vars
            new_vars = dict().fromkeys(diff, 0)
            self.variables.update(new_vars)
            self.initial_values.update(new_vars)
        
            #Indexes of vars present in chemistry domain:
            self._chemistry_ind = [(v in self._chemistry.variables) for v in self.variables]
        return
    
class PFR(Domain1d):
    def __init__(self, q=1, V=1):
        super().__init__()
        self.source = None
        self.destination = None
        self.q = q #volumetric flowrate, [l/s]
        self.V = V #total volume, [l]
        self.delay = V/q #residence time, not needed, [s]
        self._length = 1 #reactor length, [m]
        self._N = 100 #number of mesh points
        self._x_grid = np.linspace(0,self._length,self._N) # grid points
        self._dx = self.V/100 #subvolume element volume
        self._store = [] #what is this for?
        self._chemistry = None # chemistry domain (if added)
        self._chemistry_ind = False # chemistry components indexes into self.variables
        self._jac = None #Jacobian(not used so far)

    def _get_c_in(self, t):
        #used to get inlet variables values
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
        dx = self._dx
        if self._chemistry:
            self._chemistry.initial_concentrations(**{k:c_in[self._chemistry_ind][i] for i, k in enumerate(self._chemistry.variables.keys())})
            self._chemistry.time_stop = self.delay
            sol = self._chemistry.run(plot = False, output = True)
            # dmdt[self._chemistry_ind] += self.q/0.0001*([v[-1] for v in sol.y]-c_in[self._chemistry_ind])
        
        dmdt = np.zeros_like(c)
        dmdt = -self.q*(-c_in+c[:,:,0])/self._dx
        
        return dmdt
    
    def inlet(self, **kwargs):
        for k, v in kwargs.items():
            if type(v)==tuple:
                self.initial_values.update({k:v[0]})
                self.events = v[1]
            elif callable(v):
                self.initial_values.update({k:v})
            else:
                self.initial_values.update({k:v})
            self.variables.update({k:0})
        self._iv = [v for v in self.initial_values.values()]
        self._update_vars()
        # self._sort_vars()
    
    @property
    def chemistry(self):
        return self._chemistry

    @chemistry.setter
    def chemistry(self, domain):
        '''
        Adds chemistry domain into the reactor
        '''
        if type(domain) != Chemistry:
            return
        
        self._chemistry = domain
        self._update_vars()
        # self._sort_vars()
        return
    
    def _update_vars(self):
        if self._chemistry:
            ch_vars = set(self._chemistry.variables)
            self_vars = set(self.variables)
            diff = ch_vars-self_vars
            new_vars = dict().fromkeys(diff, 0)
            self.variables.update(new_vars)
            self.initial_values.update(new_vars)
        
            #Indexes of vars present in chemistry domain:
            self._chemistry_ind = [(v in self._chemistry.variables) for v in self.variables]
        return
    