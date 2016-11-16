'''
JPRelativity
Written by John Carter
Last update 08/29/2016
Can be found on GitHub at: https://github.com/jcarter2010/JCarter2010/tree/master/JPackage as part of the JPackage project
'''

import math
import os

class Constants:
    
    c = 3.0 * math.pow(10, 8) #speed of light in a vacuum [m/s]

    def __init__(self):
        self.docs = self.Load_Docs()

    def Load_Docs(self):
        docs = []
        f_path = os.path.realpath(__file__)
        path = f_path
        if os.name == 'nt':
            path = f_path[:f_path.rfind('\\')]
        else:
            path = f_path[:f_path.rfind('/')]
        with open('{}/docs/Astro_Consts_Docs'.format(path), 'r') as f_in:
            text = f_in.read()
            docs = text.split('\n')
            return(docs)

    def Var_Help(self):
        for doc in self.docs:
            print(doc)

class Equations:

    global eqs
    global docs

    def __init__(self):
        self.eqs = self.Load_Eqs()
        self.docs = self.Load_Docs()
        self.Load_Eqs()

    def Load_Docs(self):
        docs = []
        f_path = os.path.realpath(__file__)
        path = f_path
        if os.name == 'nt':
            path = f_path[:f_path.rfind('\\')]
        else:
            path = f_path[:f_path.rfind('/')]
        with open('{}/docs/Relativity_Eqs_Docs'.format(path), 'r') as f_in:
            text = f_in.read()
            docs = text.split('\n')
            return(docs)

    def Load_Eqs(self):
        eqs = []
        f_path = os.path.realpath(__file__)
        path = f_path
        if os.name == 'nt':
            path = f_path[:f_path.rfind('\\')]
        else:
            path = f_path[:f_path.rfind('/')]
        with open('{}/eqs/Relativity_Eqs'.format(path), 'r') as f_in:
            text = f_in.read()
            eqs = text.split('\n')
        return(eqs)

    def Eq_Help(self):
        for eq in self.eqs:
             print(eq)

    def Var_Help(self):
        for doc in self.docs:
            print(doc)

    def Solve_Equation(self, eq, var, val1=0, val2=0, val3=0, val4=0, val5=0):
        index = (self.eqs).index(eq)
        if index == 0:
            return(self.Gamma_Velocity_Relation(var, val1))
        if index == 1:
            return(self.Beta_Velocity_Relation(var, val1))
        if index == 2:
            return(self.Time_Dilation_Time_Gamma_Relation(var, val1, val2))
        if index == 3:
            return(self.Length_Contraction_Length_Gamma_Relation(var, val1, val2))
        if index == 4:
            return(self.Length_Contraction_Length_Velocity_Time_Relation(var, val1, val2, val3))
        if index == 5:
            return(self.Time_Dilation_Time_Velocity_Length_Relation(var, val1, val2, val3))
        if index == 6:
            return(self.Velocity_Addition_Parallel_Direction_Relation(var, val1, val2))
        if index == 7:
            return(self.Velocity_Addition_Perpendicular_Direction_Relation(var, val1, val2, val3))
        if index == 8:
            return(self.Redshift_Velocity_Relation_Parallel(var, val1))
        if index == 9:
            return(self.Redshift_Velocity_Relation_Perpendicular(var, val1))
        if index ==10:
            return(self.Energy_Mass_Relation(var, val1))

    def Gamma_Velocity_Relation(self, var, val1):
        'gamma = 1 / sqrt(1 - v^2 / c^2)'
        if var == 'gamma': #val1 = v
            return(1 / math.sqrt(1 - val1 ** 2 / Constants.c ** 2))
        if var == 'v': #val1 = gamma
            return(math.sqrt(1 - (1 / val1) ** 2 * Constants.c ** 2))

    def Beta_Velocity_Relation(self, var, val1):
        'beta = v^2 / c^2'
        if var == 'beta': #val1 = v
            return(val1 ** 2 / Constants.c ** 2)
        if var == 'v': #val1 == beta
            return(math.sqrt(val2 * Constants.c ** 2))

    def Time_Dilation_Time_Gamma_Relation(self, var, val1, val2):
        't_prime = gamma t'
        if var == 't_prime': #val1 = t, val2 = gamma
            return(val1 / val2)
        if var == 't': #val1 = t_prime, val2 = gamma
            return(val1 * val2)
        if var == 'gamma': #val1 = t_prime, val2 = t
            return(val2 / val1)

    def Length_Contraction_Length_Gamma_Relation(self, var, val1, val2):
        'l_prime = l / gamma'
        if var == 'l_prime': #val1 = l, val2 = gamma
            return(val1 / val2)
        if var == 'l': #val1 = l_prime, val2 = gamma
            return(val1 * val2)
        if var == 'gamma': #val1 = l_prime, val2 = l
            return(val2 / val1)

    def Length_Contraction_Length_Velocity_Time_Relation(self, var, val1, val2, val3):
        'x_prime = gamma(x - v t)'
        if var == 'x_prime': #val1 = x, val2 = v, val3 = t
            return(self.Gamma_Velocity_Relation('gamma', val2) * (val1 - val2 * val3))
        if var == 'x': #val1 = x_prime, val2 = v, val3 = t
            return((val1 / self.Gamma_Velocity_Relation('gamma', val2)) + val2 * val3)
        if var == 'v': #val1 = x_prime, val2 = x, val3 = t
            return((val2 - (val1 / self.Gamma_Velocity_Relation('gamma', val2))) / val3)
        if var == 't': #val1 = x_prime, val2 = x, val3 = v
            return((val2 - (val1 / self.Gamma_Velocity_Relation('gamma', val2))) / val3)

    def Time_Dilation_Time_Velocity_Length_Relation(self, var, val1, val2, val3):
        't_prime = gamma(t - v x / c^2)'
        if var == 't_prime': #val1 = t, val2 = v, val3 = x
            return(self.Gamma_Velocity_Relation('gamma', val2) * (val1 - (val2 * val3) / Constants.c ** 2))
        if var == 't': #val1 = t_prime, val2 = v, val3 = x
            return((val1 / self.Gamma_Velocity_Relation('gamma', val2)) + (val2 * val3) / Constants.c ** 2)
        if var == 'v': #val1 = t_prime, val2 = t, val3 = x
            return((val2 - (val1 / self.Gamma_Velocity_Relation('gamma', val2))) / val3 * Constants.c ** 2)
        if var == 'x': #val1 = t_prime, val2 = t, val3 = t
            return((val2 - (val1 / self.Gamma_Velocity_Relation('gamma', val2))) / val3 * Constants.c ** 2)

    def Velocity_Addition_Parallel_Direction_Relation(self, var, val1, val2):
        'V_prime_par = (V_par - v) / (1 - (V_par v / c^2))'
        if var == 'V_prime_par':
            return((val1 - val2) / (1 - (val1 * val2 / Constants.c **2)))
        if var == 'V_par': #val1 = V_prime_par, val2 = v
            return((val1 + val2) / ((val1 * val2 / Constants.c ** 2) + 1))
        if var == 'v': #val1 = V_prime_par, val2 = v_par
            return((val1 - val2) / ((val1 * val2 / Constants.c ** 2) - 1))

    def Velocity_Addition_Perpendicular_Direction_Relation(self, var, val1, val2, val3):
        'V_prime_perp = V_perp / (gamma(1 - (V_par v / c^2)))'
        if var == 'V_prime_perp': #val1 = V_perp, val2 = V_par, val3 = v
            return(val1 / (self.Gamma_Velocity_Relation('gamma', val3) * (1 - val2 * val3 / Constants.c ** 2)))
        if var == 'V_perp': #val1 = V_prime_perp, val2 = V_par, val3 = v
            return(val1 * (self.Gamma_Velocity_Relation('gamma', val3) * (1 - val2 * val3 / Constants.c ** 2)))
        if var == 'V_par': #val1 = V_prime_perp, val2 = V_perp, val3 = v
            return(((val2 / val1 / self.Gamma_Velocity_Relation('gamma', val3)) + 1) * Constants.c ** 2 / val3)

    def Redshift_Velocity_Relation_Parallel(self, var, val1):
        'z = v (sqrt(1 - beta))/(sqrt(1 + beta))'
        if var == 'z': #val1 = v
            return(val1 *  (math.sqrt(1 - self.Beta_Velocity_Relation('beta', val1)))/(sqrt(1 + self.Beta_Velocity_Relation('beta', val1))))
        if var == 'v': #val1 = z
            return(self.Beta_Velocity_Relation('v', (1 - val1 ** 2) / (1 + val1 ** 2)))

    def Redshift_Velocity_Relation_Perpendicular(self, var, val1):
        'z = gamma v'
        if var == 'z': #val1 = v
            return(self.Gamma_Velocity_Relation('gamma', val1) * val1)
        if var == 'v': #val1 = z
            return(val1 ** 2 / (1 + (val1 ** 2 / Constants.c ** 2)))

    def Energy_Mass_Relation(self, var, val1):
        'E = m c^2'
        if var == 'E': #val1 = m
            return(val1 * Constants.c ** 2)
        if var == 'm': #val1 = E
             return(val1 / Constants.c ** 2)
