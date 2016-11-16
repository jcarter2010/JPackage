import math
import os

class Constants:

    c = 3.0 * math.pow(10, 8) #speed of light in a vacuum [m/s]
    h = 6.63 * math.pow(10, -34) #Planck's constant [J*s]
    k = 1.38 * math.pow(10, -23) #Boltzmann's constant [J/K]
    G = 6.67 * math.pow(10, -11) #gravitational constant [m^3/(s^2kg)]
    g = 9.81 #gravitational acceleration near the surface of the Earth [m/s^2]
    M_Sun = 2 * math.pow(10, 30) #mass of the Sun [kg]
    L_Sun = 3.85 * math.pow(10, 26) #luminosity of the Sun [W]
    D_Sun = 1.496 * math.pow(10, 11) #Distance from the Earth to the Sun [m]
    R_Sun = 6.96 * math.pow(10, 8) #radius of the Sun [m];
    T_Sun = 5.8 * math.pow(10, 3) #temperature of the Sun [K]
    m_hydrogen = 1.67 * math.pow(10, -27) #mass of a hydrogen atom [kg]
    M_Earth = 5.97 * math.pow(10, 24) #mass of the Earth [kg]
    R_Earth = 6.37 * math.pow(10, 6) #radius of the Earth [m]
    sigma_S_B = 5.67 * math.pow(10, -8) #Stephan-Boltzmann constant [W/(m^2*K^4)]
    H_0 = 7.0 * math.pow(10, 1) #Hubble constant [km/s]

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

    def Load_Docs(self):
        docs = []
        f_path = os.path.realpath(__file__)
        path = f_path
        if os.name == 'nt':
            path = f_path[:f_path.rfind('\\')]
        else:
            path = f_path[:f_path.rfind('/')]
        with open('{}/docs/Astro_Eqs_Docs'.format(path), 'r') as f_in:
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
        with open('{}/eqs/Astro_Eqs'.format(path), 'r') as f_in:
            text = f_in.read()
            eqs = text.split('\n')
        return(eqs)

    def Eq_Help(self):
        for eq in self.eqs:
             print(eq)

    def Var_Help(self):
        for doc in self.docs:
            print(doc)

    def Solve_Equation(self, eq, var, val1, val2=0, val3=0, val4=0, val5=0):
        index = (self.eqs).index(eq)
        if index == 0:
            return(self.Force_Mass_Mass_Distance_Relation(var, val1, val2, val3))
        if index == 1:
            return(self.Brightness_Luminosity_Distance_Relation(var, val1, val2))
        if index == 2:
            return(self.Apparent_Magnitude_Absolute_Magnitude_Distance_Relation(var, val1, val2))
        if index == 3:
            return(self.Distance_Brightness_Luminosity_Relation(var, val1, val2))
        if index == 4:
            return(self.Delta_Wavelength_Wavelength_Velocity_Relation(var, val1, val2))
        if index == 5:
            return(self.Energy_Frequency_Relation(var, val1))
        if index == 6:
            return(self.Energy_Wavlength_Relation(var, val1))
        if index == 7:
            return(self.Velocity_Distance_Relation(var, val1))
        if index == 8:
            return(self.Redshift_Delta_Wavelength_Wavelength_Relation(var, val1, val2))
        if index == 9:
            return(self.Redshift_Velocity_Relation(var, val1))
        if index ==10:
            return(self.Kinetic_Energy_Mass_Velocity_Relation(var, val1, val2))
        if index == 11:
            return(self.Velocity_Temperature_Mass_Relation(var, val1, val2))
        if index == 12:
            return(self.Apparent_Magnitude_Apparent_Magnitude_Brightness_Brightness_Relation(var, val1, val2, val3))
        if index == 13:
            return(self.Blue_Visible_Apparent_Magnitude_Blue_Apparent_Magnitude_Visible_Relation(var, val1, val2, val3))
        if index == 14:
            return(self.Blue_Visible_Brightness_Blue_Brightness_Visible_Relation(var, val1, val2, val3))
        if index == 15:
            return(self.Energy_Mass_Relation(var, val1))
        if index == 16:
            return(self.Period_Semimajor_Axis_Relation(var, val1))
        if index == 17:
            return(self.Distance_Parallax_Relation(var, val1))
        if index == 18:
            return(self.Redshift_Velocity_Relation_2(var, val1))
        if index == 19:
            return(self.Time_Dilation_Velcity_Relation(var, val1, val2))
        if index == 20:
            return(self.Length_Contraction_Velocity_Relation(var, val1, val2))
        if index == 21:
            return(self.Schwarzchild_Radius_Mass_Relation(var, val1))
        if index == 22:
            return(self.Luminosity_Radius_Temperature_Relation(var, val1, val2))
        if index == 23:
            return(self.Flux_Temperature_Relation(var, val1))
        if index == 24:
            return(self.Magnification_Focal_Length_Diameter_Relation(var, val1, val2))
        if index == 25:
            return(self.Resolution_Diameter_Relation(var, val1))
        if index == 26:
            return(self.Width_Focal_Length_Apparent_Diameter_Relation(var, val1, val2))
        if index == 27:
            return(self.Limiting_Magnitude_Diameter_Relation(var, va1))
        if index == 28:
            return(self.Tidal_Force_Near_Force_Far_Force_Relation(var, val1, val2))
        if index == 29:
            return(self.Tidal_Force_Mass_Mass_Distance_Exerted_On_Radius_Relation(var, val1, val2, val3, val4))
        if index == 30:
            return(self.Max_Wavelength_Temperature_Relation(car, val1))


    def Force_Mass_Mass_Distance_Relation(self, var, val1, val2, val3):
        if var == 'F': #val1 = M, val2 = m, val3 = R
            return((Constants.G * val1 * val2) / (val3 ** 2))
        if var == 'M': #val1 = F, val2 = m, val3 = R
            return((val1 * val3 **2) / (Constants.G * val2))
        if var == 'm': #val1 = F, val2 = M val3 = R
            return((val1 * val3 **2) / (Constants.G * val2))
        if var == 'R': #val1 = F, val2 = M, val3 = m
            return(math.sqrt((Constants.G * val2 * val3) / (val1)))

    def Brightness_Luminosity_Distance_Relation(self, var, val1, val2):
        if var == 'b': #val1 = L, val2 = d
            return((val1) / (4 * math.pi * (val2 ** 2)))
        if var == 'L': #val1 = b, val2 = d
            return(val1 * 4 * math.pi * (val2 ** 2))
        if var == 'd': #val1 = b, val2 = L
            return(math.sqrt(val2 / (4 * math.pi * val1)))

    def Apparent_Magnitude_Absolute_Magnitude_Distance_Relation(self, var, val1, val2):
        if var == 'm': #val1 = M, val2 = d
            return(val1 + 5 * math.log10(val2 - 5))
        if var == 'M': #val1 = m, val2 = d
            return(val1 - 5 * math.log10(val2 - 5))
        if var == 'd': #val1 = m, val2 = M
            return(math.pow(10, (val1 - val2) / 5) + 5)

    def Distance_Brightness_Luminosity_Relation(self, var, val1, val2):
        if var == 'd': #val1 = L, val2 = b
            return(math.sqrt((val1 / Constants.L_sun) / (val2 / Constants.B_sun)) * Constants.D_sun)
        if var == 'L': #val1 = d, val2 = b
            return(math.pow(val1 / Constants.D_sun, 2) * (val2 / Constants.B_sun) * Constants.L_sun)
        if var == 'b': #val1 = d, val2 = L
            return((val2 / Constants.L_sun) / math.pow(val1 / Constants.D_sun, 2) * Constants.B_sun)

    def Delta_Wavelength_Wavelength_Velocity_Relation(self, var, val1, val2):
        if var == 'delta_lambda': #val1 = lambda_0, val2 = v
            return(val2 * val1 / Constants.c)
        if var == 'lambda_0': #val1 = delta_lambda, val2 = v
            return(val1 * Constants.c / val2)
        if var == 'v': #val1 = delta_lambda, val2 = lambda_0
            return(val1 * Constants.c / val2)

    def Energy_Frequency_Relation(self, var, val1):
        if var == 'E': #val1 = E
            return(Constants.h * val1)
        if var == 'nu': #val1 = E
            return(val1 / Constants.h)

    def Energy_Wavlength_Relation(salf, var, val1):
        if var == 'E': #val1 = lambda
            return(Constants.h * Constants.c / val1)
        if var == 'lambda': #val1 = E
            return(val1 / (Constants.h * Constants.c))

    def Velocity_Distance_Relation(self, var, val1):
        if var == 'v': #val1 = d
            return(Constants.H_0 * val1)
        if var == 'd': #val1 = v
            return(val1 / Constants.H_0)

    def Redshift_Delta_Wavelength_Wavelength_Relation(self, var, val1, val2):
        'z = (lambda - lambda_0) / lambda_0'
        if var == 'z': #val1 = lambda, val2 = lambda_0
            return((val1 - val2) / val2)
        if var == 'lambda': #val1 = z, val2 = lambda_0
            return(val2 * (val1 + 1))
        if var == 'lambda_0': #val1 = z, val2 = lambda
            return(val2 / (val1 + 1))

    def Redshift_Velocity_Relation(self, var, val1):
        'z = v / c'
        if var == 'z': #val1 = v
            return(val1 / Constants.c)
        if var == 'v': #val1 = z
            return(val1 * Constants.x)

    def Kinetic_Energy_Mass_Velocity_Relation(self, var, val1, val2):
        'KE = (1 / 2) m v^2'
        if var == 'KE': #val1 = m, val2 = v
            return((1/2) * val1 * val2 ** 2)
        if var == 'm': #val1 = KE, val2 = v
            return(2 * val1 / (val2 ** 2))
        if var == 'v': #val1 = KE, val2 = m
            return(math.sqrt(2 * val1 / val2))

    def Velocity_Temperature_Mass_Relation(self, var, val1, val2):
        'v = sqrt((3 k T) / m)'
        if var == 'v': #val1 = T, val2 = m
            return sqrt((3 * Constants.k * val1) / val2)
        if var == 'T': #val1 = v, val2 = m
            return(val1 ** 2 * val2 / (3 * Constants.k))
        if var == 'm': #val1 = v, val2 = T
            return(val1 ** 2 / (3 * Constants.k * val2))

    def Apparent_Magnitude_Apparent_Magnitude_Brightness_Brightness_Relation(self, var, val1, val2, val3):
        'm_2 - m_1 = 2.5 log(b_2 / b_1)'
        if var == 'm_1': #val1 = m_2, val2 = b_1, val3 = b_2
            return(-(2.5 * log10(val3 / val2) - val1))
        if var == 'm_2': #val1 = m_1, val2 = b_1, val3 = b_2
            return(2.5 * log10(val3 / val2) + val1)
        if var == 'b_1': #val1 = m_1, val2 = m_2, val3 = b_2
            return(val3 / math.pow(10, ((val2 - val1) / 2.5)))
        if var == 'b_2': #val1 = m_1, val2 = m_2, val3 = b_1
            return(math.pow(10, ((val2 - val1) / 2.5)) * val2)

    def Blue_Visible_Apparent_Magnitude_Blue_Apparent_Magnitude_Visible_Relation(self, var, val1, val2, val3):
        'B - V_color = m_B - m_V'
        if var == 'B': #val1 = V_color, val2 = m_B, val3 = m_V
            return(val2 - val3 + val1)
        if var == 'V_color': #val1 = B, val2 = m_B, val3 = m_V
            return(val3 - val2 + val1)
        if var == 'm_B': #val1 = B, val2 = V_color, val3 = m_V
            return(val1 - val2 + val3)
        if var == 'm_V': #val1 = B, val2 = V_color, val3 = m_B
            return(val2 - val1 + val3)

    def Blue_Visible_Brightness_Blue_Brightness_Visible_Relation(self, var, val1, val2, val3):
        'B - V_color = -2.5 log(b_B / b_V)'
        if var == 'B': #val1 = V_color, val2 = b_B, val3 = b_V
            return(-2.5 * math.log10(val2 / val1) - val1)
        if var == 'V_color': #val = B, val2 = b_B, val3 = b_V
            return(val1 - -2.5 * math.log10(val2 / val1))
        if var == 'b_B': #val1 = B, val2 = V_color, val3 = b_V
            return(math.pow(10, (val1 - val2) / -2.5) * val3)
        if var == 'b_V': #val1 = B, val2 = V_color, val3 = b_B
            return(val3 / math.pow(10, (val1 - val2) / -2.5))

    def Energy_Mass_Relation(self, var, val1):
        'E = m c^2'
        if var == 'E': #val1 = m
            return(val1 * Constants.c ** 2)
        if var == 'm': #val1 = E
            return(val1 / Constants.c ** 2)

    def Period_Semimajor_Axis_Relation(self, var, val1):
        'P^2 = a^3'
        if var == 'P': #val1 = a
            return(sqrt(val1 ** 3))
        if var == 'a': #val1 = P
            return(math.pow(val1 ** 2, 0.33333))

    def Period_Mass_Mass_Semimajor_Axis_Relation(self, var, val1, val2, val3):
        'P^2 = ((4 Pi^2) / (G (m_1 + m_2))) a^3'
        if var == 'P': #val1 = m_1, val2 = m_2, val3 = a
            return(math.sqrt(((4 * math.pi ** 2 * val3 ** 3) / (Constants.G (val1 + val2)))))
        if var == 'm_1': #val1 = P, val2 = m_2, val3 = a
            return((val1 ** 2 * Constants.G) / (4 * math.pi ** 2 * val3 ** 3) - val2)
        if var == 'm_2': #val1 = P, val2 = m_1, val3 = a
            return(va2 - ((val1 ** 2 * Constants.G) / (4 * math.pi ** 2 * val3 ** 3)))
        if var == 'a': #val1 = P, val2 = m_1, val3 = m_2
            return(math.pow((val1 ** 2) / ((4 * math.pi ** 2) / (Constants.G * (val2 + val3))), 0.333333334))

    def Distance_Parallax_Relation(self, var, val1):
        'd = 1 / p'
        if var == 'd': #val1 = p
            return(1 / val1)
        if var == 'p': #val1 = d
            return(1 / val1)

    def Redshift_Velocity_Relation_2(self, var, val1):
        'z = sqrt((c + v) / (c - v)) - 1'
        if var == 'z': #val1 = v
            return(math.sqrt((Constants.c + val1) / (Constants.c - val1)) - 1)
        if var == 'v': #val1 = z
            return((Constants.c * val1 ** 2) / (val1 * 2 - 2))

    def Time_Dilation_Velcity_Relation(self, var, val1, val2):
        't = t_0 / sqrt(1 - (v / c)^2)'
        if var == 'T': #val1 = t_0, val2 = v
            return(val1 / math.sqrt(1 - (val2 / Constants.c) ** 2))
        if var == 't_0': #val1 = t, val2 = v
            return(val1 * math.sqrt(1 - (val2 / Constants.c) ** 2))
        if var == 'v': #val1 = t_0, val2 = t
            return(math.sqrt(-((val1 / val2) ** 2 - 1) * c ** 2))

    def Length_Contraction_Velocity_Relation(self, var, val1, val2):
        'L = L_0 sqrt(1 - (v / c)^2)'
        if var == 'L': #val1 = L_0, val2 = v
            return(val1 * math.sqrt(1 - (val2 / Constants.c) ** 2))
        if var == 'L_0': #val1 = L, val2 = v
            return(math.sqrt(1 - (val2 / Constants.c) ** 2) / val1)
        if var == 'v': #val1 = L, val2 = L_0
            return((1 - (val1 / val2) ** 2) * Constants.c)

    def Schwarzchild_Radius_Mass_Relation(self, var, val1):
        'R_Sch = (2 G m) / c^2'
        if var == 'R_Sch': #val1 = m
            return((2 * Constants.G * val1) / Constants.c ** 2)
        if var == 'm': #val1 = R_Sch
            return((val1 * Constants.c ** 2) / (2 * Constants.G))

    def Luminosity_Radius_Temperature_Relation(self, var, val1, val2):
        'L = 4 Pi R^2 sigma T^4'
        if var == 'L': #val1 = R, val2 = T
            return(4 * math.pi * val1 ** 2 * Constants.sigma_S_B * val2 ** 4)
        if var == 'R': #val1 = L, val2 = T
            return(math.sqrt((val1) / (4 * math.pi * Constants.sigma_S_B * val2 ** 4)))
        if var == 'T': #val1 = L, val2 = R
            return(math.pow((val1) / (4 * math.pi * val2 ** 2 * Constants.sigma_S_B), 0.25))

    def Flux_Temperature_Relation(self, var, val1):
        'F = sigma T^4'
        if var == 'F': #val1 = T
            return(Constants.sigma_S_B * val1 ** 4)
        if var == 'T': #val1 = F
            return(math.pow(val1 / Constants.sigma_S_B, 0.25))

    def Magnification_Focal_Length_Diameter_Relation(self, var, val1, val2):
        'Magnification = f / D'
        if var == 'Magnification': #val1 = f, val2 = D
            return(val1 / val2)
        if var == 'f': #val1 = Magnification, val2 = D
            return(val1 * val2)
        if var == 'D': #val1 = Magnification, val2 = f
            return(val2 / val1)

    def Resolution_Diameter_Relation(self, var, val1):
        'Resolution = 116 / D'
        if var == 'Resolution': #val1 = D
            return(116 / val1)
        if var == 'D': #val1 = Resolution
            return(116 / val1)

    def Width_Focal_Length_Apparent_Diameter_Relation(self, var, val1, val2):
        'w = f theta Pi / 180'
        if var == 'w': #val1 = f, val2 = theta
            return(val1 * val2 * math.pi / 180)
        if var == 'f': #val1 = w, val2 = theta
            return(val1 / (val2 * math.pi / 180))
        if var == 'theta': #val1 = w, val2 = f
            return(val1 / (val2 * math.pi / 180))

    def Limiting_Magnitude_Diameter_Relation(self, var, va1):
        'm = 2.7 + 5 log(D)'
        if var == 'm': #val1 = D
            return(2.7 + 5 * math.log10(val1))
        if var == 'D': #val1 = m
            return(math.pow(10, (val1 + 2.7) / 5))

    def Tidal_Force_Near_Force_Far_Force_Relation(self, var, val1, val2):
        'F_tidal = F_near - F_far'
        if var == 'F_tidal': #val1 = F_near, val2 = F_far
            return(val1 - val2)
        if var == 'F_near': #val1 = F_tidal, val2 = F_far
            return(val1 + val2)
        if var == 'F_far': #val1 = F_tidal, val2 = F_near
            return(val2 - val1)

    def Tidal_Force_Mass_Mass_Distance_Exerted_On_Radius_Relation(self, var, val1, val2, val3, val4):
        'F_tidal = (2 G M m d) / r^3'
        if var == 'F_tidal': #val1 = M, val2 = m, val3 = d, val4 = r
            return((2 * Constants.G * val1 * val2 * val3) / (val4 ** 3))
        if var == 'M': #val1 = F_tidal, val2 = m, val3 = d, val4 = r
            return((val1 * val4 ** 3) / (2 * Constants.G * val2 * val3))
        if var == 'm': #val1 = F_tidal, val2 = M, val3 = d, val4 = r
            return((val1 * val4 ** 3) / (2 * Constants.G * val2 * val3))
        if var == 'd': #val1 = F_tidal, val2 = M, val3 = m, val4 = r
            return((val1 * val4 ** 3) / (2 * Constants.G * val2 * val3))
        if var == 'r': #val1 = F_tidal, val2 = M, val3 = m, val4 = d
            return(math.pow((2 * Constants.G * val2 * val3 * val4) / val1, 0.3333333334))

    def Max_Wavelength_Temperature_Relation(self, car, val1):
        'lambda_max = 0.0029 / T'
        if var == 'lambda_max': #val1 = T
            return(0.0029 / val1)
        if var == 'T': #val1 = lambda_max
            return(0.0029 / val1)
