import math

class Constants:
    k = 1.381 * math.pow(10, -23) #Boltzmann's constant [J/K] -- can be changed to [eV/K]
    R = 8.314 #Ideal gas constant [J/(mol*K)] -- can be changed to [(l*atm)/(mol*K)]
    h = 6.626 * math.pow(10, -34) #Planck's constant [J*s] -- can be changed to [eV*s]
    h_bar = 1.055 * math.pow(10, -34) #Planck's constant over two pi [J*s] -- can be changed to [eV*s]
    c = 3 * math.pow(10, 8) #speed of light in a vacuum [m/s]
    mu_e = 9.2848 * math.pow(10, -24) #magnetic moment of an electron [J/T]
    mu_p = 1.4106 * math.pow(10, -26) #magnetic moment of a proton [J/T]
    m_e = 9.109 * math.pow(10, -31) #mass of an electron [kg]
    m_p = 1.673 * math.pow(10, -27) #mass of a proton [kg]
    g = 9.81; #gravitational acceleration near the surface of the Earth [m/s^2]
    Na = 6.022 * math.pow(10, 23) #Avagadro's number [/mol]
    sigma_S_B = 5.67 * math.pow(10, -8) #Stephan-Boltzmann constant [W/(m^2*K^4)]

class Equations:
    '''
    Nothing here
    '''
