import math

class Constants:
    k = 1.381 * math.pow(10, -23) #Boltzmann's constant [J/K] -- can be changed to [eV/K]
	R =     8.314 #Ideal gas constant [J/(mol*K)] -- can be changed to [(l*atm)/(mol*K)]
	h =     6.626 * math.pow(10, -34) #Planck's constant [J*s] -- can be changed to [eV*s]
	h_bar = 1.055 * math.pow(10, -34) #Planck's constant over two pi [J*s] -- can be changed to [eV*s]
	c =     3 * math.pow(10, 8) #speed of light in a vacuum [m/s]
	mu_e =  9.2848 * math.pow(10, -24) #magnetic moment of an electron [J/T]
	mu_p =  1.4106 * math.pow(10, -26) #magnetic moment of a proton [J/T]
	m_e =   9.109 * math.pow(10, -31) #mass of an electron [kg]
	m_p =   1.673 * math.pow(10, -27) #mass of a proton [kg]
	g =     9.81; #gravitational acceleration near the surface of the Earth [m/s^2]
	Na =    6.022 * math.pow(10, 23) #Avagadro's number [/mol]
	sigma_S_B =  5.67 * math.pow(10, -8) #Stephan-Boltzmann constant [W/(m^2*K^4)]

    #set k to units of [J/K]
    def set_k_J_per_K():
		k = 1.381 * Math.pow(10, -23)

	#set k to units of [eV/K]
	def set_k_eV_per_K():
		k = 8.617 * Math.pow(10, -5)

	#set R to units [J/(mol*K)]
	def set_R_J_per_mol_K():
		R = 8.314

	#set R to units [(l*atm)/(mol*K)]
	def set_R_l_atm_per_mol_K():
		R = 8.206 * Math.pow(10, -2)

	#set g to units [m/s]
	def set_g_m_per_s():
		g = 9.81

	#set g to units [ft/s]
	def set_g_ft_per_s():
		g = 32.2

	#set h to units [J*s]
	def set_h_J_s():
		h = 6.626 * Math.pow(10, -34)

	#set h to units [eV*s]
	def set_h_eV_s():
		h = 4.136 * Math.pow(10, -15)

	#set h_bar to units [J*s]
	def set_h_bar_J_s():
		h_bar = 1.055 * Math.pow(10, -34)

	#set h_bar to units [eV*s]
	def set_h_bar_eV_s():
		h_bar = 6.58 * Math.pow(10, -16)

class Equations:
    
