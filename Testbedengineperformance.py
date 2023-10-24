import numpy as np
print("Conditions: Testbed performance")
#Data
m_air = 23.81 #kg/s
m_fuel = 0.4267 #kg/s
total_comp = 5.5 #-
T_a = 288 #K
P_a = 100000 #Pa
R = 287 #J/K*mol
k_a = 1.4 #-
k_g = 1.33 #-
cp_a = 1000 #J/kg K
cp_g = 1150 #J/kg K
eta_noz = 1
LHV = 43*10**6 #Lower heating value

#Estimated values (These are altered such that the gross thrust is 15167 N at sea level)
eta_comp = 0.8065         #Efficiency of the compressor is changed
eta_turb = eta_comp*0.957 #Efficiency of the turbine is assumed to have a linear relationship with that of the compressor based on the first assignment

#Total conditions are assumed to be similar
P_t2 = P_a
T_t2 = T_a

#Calculate pressure, temperature and work done at outlet of compressor
def station_3(P_t2,T_t2,k,cp_a,m_air,total_comp,eta_comp):
    P_t3 = P_t2*total_comp
    T_t3 = T_t2*(1+1/eta_comp*((P_t3/P_t2)**((k-1)/k)-1))
    W_comp = cp_a*m_air*(T_t3-T_t2)
    return P_t3,T_t3,W_comp

P_t3, T_t3, W_comp = station_3(P_t2,T_t2,k_a,cp_a,m_air,total_comp,eta_comp)

#Calculate pressure and temperature at inlet of turbine
def station_4(m_air,LHV,m_fuel,cp_g,T_t3,P_t3):
    P_t4 = P_t3
    T_t4 = (LHV*m_fuel)/(m_air*cp_g)+T_t3
    return P_t4,T_t4

P_t4,T_t4 = station_4(m_air,LHV,m_fuel,cp_g,T_t3,P_t3)

#Calculate pressure and temperature at outlet of turbine
def station_5(m_air,m_fuel,cp_g,T_t4,W_comp,eta_turb,P_t4):
    T_t5 = T_t4 - W_comp/((m_air+m_fuel)*cp_g)
    P_t5 = ((((T_t5 / T_t4 - 1) / -eta_turb) - 1) * -1) ** (k_g / (k_g - 1)) * P_t4
    return P_t5, T_t5

P_t5,T_t5 = station_5(m_air,m_fuel,cp_g,T_t4,W_comp,eta_turb,P_t4)
#Assume 100% duct efficiency
P_t7,T_t7 = P_t5,T_t5

#Calculate gross thrust and determine if nozzle is choked or unchoked
def station_8(eta_noz,k_g,P_t7,T_t7,R,m_air,m_fuel,cp_g):
    crit_ratio = (1-1/eta_noz*((k_g-1)/(k_g+1)))**(-k_g/(k_g-1))
    if P_t7/P_a>crit_ratio:
        P_8 = P_t7/crit_ratio
        T_8 = T_t7*2/(k_g+1)
        v_8 = np.sqrt(T_8*R*k_g)
        A_8 = (m_air+m_fuel) * R * T_8 / (P_8 * v_8)  # Exhaust area
        F_G_calc = (m_air+m_fuel)*v_8 + A_8 * (P_8 - P_a)
        nozzle_status = 'choked'
    else:
        v_8 = np.sqrt(2 * cp_g * eta_noz * T_t7 * (1 - (P_t7 / P_a) ** (-(k_g - 1) / k_g)))
        F_G_calc = (m_air+m_fuel)*v_8
        nozzle_status = 'not choked'
    return F_G_calc, nozzle_status

F_G_calc, nozzle_status = station_8(eta_noz,k_g,P_t7,T_t7,R,m_air,m_fuel,cp_g)

#Calculate the area ratio
def area_ratio(P_t7,P_t4,k_g):
    A_t_to_A_n = (P_t7/P_t4)**((k_g+1)/(2*k_g))
    return A_t_to_A_n

A_t_to_A_n = area_ratio(P_t7,P_t4,k_g)

print("Gross thrust at sea-level is:",F_G_calc,'N')
print("The nozzle is",nozzle_status)
print("The turbine to propulsive nozzle area ratio is:", A_t_to_A_n)
print("-----------------------------------------------------------")