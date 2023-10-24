import numpy as np
from Testbedengineperformance import P_a,T_a,T_t2,T_t4,R,k_a,cp_a,k_g,cp_g,LHV,total_comp,eta_comp,eta_turb,eta_noz,m_air #Call values from testbed code
print("Conditions: Flight performance")
#Data
M = 0.78
TIT_cruise = 1150
T_a_cruise = 288.15-0.0065*35000*0.3048   #Ambient temperature at flight altitude
P_a_cruise = 101325*(T_a_cruise/288.15)**(9.80665/(0.0065*R)) #Ambient pressure at flight condition

#Calculate pressure and temperature at inlet of compressor
def total_conditions(P_a,T_a,k,M):
    P_t = P_a*(1+(k-1)/2*M**2)**(k/(k-1))
    T_t = T_a * (1 + (k - 1) / 2 * M ** 2)
    #Inlet efficiency is assumed to be 1 such that:
    P_t2 = P_t
    T_t2 = T_t
    return P_t2, T_t2

P_t2_cruise, T_t2_cruise = total_conditions(P_a_cruise,T_a_cruise,k_a,M)

#Function to calculate compressor ratio during cruise
def compressor_ratio(T_t2,T_t4,total_comp,eta_comp,k,T_t2_cruise,TIT_cruise):
    tau_ratio_testbed = T_t2/T_t4
    tau_comp_testbed = total_comp**((k-1)/(eta_comp*k))
    tau_ratio_cruise = T_t2_cruise/TIT_cruise
    tau_cruise = 1 + tau_ratio_testbed/tau_ratio_cruise*(tau_comp_testbed-1)
    total_comp_cruise = tau_cruise**(eta_comp*k/(k-1))
    return total_comp_cruise

total_comp_cruise = compressor_ratio(T_t2,T_t4,total_comp,eta_comp,k_a,T_t2_cruise,TIT_cruise)

#Function to calculate corrected massflow and the actual massflow
def massflows(P_a,T_a,m_air,total_comp_cruise,total_comp,T_t4,T_t2,TIT_cruise,T_t2_cruise, P_t2_cruise):
    m_corrected = m_air*((total_comp_cruise/total_comp)*(np.sqrt(T_t4/T_t2)/np.sqrt(TIT_cruise/T_t2_cruise)))
    delta = P_t2_cruise / P_a
    theta = T_t2_cruise / T_a
    mflow = m_corrected * delta / np.sqrt(theta)
    return m_corrected, mflow

m_corrected, m_air_cruise = massflows(P_a,T_a,m_air,total_comp_cruise,total_comp,T_t4,T_t2,TIT_cruise,T_t2_cruise, P_t2_cruise)

#Calculate pressure, temperature and work done at outlet of compressor
def station_3(P_t2,T_t2,k,cp_a,m_air,total_comp,eta_comp):
    P_t3 = P_t2*total_comp
    T_t3 = T_t2*(1+1/eta_comp*((P_t3/P_t2)**((k-1)/k)-1))
    W_comp = cp_a*m_air*(T_t3-T_t2)
    return P_t3,T_t3,W_comp

P_t3_cruise, T_t3_cruise, W_comp_cruise = station_3(P_t2_cruise,T_t2_cruise,k_a,cp_a,m_air_cruise,total_comp_cruise,eta_comp)

#Calculate pressure and fuel flow at inlet of turbine
def station_4(m_air,LHV,cp_g,T_t3,P_t3,TIT):
    P_t4 = P_t3
    m_fuel = m_air*cp_g*(TIT-T_t3)/(LHV)
    return P_t4,m_fuel

P_t4_cruise,m_fuel_cruise = station_4(m_air_cruise,LHV,cp_g,T_t3_cruise,P_t3_cruise,TIT_cruise)

#Calculate pressure and temperature at outlet of turbine
def station_5(m_air,m_fuel,cp_g,T_t4,W_comp,eta_turb,P_t4):
    T_t5 = T_t4 - W_comp/((m_air+m_fuel)*cp_g)
    P_t5 = ((((T_t5 / T_t4 - 1) / -eta_turb) - 1) * -1) ** (k_g / (k_g - 1)) * P_t4
    #Duct efficiency is 1
    T_t7 = T_t5
    P_t7 = P_t5
    return P_t7, T_t7

P_t7_cruise,T_t7_cruise = station_5(m_air_cruise,m_fuel_cruise,cp_g,TIT_cruise,W_comp_cruise,eta_turb,P_t4_cruise)

#Calculate gross thrust and determine if nozzle is choked or unchoked
def station_8(eta_noz,k_g,P_t7,T_t7,R,m_air,m_fuel,cp_g,P_a):
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

F_G_cruise, nozzle_status = station_8(eta_noz,k_g,P_t7_cruise,T_t7_cruise,R,m_air_cruise,m_fuel_cruise,cp_g,P_a_cruise)

#Calculate the area ratio
def area_ratio(P_t7,P_t4,k_g):
    A_t_to_A_n = (P_t7/P_t4)**((k_g+1)/(2*k_g))
    return A_t_to_A_n

A_t_to_A_n = area_ratio(P_t7_cruise,P_t4_cruise,k_g)

print("Compressor ratio during cruise:", total_comp_cruise)
print("Gross thrust at flight condition is:",F_G_cruise,'N')
print("The nozzle is",nozzle_status)
print("The turbine to propulsive nozzle area ratio is:", A_t_to_A_n)
print("-----------------------------------------------------------")