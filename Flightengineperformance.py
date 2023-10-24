import numpy as np
from Testbedengineperformance import A_inlet,eta_comp,eta_turb

print("Conditions: Flight performance")

#Datasheet
total_comp = 5.5 #-
R = 287 #J/K*mol
k_a = 1.4 #-
k_g = 1.33 #-
cp_a = 1000 #J/kg K
cp_g = 1150 #J/kg K
eta_noz = 1
LHV = 43*10**6 #Lower heating value
M = 0.78 #-
h = 35000 #ft
TIT = 1150 #K

#Initial calculations
T_a = 288.15-0.0065*35000*0.3048   #Ambient temperature at flight altitude
P_a = 101325*(T_a/288.15)**(9.80665/(0.0065*R))
rho = P_a/(T_a*R) #kg/m^3
a = np.sqrt(k_a*R*T_a)
v_flight = M*a
m_air = A_inlet*rho*v_flight

def total_conditions(P_a,T_a,k,M):
    P_t = P_a*(1+(k-1)/2*M**2)**(k/(k-1))
    T_t = T_a * (1 + (k - 1) / 2 * M ** 2)
    return P_t, T_t

P_t, T_t = total_conditions(P_a,T_a,k_a,M)

def station_3(P_t,T_t,k,cp_a,m_air,eta_comp):
    P_t2 = P_t
    T_t2 = T_t
    P_t3 = P_t2*total_comp
    T_t3 = T_t2*(1+1/eta_comp*((P_t3/P_t2)**((k-1)/k)-1))
    W_comp = cp_a*m_air*(T_t3-T_t2)
    return P_t3,T_t3,W_comp

P_t3, T_t3, W_comp =station_3(P_t,T_t,k_a,cp_a,m_air,eta_comp)

def station_4(m_air,LHV,cp_g,T_t3,P_t3,TIT):
    P_t4 = P_t3
    m_fuel = m_air*cp_g*(TIT-T_t3)/(LHV)
    return P_t4,m_fuel

P_t4,m_fuel = station_4(m_air,LHV,cp_g,T_t3,P_t3,TIT)

def station_5(m_air,m_fuel,cp_g,TIT,W_comp,eta_turb,P_t4):
    T_t5 = TIT - W_comp/((m_air+m_fuel)*cp_g)
    P_t5 = ((((T_t5 / TIT - 1) / -eta_turb) - 1) * -1) ** (k_g / (k_g - 1)) * P_t4
    return P_t5, T_t5

P_t5,T_t5 = station_5(m_air,m_fuel,cp_g,TIT,W_comp,eta_turb,P_t4)
P_t7,T_t7 = P_t5,T_t5

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

print("Gross thrust at altitude:",F_G_calc)

def area_ratio(P_t7,P_t4,k_g):
    A_t_to_A_n = (P_t7/P_t4)**((k_g+1)/(2*k_g))
    return A_t_to_A_n


A_t_to_A_n = area_ratio(P_t7,P_t4,k_g)

def pressure_ratio(area_ratio,k):
    pressureratio = area_ratio**((2*k)/(k+1))
    return pressureratio

Pt_07_to_Pt_04 = pressure_ratio(A_t_to_A_n,k_g)

print("The turbine to propulsive nozzle area ratio is:", A_t_to_A_n)
print("The turbine to propulsive nozzle pressure ratio is:", Pt_07_to_Pt_04)
print("The nozzle is ",nozzle_status)
