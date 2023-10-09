import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from Flightengineperformance import total_conditions, station_3, T_a, P_a, M, m_air

R_constant = 287        #J/K*mol
k_a = 1.4               #-
k_g = 1.33              #-
cp_a = 1000             #J/kg K
cp_g = 1150             #J/kg K
def flow_angle_equations(vars, work_coef, flow_coef, R):
    alpha1, beta2 = vars
    eq1 = 1 - flow_coef * (-np.tan(beta2) + np.tan(alpha1)) - work_coef
    eq2 = (-work_coef / 2) - flow_coef * np.tan(alpha1) + 1 - R
    return [eq1, eq2]

def find_flow_angles_alpha1_beta2(work_coef, flow_coef, R):
    initial_guess = np.array([0.0, 0.0])  # Initial guess for alpha1 and beta2
    solution = fsolve(flow_angle_equations, initial_guess, args=(work_coef, flow_coef, R))

    alpha1 = solution[0]
    beta2 = solution[1]

    # checks
    print('Checks')
    if beta2 < 0.:
        print('Beta2 true')
    else:
        print('Beta2 false')

    if alpha1 > 0.:
        print('Alpha1 true')
    else:
        print('Alpha1 false')
    return alpha1, beta2

def determine_velocity_triangles(flow_coef, U, alpha1, beta2, specific_work_stage):
    # Stays constant for all stages: Vm, U, beta2, alpha1, Vt
    V_m = U * flow_coef                 # Definition flow coefficient, stays constant over rotor blade
    W_2 = V_m / np.cos(beta2)           # cos(beta2) = V_m/W_2
    V_t2 = U - abs(V_m * np.tan(beta2))
    alpha2 = np.arctan(V_t2 / V_m)      # tan(alpha2) = V_t2/V_m
    V_2 = V_m / np.cos(alpha2)          # cos(alpha2) = V_m/V_2

    V_1 = V_m / np.cos(alpha1)          # cos(alpha1) = V_m/V_1
    V_t1 = V_m * np.tan(alpha1)         # tan(alpha1) = V_t1/V_m
    beta1 = np.arctan((U - V_t1)/V_m)
    W_1 = V_m / np.cos(beta1)           # cos(beta1) = V_m/W_1

    delta_Vt = (V_t2 - V_t1)

    # checks
    if (work_coef - (1-flow_coef*np.tan(alpha1)+flow_coef*np.tan(beta2))) <= 10e-6:
        print('Flow coef true')
    else:
        print("Flow coef false")
        print('Difference: ', (work_coef-(1-flow_coef*np.tan(alpha1)+flow_coef*np.tan(beta2))))

    if abs(specific_work_stage - (U * delta_Vt)) <= 10e-6:
        print('Specific work true')
    else:
        print('Specific work false')
        print("Difference: ", abs(specific_work_stage - (U * delta_Vt)))

    if alpha2 > 0.:
        print('Alpha2 true')
    else:
        print('Alpha2 false')

    if beta1 > 0.:
        print('Beta1 true')
    else:
        print('Beta1 false')
    return V_1, W_1, beta1, V_2, W_2, alpha2, V_m

def plot_velocity_triangles(V_1, W_1, alpha1, beta1, V_2, W_2, alpha2, beta2, U):
    angs = [alpha1, beta1, alpha2, beta2, np.pi, np.pi]
    velos = [V_1, W_1, V_2, W_2, U, U]
    for i in range(len(angs)):
        dx = velos[i] * np.sin(angs[i])
        dy = velos[i] * np.cos(angs[i])
        if i < 4:
            plt.arrow(0, 0, dx, dy)
        elif i == 4:
            x0 = velos[1] * np.sin(angs[1])
            y0 = velos[1] * np.cos(angs[1])
            plt.arrow(x0, y0, dx, dy)
        elif i == 4:
            x0 = velos[3] * np.sin(angs[3])
            y0 = velos[3] * np.cos(angs[3])
            plt.arrow(x0, y0, dx, dy)

    plt.show()
    return

def calculate_number_of_stages(specific_work, work_coef, target_U):
    nr_stages = 1  # Initialize the number of stages to 1
    specific_work_stage_initial = specific_work / nr_stages
    initial_U = np.sqrt(specific_work_stage_initial / work_coef)

    # Check if the initial_U is already below the target_U
    if initial_U < target_U:
        return nr_stages

    while True:
        specific_work_stage = specific_work / nr_stages
        U = np.sqrt(specific_work_stage / work_coef)

        if U <= target_U:
            break

        nr_stages += 1

    return nr_stages, U

def calculate_thermodynamic_properties_stage(C_p, eff, k, Tt_1, Pt_1, mass_flow, V_m, V_2, W_rotor1, W_rotor2):
    # Definition of numbering: 1 before rotor, 2 between rotor and stator, 3 after stator
    # Use conservation of energy & relative enthalpy conserved over a rotor

    # after rotor
    Tt_2 = (W_rotor2**2 + W_rotor1**2)/(2*C_p) + Tt_1
    Pt_2 = Pt_1 * ((Tt_2/Tt_1 - 1) * eff + 1)**(k/(k-1))

    T_2 = Tt_2 - V_2**2/(2*C_p)
    M_2 = V_2/np.sqrt(k * R_constant * T_2)
    P_2 = Pt_2 / (1+(k-1)/2*M_2**2)**(k/(k-1))
    rho_2 = P_2/(R_constant * T_2)
    area_2 = mass_flow / (rho_2 * V_m)

    # after stator
    V_3 = V_2
    Tt_3 = Tt_2     #conservation of enthalpy
    Pt_3 = Pt_2     #isentropic flow

    T_3 = Tt_3 - V_3**2/(2*C_p)
    M_3 = V_3/np.sqrt(k * R_constant * T_3)
    P_3 = Pt_3 / (1+(k-1)/2*M_3**2)**(k/(k-1))
    rho_3 = P_3/(R_constant * T_3)
    area_3 = mass_flow / (rho_3 * V_m)

    return Tt_2, Pt_2, Tt_3, Pt_3, area_2, area_3

def calculate_blade_height(area, r_in):
    blade_height = np.sqrt(area / np.pi + r_in**2) - r_in
    return blade_height

def compressor_design(eta_comp, Pt_in, Tt_in, Tt_out, mass_flow, R, work_coef, flow_coef, r_in):
    # Determine total outlet conditions T0_out and P0_out

    #specific work
    specific_work = cp_a * (Tt_out - Tt_in)

    # Determine number of stages and peripheral speed U with limit of 600 m/s
    nr_stages, U = calculate_number_of_stages(specific_work, work_coef, 600.)

    # Determine flow angles & velocity triangles for each stage
    alpha1, beta2 = find_flow_angles_alpha1_beta2(work_coef, flow_coef, R)  # output in radians
    V_1, W_1, beta1, V_2, W_2, alpha2, V_m = determine_velocity_triangles(flow_coef, U, alpha1, beta2, specific_work/nr_stages)
    plot_velocity_triangles(V_1, W_1, alpha1, beta1, V_2, W_2, alpha2, beta2, U)

    # Run through each stage and determine the thermodynamic properties
    Tt_before = Tt_in
    Pt_before = Pt_in

    temp_array = np.zeros(2*nr_stages)
    pressure_array = np.zeros(2*nr_stages)
    stagenr_array = np.zeros(2*nr_stages)
    bladeheight_array = np.zeros(2*nr_stages)

    for stagenr in range(1, nr_stages+1):
        Tt_during, Pt_during, Tt_after, Pt_after, area_rotor, area_stator = calculate_thermodynamic_properties_stage(cp_a, eta_comp, k_a, Tt_before, Pt_before, mass_flow, V_m, V_2, W_1, W_2)

        bladeheight_rotor = calculate_blade_height(area_rotor, r_in)
        bladeheight_stator = calculate_blade_height(area_stator, r_in)

        stagenr_array[2*(stagenr-1)] = stagenr
        temp_array[2*(stagenr-1)] = Tt_during
        temp_array[2*(stagenr-1)+1] = Tt_after
        pressure_array[2*(stagenr-1)] = Pt_during
        pressure_array[2*(stagenr-1)+1] = Pt_after
        bladeheight_array[2*(stagenr-1)] = bladeheight_rotor
        bladeheight_array[2*(stagenr-1)+1] = bladeheight_stator

        Tt_before = Tt_after
        Pt_before = Tt_after


    return


# # Example usage:
work_coef = 0.38
flow_coef = 0.77
R = 0.5
eta_comp = 1.
r_in = 0.5

Pt_in, Tt_in = total_conditions(P_a, T_a, k_a, M)
Pt_out, Tt_out, W_comp = station_3(Pt_in, Tt_in ,k_a ,cp_a, m_air, eta_comp)

print('------------------------------------------')

compressor_design(eta_comp, Pt_in, Tt_in, Tt_out, m_air, R, work_coef, flow_coef, r_in)