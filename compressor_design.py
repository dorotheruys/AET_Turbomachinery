import numpy as np
from scipy.optimize import fsolve

def flow_angle_equations(vars, work_coef, flow_coef, R):
    alpha1, beta2 = vars
    eq1 = 1 - flow_coef * (-np.tan(beta2) + np.tan(alpha1)) - work_coef
    eq2 = (-work_coef / 2) - flow_coef * np.tan(alpha1) + 1 - R
    return [eq1, eq2]

def find_flow_angles_alpha1_beta2(work_coef, flow_coef, R):
    initial_guess = np.array([0.0, 0.0])  # Initial guess for alpha1 and beta2
    solution = fsolve(flow_angle_equations, initial_guess, args=(work_coef, flow_coef, R))
    return solution[0], solution[1]

def determine_velocity_triangles(flow_coef, U, alpha1, beta2):
    V_m = U * flow_coef                 # Definition flow coefficient, stays constant over rotor blade
    W_2 = V_m / np.cos(beta2)           # cos(beta2) = V_m/W_2
    V_t2 = U - V_m * np.tan(beta2)
    alpha2 = np.arctan(V_t2 / V_m)      # tan(alpha2) = V_t2/V_m
    V_2 = V_m / np.cos(alpha2)          # cos(alpha2) = V_m/V_2

    V_1 = V_m / np.cos(alpha1)          # cos(alpha1) = V_m/V_1
    V_t1 = V_m * np.tan(alpha1)         # tan(alpha1) = V_t1/V_m
    beta1 = np.arctan((U - V_t1)/V_m)
    W_1 = V_m / np.cos(beta1)           # cos(beta1) = V_m/W_1

    return V_1, W_1, beta1, V_2, W_2, alpha2
def plot_velocity_triangles(V_1, W_1, alpha1, beta1, V_2, W_2, alpha2, beta2):
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

def calculate_thermodynamic_properties_stage():
    # Definition of numbering: 0 before rotor, 1 between rotor and stator, 2 after stator

    return

def compressor_design(eta_comp, P0_in, T0_in, rot_speed, pressure_ratio, mass_flow, R, work_coef, flow_coef):
    # Determine total outlet conditions T0_out and P0_out

    #specific work
    specific_work = cp_a * (T0_out - T0_in)

    # Determine number of stages and peripheral speed U with limit of 600 m/s
    nr_stages, U = calculate_number_of_stages(specific_work, work_coef, 600.)

    # Determine flow angles & velocity triangles
    alpha1, beta2 = find_flow_angles_alpha1_beta2(work_coef, flow_coef, R)  # output in radians
    V_1, W_1, beta1, V_2, W_2, alpha2 = determine_velocity_triangles(flow_coef, U, alpha1, beta2)



    return


# # Example usage:
# work_coef = 0.5
# flow_coef = 0.2
# R = 0.1
#
# alpha1, beta2 = find_flow_angles_alpha1_beta2(work_coef, flow_coef, R)
# print(f"Alpha1: {alpha1}, Beta2: {beta2}")

compressor_design()