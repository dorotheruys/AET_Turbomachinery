import numpy as np
import matplotlib.pyplot as plt
from Flightengineperformance import T_a, P_a, M

R_constant = 287        #J/K*mol
k_a = 1.4               #-
k_g = 1.33              #-
cp_a = 1000             #J/kg K
cp_g = 1150             #J/kg K
m_air = 23.81           #kg/s
total_PR = 5.5          #[-]

def total_conditions(P_a, T_a, k, M):
    P_t = P_a*(1+(k-1)/2*M**2)**(k/(k-1))
    T_t = T_a * (1 + (k - 1) / 2 * M ** 2)
    return P_t, T_t
def station_3(P_t,T_t,k,cp_a,m_air,efficiency_comp, total_PR):
    P_t2 = P_t
    T_t2 = T_t
    P_t3 = P_t2*total_PR
    T_t3 = T_t2*(1+1/efficiency_comp*((P_t3/P_t2)**((k-1)/k)-1))
    W_comp = cp_a*m_air*(T_t3-T_t2)
    return P_t3,T_t3,W_comp
def flow_angle_equations(vars, work_coef, flow_coef, R):
    alpha1, beta2 = vars
    eq1 = 1 - flow_coef * (np.tan(beta2) - np.tan(alpha1)) - work_coef
    eq2 = (-work_coef / 2) - flow_coef * np.tan(alpha1) + 1 - R
    return [eq1, eq2]

def find_flow_angles_alpha1_beta2(work_coef, flow_coef, R):
    # initial_guess = np.array([1., -1.])  # Initial guess for alpha1 and beta2
    # solution = fsolve(flow_angle_equations, initial_guess, args=(work_coef, flow_coef, R))
    #
    # alpha1 = solution[0]
    # beta2 = solution[1]

    tanalpha1 = -1/flow_coef * (R + work_coef/2 - 1)
    alpha1 = np.arctan(tanalpha1)
    tanbeta2 = 1/flow_coef * (work_coef - 1 + flow_coef * np.tan(alpha1))
    beta2 = np.arctan(tanbeta2)

    return alpha1, beta2

def determine_velocity_triangles(flow_coef, work_coef, U, alpha1, beta2, specific_work):

    # Stays constant for all stages: Vm, U, beta2, alpha1, Vt
    V_m = U * flow_coef                 # Definition flow coefficient, stays constant over rotor blade

    # W_2test = V_m / np.cos(beta2)           # cos(beta2) = V_m/W_2
    # V_t2test = U + (V_m * np.tan(beta2))
    # alpha2_geo = np.arctan(V_t2test / V_m)      # tan(alpha2) = V_t2/V_m
    # V_2test = V_m / np.cos(alpha2_geo)          # cos(alpha2) = V_m/V_2
    #
    # V_1test = V_m / np.cos(alpha1)          # cos(alpha1) = V_m/V_1
    # V_t1test = V_m * np.tan(alpha1)         # tan(alpha1) = V_t1/V_m
    # beta1_geo = np.arctan((-U + V_t1test)/V_m)
    # W_1test = V_m / np.cos(beta1_geo)           # cos(beta1) = V_m/W_1

    beta1 = np.arctan(np.tan(alpha1)-1/flow_coef)
    alpha2 = np.arctan(np.tan(beta2)+1/flow_coef)

    print(f"For \u03C8 = {work_coef},\u03C6 = {flow_coef} and R = {R}: ")
    print("Alpha 1 = ", np.rad2deg(alpha1))
    print("Alpha 2 = ", np.rad2deg(alpha2))
    print("Beta 1 = ", np.rad2deg(beta1))
    print("Beta 2 = ", np.rad2deg(beta2))
    print("U = ", U)

    W_1 = V_m / np.cos(beta1)
    V_t1 = V_m * np.tan(alpha1)
    V_1 = V_m / np.cos(alpha1)

    V_2 = V_m / np.cos(alpha2)
    V_t2 = V_m * np.tan(alpha2)
    W_2 = V_m / np.cos(beta2)

    delta_Vt = (V_t2 - V_t1)

    # checks
    # if (work_coef - (1-flow_coef*np.tan(alpha1)+flow_coef*np.tan(beta2))) <= 10e-5:
    #     print('Flow coef true')
    # else:
    #     print("Flow coef false")
    #     print('Difference: ', (work_coef-(1-flow_coef*np.tan(alpha1)+flow_coef*np.tan(beta2))))
    #
    # specific_work_diff = specific_work - (U * delta_Vt)
    # if abs(specific_work_diff) <= 10e-5:
    #     print('Specific work true')
    # else:
    #     print('Specific work false')
    #     print("Difference: ", specific_work_diff)
    return V_1, W_1, beta1, V_2, W_2, alpha2, V_m

def plot_velocity_triangles(V_1, W_1, alpha1, beta1, V_2, W_2, alpha2, beta2, U, pos):
    angs = [alpha1, beta1, alpha2, beta2, np.pi/2, np.pi/2]
    velos = [V_1, W_1, V_2, W_2, U, U]
    labels = ['V1', 'W1', 'V2', 'W2', 'U1', 'U2']
    colours = ['red', 'blue', 'orange', 'darkblue', 'lightgreen', 'green']
    for i in range(len(angs)):
        dx = velos[i] * np.sin(angs[i])
        dy = -1 * velos[i] * np.cos(angs[i])
        labelnow = labels[i]
        colournow = colours[i]

        # Make plotting factor -1 to have the right direction of U
        if alpha1 < beta1:
            plotting_factor = -1
        else:
            plotting_factor = 1

        if i < 4:
            pos.plot([0,dx], [0, dy], label=labelnow, color=colournow)
        elif i == 4:
            x0 = velos[1] * np.sin(angs[1])
            y0 = -1 * velos[1] * np.cos(angs[1])
            pos.plot([x0, x0 + plotting_factor*(dx)], [y0, y0 + dy], label=labelnow, color=colournow)
        elif i == 5:
            x0 = velos[3] * np.sin(angs[3])
            y0 = -1 * velos[3] * np.cos(angs[3])
            pos.plot([x0, x0 + plotting_factor*(dx)], [y0, y0 + dy], label=labelnow, color=colournow)

    pos.legend()
    pos.grid()
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

def calculate_thermodynamic_properties_stage(C_p, eff, k, Tt_1, Pt_1, mass_flow, V_m, V_2, W_rotor1, W_rotor2, specific_work_stage):
    # Definition of numbering: 1 before rotor, 2 between rotor and stator, 3 after stator
    # Use conservation of energy & relative enthalpy conserved over a rotor

    # after rotor
    #Tt_2 = (W_rotor1**2 - W_rotor2**2)/(2*C_p) + Tt_1
    Tt_2 = specific_work_stage/C_p + Tt_1
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

    # entropy, enthalpy
    delta_S = cp_a * np.log(Tt_3 / Tt_1) - R_constant * np.log(Pt_3 / Pt_1)
    delta_h = cp_a * (Tt_3 - Tt_1)

    return Tt_2, Pt_2, P_2, Tt_3, Pt_3, P_3, area_2, area_3, delta_S, delta_h

def calculate_blade_height(area, r_in):
    blade_height = np.sqrt(area / np.pi + r_in**2) - r_in
    return blade_height

def plot_meridional_gaspath(stage_array, bladeheight_array, inner_radius, blade_chord):
    spacing = 0.005
    last_stage = stage_array[-1]
    for i in np.arange(1, last_stage + 1, 1):
        i = int(i)
        start_x = (i - 1) * blade_chord + (i - 1) * spacing
        start_y = inner_radius + bladeheight_array[2*i-1]
        vertices_x = [start_x, start_x + blade_chord, start_x + blade_chord, start_x, start_x]
        vertices_y = [start_y, start_y,-1 * start_y, -1 * start_y, start_y]

        plt.plot(vertices_x, vertices_y, '--bo')
    plt.show()
    return

def plot_h_s_diagram(entropy, enthalpy, pos):
    pos.plot(entropy, enthalpy)
    pos.xlabel(r"$\Delta$ Entropy [J/kg]")
    pos.ylabel(r"$\Delta$ Enthalpy [J]")
    pos.grid()
    return
def compressor_design(eta_comp, total_PR, mass_flow, R, work_coef, flow_coef, r_in):
    # Determine total outlet conditions T0_out and P0_out
    Pt_in, Tt_in = total_conditions(P_a, T_a, k_a, M)
    Pt_out, Tt_out, W_comp = station_3(Pt_in, Tt_in, k_a, cp_a, m_air, eta_comp, total_PR)

    #specific work
    specific_work = cp_a * (Tt_out - Tt_in)

    # Determine number of stages and peripheral speed U with limit of 600 m/s
    nr_stages, U = calculate_number_of_stages(specific_work, work_coef, 600.)

    # Determine flow angles & velocity triangles for each stage
    alpha1, beta2 = find_flow_angles_alpha1_beta2(work_coef, flow_coef, R)  # output in radians
    V_1, W_1, beta1, V_2, W_2, alpha2, V_m = determine_velocity_triangles(flow_coef, work_coef, U, alpha1, beta2, specific_work/nr_stages)

    # Run through each stage and determine the thermodynamic properties
    Tt_before = Tt_in
    Pt_before = Pt_in
    entropy = 0.
    enthalpy = 0.

    temp_array = np.zeros(2*nr_stages+1)
    total_pressure_array = np.zeros(2*nr_stages+1)
    static_pressure_array = np.zeros(2*nr_stages+1)
    PR_array = np.zeros(2*nr_stages+1)
    stagenr_array = np.zeros(2*nr_stages+1)
    bladeheight_array = np.zeros(2*nr_stages+1)
    entropy_array = np.zeros(nr_stages+1)
    enthalpy_array = np.zeros(nr_stages+1)

    stagenr_array[0] = 0
    temp_array[0] = Tt_in/Tt_in
    total_pressure_array[0] = Pt_in/Pt_in
    bladeheight_array[0] = 0
    enthalpy_array[0] = enthalpy
    entropy_array[0] = entropy

    for stagenr in range(1, nr_stages+1):
        Tt_during, Pt_during, P_during, Tt_after, Pt_after, P_after, area_rotor, area_stator, delta_S, delta_h = calculate_thermodynamic_properties_stage(cp_a, eta_comp, k_a, Tt_before, Pt_before, mass_flow, V_m, V_2, W_1, W_2, specific_work/nr_stages)

        entropy += delta_S
        enthalpy += delta_h

        bladeheight_rotor = calculate_blade_height(area_rotor, r_in)
        bladeheight_stator = calculate_blade_height(area_stator, r_in)

        print(f"Condtions at stage {stagenr}: ")
        print(f"Total temperature after rotor [K] = {Tt_during}")
        print(f"Total pressure after rotor [Pa] = {Pt_during}")


        stagenr_array[2*(stagenr-1)+1] = stagenr - 0.5
        stagenr_array[2*(stagenr-1)+1+1] = stagenr
        temp_array[2*(stagenr-1)+1] = Tt_during/Tt_in
        temp_array[2*(stagenr-1)+1+1] = Tt_after/Tt_in
        total_pressure_array[2*(stagenr-1)+1] = Pt_during/Pt_in
        total_pressure_array[2*(stagenr-1)+1+1] = Pt_after/Pt_in
        static_pressure_array[2*(stagenr-1)+1] = P_during/Pt_in
        static_pressure_array[2*(stagenr-1)+1+1] = P_after/Pt_in
        PR_array[2*(stagenr-1)+1] = Pt_after/Pt_before
        PR_array[2*(stagenr-1)+1+1] = Pt_after/Pt_before
        bladeheight_array[2*(stagenr-1)+1] = bladeheight_rotor
        bladeheight_array[2*(stagenr-1)+1+1] = bladeheight_stator

        entropy_array[stagenr] = entropy
        enthalpy_array[stagenr] = enthalpy

        Tt_before = Tt_after
        Pt_before = Pt_after

    print('Number of stages: ', nr_stages)
    print('Blade Heights: ', bladeheight_array)
    print('Pt/Pt0: ', total_pressure_array)
    print('Tt/Tt0: ', temp_array)

    fig,axs = plt.subplots(2, 2)
    axs[0, 0].plot(stagenr_array, total_pressure_array)
    axs[0, 0].set_xticks([i for i in np.arange(0,int(len(stagenr_array)/2),0.5)])
    axs[0, 0].set_yticks([i for i in np.arange(min(total_pressure_array), max(total_pressure_array)+0.5, 0.5)])
    axs[0, 0].grid()
    axs[0, 0].set_xlabel('Stage Number')
    axs[0, 0].set_ylabel(r'$P_{t}/P_{t0}$ [-]')

    axs[0, 1].plot(stagenr_array, static_pressure_array)
    axs[0, 1].set_xticks([i for i in np.arange(0,int(len(stagenr_array)/2 ),0.5)])
    axs[0, 1].set_yticks([i for i in np.arange(min(static_pressure_array), max(static_pressure_array), 0.25)])
    axs[0, 1].grid()
    axs[0, 1].set_xlabel('Stage Number')
    axs[0, 1].set_ylabel(r'$P/P_{t0}$ [-]')

    axs[1, 0].plot(stagenr_array, temp_array)
    axs[1, 0].set_xticks([i for i in np.arange(0,int(len(stagenr_array)/2),0.5)])
    axs[1, 0].set_yticks([i for i in np.arange(min(temp_array), max(temp_array), 0.1)])
    axs[1, 0].grid()
    axs[1, 0].set_xlabel('Stage Number')
    axs[1, 0].set_ylabel(r'$T_{t}/T_{t0}$ [-]')

    axs[1, 1].plot([stagenr_array[2], stagenr_array[4]], [PR_array[2], PR_array[4]])
    axs[1, 1].set_xticks([i for i in np.arange(1,2.5, 0.5)])
    axs[1, 1].set_yticks([i for i in np.arange(2, 2.6, 0.1)])
    axs[1, 1].grid()
    axs[1, 1].set_xlabel('Stage Number')
    axs[1, 1].set_ylabel(r'$\beta_{stage}$ [-]')

    fig.set_size_inches(9.5, 6.5)
    fig.savefig("thermodynamic-properties-stage-flight.svg", format="svg", bbox_inches='tight')
    #plt.show()

    fig2 = plt.figure(2)
    plot_velocity_triangles(V_1, W_1, alpha1, beta1, V_2, W_2, alpha2, beta2, U, plt)
    fig2.savefig("velocity-triangles-flight.svg", format="svg", bbox_inches='tight')
    #plt.show()

    fig3 = plt.figure(3)
    plot_h_s_diagram(entropy_array, enthalpy_array, plt)
    fig3.savefig("h-s-diagram-flight.svg", format="svg", bbox_inches='tight')
    plt.show()


    # axs[1, 2].plot(stagenr_array, bladeheight_array, label='Blade Height')
    # axs[1, 2].set_xticks([i for i in range(1,int(len(stagenr_array)/2 + 1),1)])
    # axs[1, 2].set_yticks([i for i in np.arange(min(bladeheight_array), max(bladeheight_array), 0.01)])
    # axs[1, 2].grid()
    # axs[1, 2].set_xlabel('Stage Number')
    # axs[1, 2].set_ylabel('Blade Height [m]')

    return stagenr_array, bladeheight_array


# # Example usage:
work_coef = 0.38
flow_coef = 0.77
eta_comp = 0.91

# work_coef = 0.1
# flow_coef = 0.2
# eta_comp = 1.

R = 0.5
r_in = 0.05
blade_chord = 0.04      #[m]

stagenr_array, bladeheight_array = compressor_design(eta_comp, total_PR, m_air, R, work_coef, flow_coef, r_in)
plot_meridional_gaspath(stagenr_array, bladeheight_array, r_in, blade_chord)
plt.show()