import numpy as np

def plot_velocity_triangles():

    return
def compressor_design(eta_comp, P0_in, T0_in, rot_speed, pressure_ratio, mass_flow, R, work_coef, flow_coef):
    # Determine total outlet conditions T0_out and P0_out

    #specific work
    specific_work = cp_a * (T0_out - T0_in)

    # Determine peripheral speed U
    U = 0.  #[m/s]
    nr_stages = 1
    while U <= 600.:
        specific_work_stage = specific_work/nr_stages
        U = np.sqrt(specific_work_stage/work_coef)
        nr_stages += 1

    # Determine axial velocity component V_m
    V_m = U * flow_coef

    return
