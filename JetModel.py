import math

def jet_model(rpm, egt, eta_c, eta_t, t_0=21, p_0=101325, extra_verbose=False):
    """
    Models performance of the JetJoe 1400 model turbojet. 

    Args:
    rpm: Revolutions per minute of the motor
    egt: Exhaust gas temperature in celsius 
    eta_c: Compressor adiabatic efficiency
    eta_t: Turbine adiabatic efficiency
    t_0 (optional): Ambient temperature in celsius
    p_0 (optional): Ambient pressure in pascals
    extra_verbose (optional): Boolean flag to print extra information
    """
    try:
        # CONVERSION
        IN_TO_M = 1 / 39.37
        C_TO_K = 273.15
        N_TO_LBS = 4.448

        # CONSTANTS
        ALPHA_1 = 62  # deg
        BETA_1 = 27  # deg
        R_1 = 0.671 * IN_TO_M  # m
        R_2 = 1.081 * IN_TO_M  # m
        R_AVG = (R_1 + R_2) / 2  # m
        A_FAN = (math.pi * R_2**2) - (math.pi * R_1**2)  # m^2
        H_FUEL = 43000000  # J/kg
        C_PH = 1244  # J/K-kg
        C_PC = 1004  # J/K-kg
        GAMMA_C = 1.4
        GAMMA_H = 1.3
        R = 287  # J/K-kg

        # INPUTS
        RPM = rpm  # RPM
        EGT = egt + C_TO_K  # K

        # PARAMETERS
        T_0 = t_0 + C_TO_K  # K
        P_0 = p_0  # Pa
        ETA_C = eta_c
        ETA_T = eta_t

        # ALIASES
        T_t5 = EGT  # K
        T_t2 = T_0  # K
        P_t2 = P_0  # Pa

        # Fuel mass fraction
        f = (C_PH * T_t5 - C_PC * T_t2) / (H_FUEL - C_PH * T_t5)

        # Blade rotation rate
        Omega = RPM * (2 * math.pi / 60)

        # Blade speed
        u = Omega * R_AVG

        # Work coefficient
        tan_ratio = math.tan(math.radians(ALPHA_1)) / math.tan(math.radians(BETA_1))
        c_theta1 = (tan_ratio * u) / (1 + tan_ratio)
        lam = c_theta1 / u

        # Specific work of turbine
        w_T = lam * u**2

        # Stagnation temp at turbine entrance
        T_t4 = w_T / C_PH + T_t5

        # Compressor exit temp
        T_t3 = (f + 1) * (C_PH / C_PC) * T_t4 - (f + 1) * (C_PH / C_PC) * T_t5 + T_t2

        # Compressor pressure ratio
        PR_c = math.pow(ETA_C * (T_t3 / T_t2 - 1) + 1, GAMMA_C / (GAMMA_C - 1))

        # Turbine pressure ratio
        PR_t = math.pow(1 - ((1 - T_t5 / T_t4) / ETA_T), GAMMA_H / (GAMMA_H - 1))

        # Compressor exit pressure
        P_t3 = PR_c * P_t2

        # Turbine exit pressure
        P_t5 = PR_t * P_t3

        # Choked NGV mass flow
        m = (
            P_t3
            * math.sqrt(GAMMA_H / (R * T_t4))
            * math.pow(1 + (GAMMA_H - 1) / 2, -(GAMMA_H + 1) / (2 * (GAMMA_H - 1)))
            * A_FAN
        )

        # Exit mach number given isentropic expansion
        M_6 = math.sqrt(
            (2 * math.pow(P_t5 / P_0, (GAMMA_H - 1) / GAMMA_H) - 2) / (GAMMA_H - 1)
        )

        # Static exit temperature
        T_6 = T_t5 / (1 + (GAMMA_H - 1) / 2 * M_6**2)

        # Speed of sound at exit
        a_6 = math.sqrt(GAMMA_H * R * T_6)

        # Exit jet speed
        c_6 = M_6 * a_6

        # Thrust
        F = m * c_6

        # Specific thrust
        ST = (f + 1) * c_6 / math.sqrt(GAMMA_C * R * T_0)

        # Specific fuel consumption
        SFC = f*m / F

        divider = "------------------------------------------------------------------------------------------------"

        if extra_verbose:
            print(divider)
            print("MODEL PARAMETERS")
            print(divider)
            print("{0:<40} {1:<40}".format(f"Revolutions: {RPM} RPM", f"Exhaust Gas Temp: {EGT} K"))
            print("{0:<40} {1:<40}".format(f"Ambient Temp: {T_0} K", f"Ambient Pressure: {P_0} Pa"))
            print("{0:<40} {1:<40}".format(f"Combustor Efficiency: {ETA_C}", f"Turbine Efficiency: {ETA_T}"))
            print(divider)
            print("CYCLE ANALYSIS")
            print(divider)
            print("{0:<25} {1:<25} {2:<25} {3:<25}".format(f"T_t2: {T_t2:.2f} K", f"T_t3: {T_t3:.2f} K", f"T_t4: {T_t4:.2f} K", f"T_t5: {T_t5:.2f} K"))
            print("{0:<25} {1:<25} {2:<25} {3:<25}".format(f"P_t2: {P_t2:.2f} Pa", f"P_t3: {P_t3:.2f} Pa", f"P_t4: {P_t3:.2f} Pa", f"P_t5: {P_t5:.2f} Pa"))
        
        print(divider)
        print("PERFORMANCE VALUES")
        print(divider)
        print("{0:<40} {1:<40}".format(f"PR_c: {PR_c:.2f}", f"PR_t: {PR_t:.2f}"))
        print("{0:<40} {1:<40}".format(f"Mass Flow: {m:.2f} kg/s", f"Turbine Inlet Temp: {T_t4:.2f} K"))
        print("{0:<40} {1:<40}".format(f"Specific Thrust: {ST:.2f}", f"Specific Fuel Consumption {SFC:.2f}"))
        print("{0:<40} {1:<40}".format(f"Thrust: {F:.2f} N ({F/N_TO_LBS:.2f} lbs)", f"Exit Mach {M_6:.2f}"))
        print(divider)

    except:
        print("Model broke :(")


if __name__ == "__main__": 

    # Nominal operating conditions
    jet_model(rpm=160000, egt=680, eta_c=0.6, eta_t=0.65, t_0=21, p_0=101325, extra_verbose=True)
