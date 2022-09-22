from unittest import result
import numpy as np
import scipy as sp
import pandas as pd

"""_summary_
Modelling equation 8 from 'Scaling up silicon photonic based accelerators: challenges and opportunities'

1) For a given bit width no_of_bits and Data rate (Gbps)
2) Find the Photodetector sensitivity (dBm) from Eq 8
3) Fix M and Find the N value such that 10dBm- PLaser(Eq 13) = Min Error. This gives the N value for a given Bitwidth and data rate for AMM
"""


# Finding the Photodetector sensitivity (dBm) from Eq 8


def findPDSensitivity(no_of_bits, DR):
    e = 1.6*10**(-19)  # columbs
    KT = 0.0259*e  # columbs*Volt
    R = 1.2  # A/W
    Id = 35*10**(-9)  # A
    RL = 50  # ohm
    DR = DR*10**(9)  # Bits/s
    RIN = 10**(-140/10)  # power ratio/ Hz

    Pd_range = np.arange(-35, 0, 0.1)
    error_list = []
    pd_list = []
    for pd_dbm in Pd_range:
        Pd = 10**((pd_dbm-30)/10)  # W
        A = R*Pd
        B = 2*e*(R*Pd+Id)
        C = (4*KT)/RL
        D = (R**2)*(Pd**2)*RIN
        E = 2*e*Id + C
        F = DR/np.sqrt(2)
        no_of_bits_hat = (
            1/6.02)*(20*np.log10(A/((np.sqrt(B+C+D)+np.sqrt(E))*np.sqrt(F)))-1.76)
        error = abs(no_of_bits_hat - no_of_bits)
        error_list.append(error)
        pd_list.append(pd_dbm)

    min_error = min(error_list)
    min_error_idx = error_list.index(min_error)
    pd = pd_list[min_error_idx]
    print("*******Calculated PD Sensitivity*****", pd)
    return pd
# """_summary_
# Modelling equation 13 from 'Scaling up silicon photonic based accelerators: challenges and opportunities'
# """


def findNandPLaser(pd):
    # Calculating the N value such that 10dBm- PLaser(Eq 13) = Min Error. This gives the N value for a given Bitwidth and data rate for MAM
    error_list = []
    N_list = []
    PLaser_list = []
    N_range = np.arange(1, 160, 1)
    for N in N_range:
        M = N
        Pd_dbm = pd  # dBm
        Pd = 10**((Pd_dbm)/10)  # mW

        # dB/mm #! 0.1 distance between MRM and MRR
        eta_WG = 10**((0.3*(0.02+0.1)*N)/10)
        eta_SMF = 1
        eta_EC = 10**(1.6/10)
        eta_WPE = 0.1
        IL_MRM = 10**(4/10)
        IL_MRR = 10**(0.01/10)
        OBL_MRM = 10**(0.01/10)
        OBL_MRR = 10**(0.01/10)
        # ! 2 dB additional network penalty for AMM
        IL_Penalty = 10**((4.8+1)/10)
        EL_splitter = 10**((0.01)/10)
        d_MRR = 0.02  # mm

        eta_EC = 1/eta_EC
        IL_MRM = 1/IL_MRM
        OBL_MRM = 1/OBL_MRM
        EL_splitter = 1/EL_splitter
        IL_MRR = 1/IL_MRR
        OBL_MRR = 1/OBL_MRR
        IL_Penalty = 1/IL_Penalty
        # print(eta_EC)
        # print(IL_MRM)
        # print(IL_MRR)
        # print(OBL_MRM)
        # print(IL_MRR)
        # print(IL_Penalty)
        # print(EL_splitter)

        A = M*eta_WG
        B = eta_SMF * eta_EC * IL_MRM * \
            (OBL_MRM**(N-1))*(EL_splitter**(np.log2(M)))
        C = Pd
        D = eta_WPE*IL_MRR*(OBL_MRR**(N-1))*IL_Penalty
        PLaser = (A/B)*(C/D)
        PLaser_dbm = 10*np.log10(PLaser)
        PLaser_dbm = PLaser_dbm/2
        # print('Plaser :', PLaser_dbm)
        error_list.append(abs(10-PLaser_dbm))
        N_list.append(N)
        PLaser_list.append(PLaser_dbm)

    min_error_idx = error_list.index(min(error_list))
    print("Minimum Error ", min(error_list))
    N = N_list[min_error_idx]
    PLaser = PLaser_list[min_error_idx]
    print("Calculated N value", N)
    return N, PLaser


def findOpticalRecievedPower(N, PLaser):
    # # Optical Power Calculation
    Psmf_att = 0
    Pec_il = 1.6
    Psi_att = 0.3
    Pmrm_ip_il = 4
    Pmrm_ip_obl = 0.01
    Psplitter_il = 0.01
    Pmrr_w_il = 0.01
    Pmrr_w_obl = 0.01
    p_penalty = 4.8
    N = N
    M = N
    Pout = PLaser - Psmf_att - Pec_il - Psi_att - Pmrm_ip_il - (N-1)*Pmrm_ip_obl - (
        10*np.log10(M)+Psplitter_il*np.log2(M)) - Pmrr_w_il - (N-1)*Pmrr_w_obl - p_penalty

    return Pout


# inputs parameters
bits_range = [1]
DR_range = [1, 3, 5, 10, 20 , 30, 40 , 50]
result_list = []
for no_of_bits in bits_range:
    for DR in DR_range:
        result = {}
        Pd = findPDSensitivity(no_of_bits, DR)
        N, PLaser = findNandPLaser(Pd)
        Pout = findOpticalRecievedPower(N, PLaser)
        # result['PD_Sensitivity'] = Pd
        result['N'] = N
        # result['PLaser'] = PLaser
        result['Pout'] = Pout
        result['no_of_bits'] = no_of_bits
        result['DR'] = DR
        result_list.append(result)
df = pd.DataFrame(result_list)
df.to_csv('Result/Revision/AMM_BNN_Recieved_Power2.csv')
