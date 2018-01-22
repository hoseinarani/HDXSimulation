import os, sys, re
from collections import defaultdict
import numpy as np
import glob
import sys
import math
cur_dir = os.getcwd()
constR = 1.987/1000. #kcal/mol*K
AA_3to1 = {'GLY':'G', 'ALA':'A', 'VAL':'V', 'LEU':'L', 'ILE':'I',
'MET':'M', 'PHE':'F', 'TRP':'W', 'PRO':'P', 'SER':'S',
'THR':'T', 'CYS':'C', 'TYR':'Y', 'ASN':'N', 'GLN':'Q',
'ASP':'D', 'GLU':'E', 'LYS':'K', 'ARG':'R', 'HIS':'H',
'NTR':'X', 'CTR':'Z'}
AA_1to3 = {'G':'GLY', 'A':'ALA', 'V':'VAL', 'L':'LEU', 'I':'ILE',
'M':'MET', 'F':'PHE', 'W':'TRP', 'P':'PRO', 'S':'SER',
'T':'THR', 'C':'CYS', 'Y':'TYR', 'N':'ASN', 'Q':'GLN',
'D':'ASP', 'E':'GLU', 'K':'LYS', 'R':'ARG', 'H':'HIS',
'X':'NTR', 'Z':'CTR'}

def sequence_kint(seq):
     timeUnit='min'
     myTempC=0
     pD_read = 7.4
     myTemp = float(myTempC) + 273.0
     pD_read = float(pD_read)
     if re.search('sec', timeUnit):
          timeUnitFactor = 1/60.
     else:
          timeUnitFactor = 1.0
     log_AB_LR_dic = {#residue: [AL, AR, BL, BR ],
     'ALA': [0.00, 0.00, 0.00, 0.00],
     'ARG': [-0.59, -0.32, 0.08, 0.22],
     'ASN': [-0.58, -0.13, 0.49, 0.32],
     'ASPd': [0.9, 0.58, -0.30, -0.18], #D(COO-)
     'ASPp': [-0.9, -0.12, 0.69, 0.6], #D(COOH)
     'CYS': [-0.54, -0.46, 0.62, 0.55],
     'GLY': [-0.22, 0.22, 0.27, 0.17],
     'GLN': [-0.47, -0.27, 0.06, 0.20],
     'GLUd': [-0.9, 0.31, -0.51, -0.15], #E(COO-)
     'GLUp': [-0.6, -0.27, 0.24, 0.39], #E(COOH)
     'HISd': ['NA', 'NA', -1.0, 0.14],
     'HISp': [-0.8, -0.51, 0.8, 0.83],
     'ILE': [-0.91, -0.59, -0.73, -0.23],
     'ILE': [-0.91, -0.59, -0.73, -0.23],
     'LEU': [-0.57, -0.13, -0.58, -0.21],
     'LYS': [-0.56, -0.29, -0.04, 0.12],
     'MET': [-0.64, -0.28, -0.01, 0.11],
     'PHE': [-0.52, -0.43, -0.24, 0.06],
     'PRO': ['NA', -0.19, 'NA', -0.24], #P-trans
      # 'PROc': ['NA', -0.85, 'NA', 0.60], #P-cis
     'SER': [-0.44, -0.39, 0.37, 0.30],
     'THR': [-0.79, -0.47, -0.07, 0.20],
     'TRP': [-0.40, -0.44, -0.41, -0.11],
     'TYR': [-0.41, -0.37, -0.27, 0.05],
     'VAL': [-0.74, -0.30, -0.70, -0.14],
     'NTR': ['NA', -1.32, 'NA', 1.62],
     'CTRd': [0.96, 'NA', -1.8, 'NA'],
     'CTRp': [0.05, 'NA', 'NA', 'NA'],
     }
     # print(log_AB_LR_dic['NTR'])
     #titrable_res_dic = {'ASP':4.08, 'GLU':4.53, 'HIS':7.02, 'CTR':7.0} # Bai et al. 1993
     titrable_res_dic = {'ASP':4.50, 'GLU':4.50, 'HIS':6.75, 'CTR':3.720} #, 'NTR':8.02} # SPHERE
     # activation energies (kcal/mol)
     #Ea_A = 14.0; Ea_B = 17.0; Ea_W = 19.0 # Bai 1993     
     Ea_A = 15.0; Ea_B = 2.6; Ea_W = 13.0 # SPHERE
     """
     Table I from Bai et al. (1993)
     """
     # standard reference rate constants at 293K
     ###########################################
     ## at high salt condition (0.5 KCl) oligo-peptide (from TableIII of Bai 1993)
     """
     kA_ref = (timeUnitFactor)*10**1.56 #1/(M*min)
     kB_ref = (timeUnitFactor)*10**10.20 #1/(M*min)
     kW_ref = (timeUnitFactor)*10**-2.3 #1/min
     """
     ## at high salt condition (0.5 KCl) poly-peptide (from TableIII of Bai 1993)
     """
     kA_ref = (timeUnitFactor)*10**1.19 #1/(M*min)
     kB_ref = (timeUnitFactor)*10**9.900 #1/(M*min)
     kW_ref = (timeUnitFactor)*10**-2.5 #1/min
     """
     # at normal lower salt condition oligo-peptide (from TableIII of Bai 1993)
     #kA_ref = (timeUnitFactor)*10.**2.04 #1/(M*min)
     #kB_ref = (timeUnitFactor)*10.**10.36 #1/(M*min)
     #kW_ref = (timeUnitFactor)*10.**-1.5 #1/min
     # logkA_ref = math.log10((timeUnitFactor)*10**2.04) #1/(M*min)
     # logkB_ref = math.log10((timeUnitFactor)*10**10.36) #1/(M*min)
     # logkW_ref = math.log10((timeUnitFactor)*10**-1.5) #1/min

     # at normal lower salt condition poly-peptide (from TableIII of Bai 1993)
     kA_ref = (timeUnitFactor)*10**1.62 #1/(M*min)
     kB_ref = (timeUnitFactor)*10**10.05 #1/(M*min)
     kW_ref = (timeUnitFactor)*10**-1.5 #1/min
     logkA_ref = math.log10(kA_ref)
     logkB_ref = math.log10(kB_ref)
     logkW_ref = math.log10(kW_ref)

     ###########################################
     # pK(D, 20C) = 15.65
     #pkD = 15.65
     pkD = 15.05 # SPHERE; Connelly 1993
     ##############################
     pD_corr =pD_read + 0.4
     pOD_corr = pkD - pD_corr
     ##############################
     #Dcation = 10**-(pD_corr) # pD = -log10[D+]
     #ODanion = 10**-(pOD_corr)
     k_int=[]
     seq_NH = 'x'.join(seq)
     for k, v in enumerate(seq_NH):
         logk_acid = 0; logk_base =0; logk_wat = 0; logk_int = 0
         if v == 'x':
                res_R = AA_1to3[seq_NH[k - 1]]
                res_L = AA_1to3[seq_NH[k + 1]]
                special_res = titrable_res_dic.keys()
                for s_res in special_res:
                        if res_R == s_res:
                            if pD_corr <= titrable_res_dic[s_res]:
                                res_R = res_R + 'p'
                            else:   
                                res_R = res_R + 'd'
                        if res_L == s_res:
                            if pD_corr <= titrable_res_dic[s_res]:
                                res_L = res_L + 'p'
                            else:
                                res_L = res_L + 'd'
                if log_AB_LR_dic[res_L][0] == 'NA':
                      #A_L = 1.0
                      logA_L = 0.0
                else:
                      #A_L = 10**(log_AB_LR_dic[res_L][0])
                      logA_L = log_AB_LR_dic[res_L][0]
                      #
                if log_AB_LR_dic[res_R][1] == 'NA':
                      #A_R = 1.0
                      logA_R = 0.0
                else:
                      #A_R = 10**(log_AB_LR_dic[res_R][1])
                      logA_R = log_AB_LR_dic[res_R][1]
                      #
                if log_AB_LR_dic[res_L][2] == 'NA':
                      #B_L = 1.0
                      logB_L = 0.0
                else:
                      #B_L = 10**(log_AB_LR_dic[res_L][2])
                      logB_L = log_AB_LR_dic[res_L][2]
                      #
                if log_AB_LR_dic[res_R][3] == 'NA':
                      #B_R = 1.0
                      logB_R = 0.0
                else:
                      #B_R = 10**(log_AB_LR_dic[res_R][3])
                      logB_R = log_AB_LR_dic[res_R][3]
                #kint = krc_293K * math.exp(-Ea*(1/float(myTemp)-1/293.0)/constR)
                TempFactor = (1.0/myTemp - 1.0/293.0)/constR
                #k_acid = (kA_ref * A_L * A_R * Dcation) * math.exp(-Ea_A*TempFactor)
                #k_base = (kB_ref * B_L * B_R * ODanion) * math.exp(-Ea_B*TempFactor)
                #k_wat = (kW_ref * B_L * B_R) * math.exp(-Ea_W*TempFactor)
                #k_int = (k_acid + k_base + k_wat)
                logk_acid = (logkA_ref + logA_L + logA_R - pD_corr)
                logk_base = (logkB_ref + logB_L + logB_R - pOD_corr)
                logk_wat = (logkW_ref + logB_L + logB_R)
                #
                k_acid = (10**(logk_acid)) * (math.exp(-Ea_A*TempFactor))
                k_base = (10**(logk_base)) * (math.exp(-Ea_B*TempFactor))
                k_wat = (10**(logk_wat)) * (math.exp(-Ea_W*TempFactor))
                k_int.append(k_acid + k_base + k_wat)
     return k_int

# pD_corr = 7.4
# myTempC = 0
# seq='VSQEEVKKWAESLENLINHECGLAAFKAFLKSEYSEENIDFWISCEEYKKIKS' \
#     'PSKLSPKAKKIYNEFISVQATKEVNLDSCTREETSRNMLEPTITCFDEAQKKIFNL' \
#     'MEKDSYRRFLKSRFYLDLT'
# timeUnit = 'sec'
# a=sequence_kint(seq)
# print(a)
# from calc_kint import sequence_kint
