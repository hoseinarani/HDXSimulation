# -*- coding: utf-8 -*-
"""
Created on Sun Apr  2 14:40:34 2017

@author: HM
"""

from math import *
import numpy as np
import itertools as it
import os
import matplotlib.pyplot as plt
from calc_kint import sequence_kint
import protein_info
os.getcwd()
pr=protein_info.get_protein_info()
import plotandsavedata
import time
from multiprocessing import pool

class model(object):
    def __init__(self):
        self.year=[]
        self.name = []
        self.par_name = []
        self.par_val = []
mod = model()
mod.name='persson_optimized'

mod.year='2015'
mod.par_name=['','']
mod.par_val=[]



path=os.getcwd()+'/models'
water = np.loadtxt( os.path.join(path,'final_water.dat'), skiprows=0)
polar = np.loadtxt( os.path.join(path,'final_polar.dat'), skiprows=0)

def calc_DI_MD(water, polar):
    w_ther=2.4298114
    p_ther=2.7297818
    # w_ther=2.6
    # p_ther=2.6
    water_tf = (water[:, 1:] < w_ther) * 1
    polar_tf = (polar[:, 1:] < p_ther) * 1
    resnumber = int(((np.size(water, 1) - 1) / 3))
    nf = int((np.size(water, 0)))
    water_number = np.zeros([nf, resnumber])
    polar_number = np.zeros([nf, resnumber])
    for i in range(1, resnumber + 1):
        water_number[:, i - 1] = np.sum(water_tf[:, i * 3 - 3:i * 3 - 1], axis=1).astype('float')
        polar_number[:, i - 1] = np.sum(polar_tf[:, i * 3 - 3:i * 3 - 1], axis=1).astype('float')
    open = ((water_number >= 2) * 1 & (polar_number <= 1) * 1)
    to = np.zeros(np.size(open, 1))
    tc = np.zeros(np.size(open, 1))
    for i in range(0, np.size(open, 1)):
        to[i] = open[:, i].sum()
        tc[i] = np.size(open, 0) - open[:, i].sum()
    PF = tc / (to+1e-6)
    indexofnan = [index for index, value in enumerate(PF) if np.isnan(value)]
    PF[indexofnan] = np.nanmax(PF)
    t = np.matrix([1, 3, 10, 30, 100, 300, 1000])
    kint = np.array(sequence_kint(pr.seq))
    #   PF=PF[1:]# remove the N-terminous residue from PF
    ind = []  # finding index of Prolins
    for k, l in enumerate(list(pr.seq)):
        if l == 'P':
            ind = k
            PF = np.insert(PF, ind, 'NAN')

    DI = np.zeros((7, pr.Cter_resid - pr.Nter_resid))
    for k, t in np.ndenumerate(t):
        DI[k[1], :] = (1 - np.exp(-kint / (PF + 1) * t))
    return DI

tic = time.clock()
DI_MD = calc_DI_MD(water, polar)
toc = time.clock()
print("Running time:", toc - tic)

plotandsavedata.func(DI_MD,pr,mod)

