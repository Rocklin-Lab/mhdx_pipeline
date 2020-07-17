import time
import random
import pickle
import Bio.PDB
import simanneal
import matplotlib
import importlib.util

import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.gridspec as gridspec

from Bio.PDB import *
from simanneal import Annealer
from scipy.stats import gaussian_kde
from matplotlib import pyplot as plt
from scipy.special import expit, logit
from scipy.signal import argrelextrema
from scipy.optimize import basinhopping
from sklearn.metrics import mean_squared_error
from scipy.optimize import fmin_powell, fmin_bfgs

molmass_spec = importlib.util.spec_from_file_location("molmass.py", "scripts/auxiliary/molmass.py")
molmass = importlib.util.module_from_spec(molmass_spec)
molmass_spec.loader.exec_module(molmass)


def distance_to_ss(x):
    out = []
    for i in range(len(x)):
        if x[i] == 'L':
            for j in range(1,len(x)):
                region = x[max([0, i-j]):min([len(x), i+j+1])]
                if 'E' in region or 'H' in region:
                    out.append(-j+1)
                    break
            if j == len(x) - 1: out.append(-j+1)
                
        else:
            for j in range(1,len(x)):
                region = x[max([0, i-j]):min([len(x), i+j+1])]
                if 'L' in region:
                    out.append(j)
                    break
            if j == len(x) - 1: out.append(j)
    return np.array(out)

def distance_to_nonpolar(x):
    out = []
    for i in range(len(x)):
        if x[i] in 'DEGHKNPQRST':
            for j in range(1,len(x)):
                region = x[max([0, i-j]):min([len(x), i+j+1])]
                for aa in 'ACFILMVWY':
                    if aa in region:
                        out.append(-j+1)
                        break
                if aa in region: break
            if j == len(x) - 1: out.append(-j+1)
                
        else:
            for j in range(1,len(x)):
                region = x[max([0, i-j]):min([len(x), i+j+1])]
                for aa in 'DEGHKNPQRST':
                    if aa in region:
                        out.append(j)
                        break
                if aa in region: break
            if j == len(x) - 1: out.append(j)
    return np.array(out)

def np_d(x, j):
    out = []
    for i in range(len(x)):
        out.append(0)
        region = x[max([0, i-j]):min([len(x), i+j+1])]
        for r in region:
            if r in 'ACFILMVWY':
                out[-1] += 1
    return np.array(out)




def PoiBinPdf(success_probabilities):

    number_trials = success_probabilities.size

    omega = 2 * np.pi / (number_trials + 1)

    chi = np.empty(number_trials + 1, dtype=complex)
    chi[0] = 1
    half_number_trials = int(
        number_trials / 2 + number_trials % 2)
    # set first half of chis:
    
    #idx_array = np.arange(1, half_number_trials + 1)
    exp_value = np.exp(omega * np.arange(1, half_number_trials + 1) * 1j)
    xy = 1 - success_probabilities + \
        success_probabilities * exp_value[:, np.newaxis]
    # sum over the principal values of the arguments of z:
    argz_sum = np.arctan2(xy.imag, xy.real).sum(axis=1)
    # get d value:
    #exparg = np.log(np.abs(xy)).sum(axis=1)
    d_value = np.exp(np.log(np.abs(xy)).sum(axis=1))
    # get chi values:
    chi[1:half_number_trials + 1] = d_value * np.exp(argz_sum * 1j)

    
    # set second half of chis:
    chi[half_number_trials + 1:number_trials + 1] = np.conjugate(
        chi[1:number_trials - half_number_trials + 1] [::-1])
    chi /= number_trials + 1
    xi = np.fft.fft(chi)
    return xi.real



Rconstant=0.0019872036
xs = np.logspace(-1,5,5000)
def hx(xs, k):
    return np.array([1.0-np.exp(-k*t) for t in xs])

def hx_of_t(t, k, backexchange=0.9, D2O_purity = 0.95, D2O_fraction = 0.9):
    return (1.0-np.exp(-k*t)) * (D2O_fraction * D2O_purity * backexchange) #D20 fraction * D20 purity * backexchange

def hx_of_t_alt(t, k, backexchange=0.001, D2O_purity = 0.95, D2O_fraction = 0.9):
    return (1.0-np.exp(-k*t)) * (np.exp(-k*backexchange)) * (D2O_fraction * D2O_purity) #D20 fraction * D20 purity * backexchange

def rate_with_fe(xs, rates, fes, backexchange=0.9, T=295): #295 C = 71.33 F):
    if type(fes) == type([]): fes == np.array(fes)
    fractions = np.exp(-fes / (Rconstant * T)) / (1.0 + np.exp(-fes / (Rconstant * T)))
    rates = rates * fractions
    sumrate=hx(xs,rates[0])
    for i in range(1,len(rates)):
        sumrate += hx(xs, rates[i])
        
    sumrate = sumrate * backexchange
    sumrate = sumrate * 0.9  #protein dilution into D2O
    sumrate = sumrate * 0.95 #D2O purity
    
    return sumrate

def rates_from_fes(p, fes, T):
    fes_temp=np.array(fes, dtype=np.float)
    fes_temp_indexing=np.logical_not(np.isnan(fes_temp))
    out=np.array(fes)
    
    fractions = np.exp(-fes_temp / (Rconstant * T)) / (1.0 + np.exp(-fes_temp / (Rconstant * T)))
    out[fes_temp_indexing] = p.rates[fes_temp_indexing] * fractions[fes_temp_indexing]
    return out

def fes_from_rates(p, meas_rates,unit='sec',fill_blank=False,fast_factor=2.0,min_fe=-10,use_approx_rates=False):
    meas_rates_temp=np.array(meas_rates, dtype=np.float)
    if fill_blank:
        max_rate=max(meas_rates_temp[np.isfinite(meas_rates_temp)])
        meas_rates_temp[np.isnan(meas_rates_temp)] = max_rate * fast_factor
    
    if unit == 'min':
        meas_rates_temp = meas_rates_temp / 60.0
    

    meas_rates_indexing=np.logical_not(np.isnan(meas_rates_temp))
    meas_rates_temp=meas_rates_temp[meas_rates_indexing]
    
    intrinsic_rates = p.rates if not use_approx_rates else p.approx_rates

    Kop = meas_rates_temp / (intrinsic_rates[meas_rates_indexing] - meas_rates_temp)
    
    out = np.array(meas_rates, dtype=np.float)
    fes = -Rconstant * p.temp * np.log(Kop)
    fes[np.isnan(fes)] = min_fe
    fes[fes < min_fe] = min_fe
    out[meas_rates_indexing] = fes
    out[intrinsic_rates == 0] = np.nan
    
    return out

def fe_from_rate(intrinsic_rate, meas_rate, temp, min_fe=-10):
    if meas_rate >= intrinsic_rate: return min_fe
    
    Kop = meas_rate / (intrinsic_rate - meas_rate)

    out= -Rconstant * temp * np.log(Kop)
    return max([out, min_fe])

def poibin_dist(t, rates, fes, backexchange=0.9, T=295):
    if type(fes) == type([]): fes = np.array(fes)
    fractions = np.exp(-fes / (Rconstant * T)) / (1.0 + np.exp(-fes / (Rconstant * T)))
    probabs = hx_of_t(t, rates*fractions, backexchange=backexchange)
    probabs[probabs != probabs] = 0
    #print rates, fes, fractions, probabs
    return PoiBinPdf(probabs)

def poibin_dist_rates(t, rates, backexchange=0.9, T=295):
    #if type(rates) == type([]): rates = np.array(rates)
    probabs = hx_of_t(t, np.exp(rates), backexchange=backexchange)
    #print rates, fes, fractions, probabs
    return PoiBinPdf(probabs)

def isotope_poibin_dist(isotope_dist, t, rates, fes, nbins = 50, backexchange=0.9):
    pbdist=poibin_dist(t, rates, fes, backexchange=backexchange)
    iso_binom=np.convolve(pbdist, isotope_dist)
    return iso_binom[0:nbins] / max(iso_binom[0:nbins])

def isotope_poibin_dist_rates(isotope_dist, t, rates, nbins = 50, backexchange=0.9):
    pbdist = poibin_dist_rates(t, rates, backexchange=backexchange)
    iso_binom=np.convolve(pbdist, isotope_dist)
    return iso_binom[0:nbins] / max(iso_binom[0:nbins])

def isotope_poibin_dist_rates_at_timepoints(timepoints, rates, nbins, isotope_dist, backexchange):
    out=np.zeros((len(timepoints), nbins))
    for i, timepoint in enumerate(timepoints):
        out[i,:] = isotope_poibin_dist_rates(isotope_dist,
                                           timepoint,
                                           rates,
                                           nbins = nbins,
                                           backexchange=backexchange[i])
    return np.ravel(out)

def get_concat_data(timepoints, data):
    if len(np.ravel(timepoints)) == 1:
        timepoints, data = [timepoints], [data]

    concat_data=[]
    for timepoint, dataset in zip(timepoints, data):
        concat_data  += list(dataset['major_species_integrated_intensities'] / max(dataset['major_species_integrated_intensities']))
    return np.array(concat_data)

def rmse(predictions, targets):
    #print np.sqrt(((predictions - targets) ** 2).mean())
    return np.sqrt(((predictions - targets) ** 2).mean())

def smooth_dist_cutoff(x):
    return 1.0 / (1.0 + (np.exp(x - 10)))

def phi_psi_omega_to_abego(phi, psi, omega):
    if psi == None: return 'O'
    if omega == None: omega = 180
    if phi == None: phi=90
    phi = 180 * phi / np.pi
    psi = 180 * psi / np.pi
    omega = 180 * omega / np.pi




    if abs(omega) < 90:
        return 'O'
    elif phi > 0:
        if -100.0 <= psi < 100:
            return 'G'
        else:
            return 'E'
    else:
        if -75.0 <= psi < 50:
            return 'A'
        else:
            return 'B'
    return 'X'

def abego_string(structure):
    model = structure[0]

    residues = [res for res in model.get_residues()]
    phi_psi_omega=[]
    for i in range(len(residues)):
        if i +1 == len(residues):
            omega = None
            psi = None
        else:
            if i > 0: aminus1 = residues[i-1]['C'].get_vector()
            a0 = residues[i]['N'].get_vector()
            a1 = residues[i]['CA'].get_vector()
            a2 = residues[i]['C'].get_vector()
            a3 = residues[i+1]['N'].get_vector()
            a4 = residues[i+1]['CA'].get_vector()
            omega = Bio.PDB.calc_dihedral(a1,a2,a3,a4)
            psi = Bio.PDB.calc_dihedral(a0, a1, a2, a3)
        if i > 0:
            phi = Bio.PDB.calc_dihedral(aminus1, a0, a1, a2)
        else:
            phi = None
        phi_psi_omega.append((phi, psi, omega))

    out = ''
    for x in phi_psi_omega:
        out += phi_psi_omega_to_abego(x[0], x[1], x[2])
    return out

def dssp(prot):
    abego_string = ('G' * len(prot.Nterm)) +prot.abego_string
    hbond_partners = prot.hbond_partners
 
    out = 'L'

    for i in range(1,len(abego_string)-1):
        if abego_string[i] in 'GO':
            out += 'L'
        elif abego_string[i] == 'A':
            nextaa='L'
            for partner in [-5, -4, -3, 3, 4, 5]:
                if i+partner < 0 or i+partner >= len(abego_string): continue
                if partner < 0:
                    if hbond_partners[i] == i+partner:
                        seg=abego_string[i+partner:i+1]
                        if seg.count('A') == len(seg):
                            nextaa = 'H'
                else:
                    if hbond_partners[i+partner] == i:
                        seg=abego_string[i:i+partner+1]
                        if seg.count('A') == len(seg):
                            nextaa = 'H'
            out += nextaa
        elif abego_string[i] in 'BE':
            nextaa='L'
            if hbond_partners[i] != None and i in hbond_partners:
                if abego_string[hbond_partners[i]] in 'BE' and abego_string[hbond_partners.index(i)] in 'BE':
                    nextaa='E'

            out += nextaa
    out += 'L'
    
    newE = []
    for i in range(1,len(abego_string)-1):
        if abego_string[i] in 'BE' and out[i] == 'L':
            if out[i-1] == 'E' and out[i+1] == 'E':
                newE.append(i)
            
            elif 'E' in [out[i-1], out[i+1]]:
                if hbond_partners[i] != None:
                    if abego_string[hbond_partners[i]] in 'BE':
                        newE.append(i)
                elif i in hbond_partners:
                    if abego_string[hbond_partners.index(i)] in 'BE':
                        newE.append(i)
    for i in newE:
        out = out[0:i] + 'E' + out[i+1:]

    #for i in range(len(abego_string)):
    #    print i, abego_string[i], hbond_partners[i], out[i], prot.sequence[i]
                
    
    
    return out

def plot_decay(p, xs, fes=None, rates=None, color=None, backexchange=0.9,label=''):
    if fes == None: fes = np.array([-5 for x in p.sequence])
    if rates == None: rates=p.rates
    plt.semilogx(xs, rate_with_fe(xs, rates, fes,backexchange=backexchange),color=color,label=label)

class hxprot:
    def calc_pred_kc(self, pD, temp):
        bai_data={}
        #pD = 5.93
        #Temp = 298

        bai_data['ka'], bai_data['kb'], bai_data['kw'] = 6.95E-01, 1.87E+08, 5.27E-04 
        bai_data['R']=1.987
        bai_data['pKc_D'] = np.log10( (10**-4.48) * np.exp(-1000 * (((1.0 / temp) - (1.0 / 278)) / bai_data['R'])))*-1
        bai_data['pKc_E'] = np.log10( (10**-4.93) * np.exp(-1083 * (((1.0 / temp) - (1.0 / 278)) / bai_data['R'])))*-1
        bai_data['pKc_H'] = np.log10( (10**-7.42) * np.exp(-7500 * (((1.0 / temp) - (1.0 / 278)) / bai_data['R'])))*-1


        
        bai_raw_data="""A   0.000000000 0.000000000     0.000000000 0.000000000
        C   -0.540000000    -0.460000000        0.620000000 0.550000000
        C2  -0.740000000    -0.580000000        0.550000000 0.460000000
        D   0.886779836 0.569277924     0.136651587 -0.118131722
        D+  -0.900000000    -0.120000000        0.690000000 0.600000000
        E   -0.866486913    0.283403579     -0.068668145    -0.071224325
        E+  -0.600000000    -0.270000000        0.240000000 0.390000000
        F   -0.520000000    -0.430000000        -0.235859464    0.063131587
        G   -0.220000000    0.218176047     0.267251569 0.170000000
        H   -0.655259338    -0.443090948        0.770757577 0.803458167
        I   -0.910000000    -0.590000000        -0.730000000    -0.230000000
        K   -0.560000000    -0.290000000        -0.040000000    0.120000000
        L   -0.570000000    -0.130000000        -0.576252728    -0.210000000
        M   -0.640000000    -0.280000000        -0.008954843    0.110000000
        N   -0.580000000    -0.130000000        0.490000000 0.320000000
        P   0.000000000 -0.194773472        0.000000000 -0.240000000
        Pc  0.000000000 -0.854416534        0.000000000 0.600000000
        Q   -0.470000000    -0.270000000        0.060000000 0.200000000
        R   -0.590000000    -0.320000000        0.076712254 0.220000000
        S   -0.437992278    -0.388518935        0.370000000 0.299550286
        T   -0.790000000    -0.468073126        -0.066257980    0.200000000
        V   -0.739022273    -0.300000000        -0.701934483    -0.140000000
        W   -0.400000000    -0.440000000        -0.410000000    -0.110000000
        Y   -0.410000000    -0.370000000        -0.270000000    0.050000000
        NT  0.000000000 -1.320000000        0.000000000 1.620000000
        CT  0.928161740 0.000000000     -1.800000000    0.000000000
        NMe 0.511939295 0.000000000     -0.577243551    0.000000000
        Ac  0.000000000 0.293000000     0.000000000 -0.197000000""".split('\n')

        bai_data['acid_left'] = {x.split()[0]: float(x.split()[1]) for x in bai_raw_data}
        bai_data['acid_right'] = {x.split()[0]: float(x.split()[2]) for x in bai_raw_data}
        bai_data['base_left'] = {x.split()[0]: float(x.split()[3]) for x in bai_raw_data}
        bai_data['base_right'] = {x.split()[0]: float(x.split()[4]) for x in bai_raw_data}

        bai_data['acid_left']['D'] = np.log10(  10**(-0.9   -  pD  )  /(10**(- bai_data['pKc_D']  )+10**(-  pD  ))+10**(0.9-  bai_data['pKc_D']  )/(10**(-  bai_data['pKc_D']  )+10**(-  pD)))
        bai_data['acid_left']['E'] = np.log10(  10**(-0.6   -  pD  )  /(10**(- bai_data['pKc_E']  )+10**(-  pD  ))+10**(-0.9-  bai_data['pKc_E']  )/(10**(-  bai_data['pKc_E']  )+10**(-  pD)))
        bai_data['acid_left']['H'] = np.log10(  10**(-0.8   -  pD  )  /(10**(- bai_data['pKc_H']  )+10**(-  pD  ))+10**(   0-  bai_data['pKc_H']  )/(10**(-  bai_data['pKc_H']  )+10**(-  pD)))
        bai_data['acid_left']['CT'] = np.log10(  10**(-0.05 -  pD  )  /(10**(- bai_data['pKc_E']  )+10**(-  pD  ))+10**(0.96 -  bai_data['pKc_E']  )/(10**(-  bai_data['pKc_E']  )+10**(-  pD)))
        
        bai_data['acid_right']['D'] =np.log10(   10**(-0.12-   pD  ) / (10**(- bai_data['pKc_D']  )+10**(-  pD  ))+10**(0.58-  bai_data['pKc_D']   )/(10**(-  bai_data['pKc_D']  )+10**(-  pD)))
        bai_data['acid_right']['E'] =np.log10(   10**(-0.27-   pD  ) / (10**(- bai_data['pKc_E']  )+10**(-  pD  ))+10**(0.31-  bai_data['pKc_E']   )/(10**(-  bai_data['pKc_E']  )+10**(-  pD)))
        bai_data['acid_right']['H'] =np.log10(   10**(-0.51-   pD  ) / (10**(- bai_data['pKc_H']  )+10**(-  pD  ))+10**(   0-  bai_data['pKc_H']   )/(10**(-  bai_data['pKc_H']  )+10**(-  pD)))
        
        bai_data['base_left']['D'] = np.log10(   10**(0.69-    pD  ) / (10**(-  bai_data['pKc_D'] )+10**(-  pD  ))+10**( 0.10 -  bai_data['pKc_D']   )/(10**(-  bai_data['pKc_D']  )+10**(-  pD)))
        bai_data['base_left']['E'] = np.log10(   10**(0.24-    pD  ) / (10**(-  bai_data['pKc_E'] )+10**(-  pD  ))+10**(-0.11 -  bai_data['pKc_E']   )/(10**(-  bai_data['pKc_E']  )+10**(-  pD)))
        bai_data['base_left']['H'] = np.log10(   10**(0.80-    pD  ) / (10**(-  bai_data['pKc_H'] )+10**(-  pD  ))+10**(-0.10 -  bai_data['pKc_H']   )/(10**(-  bai_data['pKc_H']  )+10**(-  pD)))
        
        bai_data['base_right']['D'] = np.log10(   10**(0.60-    pD  ) / (10**(-  bai_data['pKc_D'] )+10**(-  pD  ))+10**(-0.18 -  bai_data['pKc_D']   )/(10**(-  bai_data['pKc_D']  )+10**(-  pD)))
        bai_data['base_right']['E'] = np.log10(   10**(0.39-    pD  ) / (10**(-  bai_data['pKc_E'] )+10**(-  pD  ))+10**(-0.15 -  bai_data['pKc_E']   )/(10**(-  bai_data['pKc_E']  )+10**(-  pD)))
        bai_data['base_right']['H'] = np.log10(   10**(0.83-    pD  ) / (10**(-  bai_data['pKc_H'] )+10**(-  pD  ))+10**( 0.14 -  bai_data['pKc_H']   )/(10**(-  bai_data['pKc_H']  )+10**(-  pD)))
        
        
        seq=self.sequence
        bai_data['conc_Dplus'] = 10**(-1*pD)
        bai_data['conc_ODminus'] =10**(pD - 15.05) #pK(D,20c)
        bai_data['Fta']= np.exp(-1*14000*((1.0/temp)-(1.0/293))/bai_data['R'])
        bai_data['Ftb']= np.exp(-1*17000*((1.0/temp)-(1.0/293))/bai_data['R'])
        bai_data['Ftw']= np.exp(-1*19000*((1.0/temp)-(1.0/293))/bai_data['R'])
        
        pred_kc=[]
        for i, aa in enumerate(self.sequence):
            if i ==0 or aa in ['P']:
                pred_kc.append(0.0)
            else:
                Fa = 10.0 ** (bai_data['acid_left'][aa] + bai_data['acid_right'][seq[i-1]] )
                Fb = 10.0 ** (bai_data['base_left'][aa] + bai_data['base_right'][seq[i-1]] )
                if i == 1:
                    Fa *= 10.0 ** bai_data['acid_right']['NT']
                    Fb *= 10.0 ** bai_data['base_right']['NT']
                if i == len(seq) - 1:
                    Fa *= 10.0 ** bai_data['acid_left']['CT']
                    Fb *= 10.0 ** bai_data['base_left']['CT']

                Fa_ka_D_Fta = Fa * bai_data['ka'] * bai_data['conc_Dplus'] * bai_data['Fta']
                Fb_kb_OD_Ftb = Fb * bai_data['kb'] * bai_data['conc_ODminus'] * bai_data['Ftb']
                Fb_kw_Ftw = Fb * bai_data['kw'] * bai_data['Ftw']
                pred_kc.append(Fa_ka_D_Fta + Fb_kb_OD_Ftb + Fb_kw_Ftw)


        return np.array(pred_kc)

    def calc_fe_from_rate_grid(self, meas_rates):
        meas_rates=np.array(sorted(meas_rates))
        out=np.zeros((len(meas_rates), len(meas_rates)))
        #out[pos_index][rate_index] == fe
        int_rates=self.rates[self.allowed]
        for pos_index, int_rate in enumerate(int_rates):
            for rate_index, meas_rate in enumerate(meas_rates):
                out[pos_index][rate_index] = fe_from_rate(int_rate, meas_rate, self.temp)
        return out
    
    def calc_approx_fe_distribution(self, meas_rates):
        meas_rates=np.array(sorted(meas_rates))
        out=np.zeros((len(meas_rates), len(meas_rates)))
        #out[pos_index][rate_index] == fe
        int_rates=self.rates[self.allowed]



        for pos_index, int_rate in enumerate(int_rates):
            for rate_index, meas_rate in enumerate(meas_rates):
                out[pos_index][rate_index] = fe_from_rate(int_rate, meas_rate, self.temp)

        return out

    def __init__(self, pdbfilename='', seq='', Nterm='', Cterm='', pD=5.93, temp=295, hbond_length=2.7, hbond_angle=120):
        
        self.Nterm = Nterm
        self.Cterm = Cterm
        self.pD = pD
        self.temp = temp

        if pdbfilename != '':
            structure = Bio.PDB.PDBParser(QUIET=True).get_structure(pdbfilename,pdbfilename)
            model = structure[0]
            polypeptides = Bio.PDB.PPBuilder().build_peptides(model)
        
            
            residues = [res for res in model.get_residues()]
            seq = ''.join([Polypeptide.three_to_one(aa.get_resname()) for aa in residues])
        
        self.sequence    =   Nterm + seq + Cterm
        self.allowed     =   np.array([False for x in Nterm] + [False if x == 'P' else True for x in seq] + [False for x in Cterm])
        self.allowed[0:2]  = [False, False]
        self.rates       =   self.calc_pred_kc(pD, temp) #exchange per second
        self.rates[1]    = 0 #the second residue back-exchanges too rapidly to measure
        self.h_indices   = np.array([x for x in range(len(self.rates)) if self.rates[x] > 0])
        self.blank_fes = np.zeros(len(self.sequence))
        self.distance_to_nonpolar = distance_to_nonpolar(self.sequence)
        self.active_index=np.array([sum(self.allowed[0:x]) for x in range(len(self.allowed))], dtype=object)
        for i in range(len(self.active_index)-1):
            if self.active_index[i+1] == self.active_index[i]:
                self.active_index[i] = None
        self.idx_to_full=np.array([list(self.active_index).index(i) for i in range(sum(self.allowed))])
        

        if pdbfilename != '':
            hbonds=[]
            hbond_dists=[None for x in residues]
            hbond_angles=[None for x in residues]
            hbond_partners=[None for x in residues]
            all_hbonds=[]
            sc_hbonds=[]
            atoms_in_radius=[0 for x in Nterm]
            self.pairlist=[]
            for i, res in enumerate(residues):
                n=0
                atoms_in_radius.append(sum([smooth_dist_cutoff(np.linalg.norm(res['N'].coord.astype(float) - j.coord.astype(float)))  
                                for j in model.get_atoms() if 'H' not in j.get_name()]))
                  
                if 'H' not in res:
                    continue
                
                for j, res2 in enumerate(residues):
                    if j == i: continue
                    #print i+1, j+1,
                    new_hbond_dist = np.linalg.norm(res['H'].coord.astype(float) - res2['O'].coord.astype(float))
                    new_hbond_angle = np.rad2deg(min([Bio.PDB.calc_angle(res['N'].get_vector(),
                                                                           res['H'].get_vector(),
                                                                           res2['O'].get_vector()),
                                                   Bio.PDB.calc_angle(res['H'].get_vector(),
                                                                      res2['O'].get_vector(),
                                                                      res2['C'].get_vector())]))
                    if (new_hbond_dist < hbond_length) and (new_hbond_angle > hbond_angle):
                        hbonds.append(i)
                        all_hbonds.append(i)
                        all_hbonds.append(j)
                        hbond_partners[i] = j+len(Nterm)
                        
                        if hbond_dists [i] != None:
                            if new_hbond_dist < hbond_dists[i]:
                                hbond_dists[i] = new_hbond_dist
                                hbond_angles[i] = new_hbond_angle
                        
                        #print i, j, res['H'] - res2['O']
                        #print 'true'
                    for O_atom in [x for x in res2 if x.name[0] == 'O' and x.name != 'O']:
                        if np.linalg.norm(res['H'].coord.astype(float) - O_atom.coord.astype(float)) < hbond_length:
                            sc_hbonds.append(i) 
                    
                    if i < j:
                        if 'H' in res and 'H' in res2:
                            d = np.linalg.norm(res['H'].coord.astype(float) - res2['H'].coord.astype(float))
                            if d < 13:
                                self.pairlist.append((i+len(self.Nterm), j+len(self.Nterm), d))
                    
                    
            atoms_in_radius += [0 for x in Cterm]        
                    
            hbonds=sorted(set(hbonds))
            all_hbonds=sorted(set(all_hbonds))

            hbond_depth=[1 if x in all_hbonds else 0 for x in range(len(residues))]
            for i in range(1, len(residues)-1):
                if hbond_depth[i-1:i+2] == [1,0,1]:
                    hbond_depth[i] = 1
                    all_hbonds.append(i)

            n = 1
            while hbond_depth.count(n) > 1:
                for i in range(n, len(residues)-n):
                    #print i
                    if hbond_depth[i-1:i+2] in [[n, n, n], [n+1, n, n]] :
                        hbond_depth[i] += 1
                n+=1

            hbond_depth=np.array(hbond_depth) * np.array([1 if x in hbonds else 0 for x in range(len(residues))])

            



            self.abego_string = abego_string(structure)
            self.burial      =   np.array(atoms_in_radius)
            self.hbond_depth =   np.array([0 for x in Nterm] + list(hbond_depth) + [0 for x in Cterm])
            self.hbonds      =   np.array([0.0 if x == 0 else 1 for x in self.hbond_depth])
            self.hbond_dists =   np.array([None for x in Nterm] + list(hbond_dists))
            self.hbond_angles =   np.array([None for x in Nterm] + list(hbond_angles))
            self.hbond_partners = [None for x in Nterm] + hbond_partners + [None for x in Cterm]
            

            self.sc_hbonds   =   np.array([0 for x in Nterm] + [1 if x in sc_hbonds else 0 for x in range(len(residues))] + [0 for x in Cterm])
            self.hbond_fes = np.array([-1.0 if i == 0 else 2.0 for i in self.hbonds])
            self.dssp = dssp(self)
            self.distance_to_ss = distance_to_ss(self.dssp)
            self.approx_rates = self.calculate_approx_rates()


    def calculate_approx_rates(self):
        hbond_rates = self.rates[(self.rates != 0) & (self.hbonds == 1)]
        nohbond_rates = self.rates[(self.rates != 0) & (self.hbonds == 0)]
        return np.array([np.median(hbond_rates) for x in hbond_rates] + [np.median(nohbond_rates) for x in nohbond_rates])


    def random_sa_state(self):
        out=list(range(sum(self.allowed)))
        random.shuffle(out)
        return out
    
    def calculate_isotope_dist(self, undeut, n_isotopes=None, use_empirical=False):
        theo_isotope_dist = self.calculate_theo_isotope_dist(n_isotopes=n_isotopes)
        emp_isotope_dist = self.calculate_empirical_isotope_dist(undeut, n_isotopes=n_isotopes)

        minlen = min([len(theo_isotope_dist), len(emp_isotope_dist)])

        self.isotope_dist = emp_isotope_dist if use_empirical else theo_isotope_dist
        return np.linalg.norm(np.dot(theo_isotope_dist[0:minlen], emp_isotope_dist[0:minlen])) / np.linalg.norm(theo_isotope_dist) / np.linalg.norm(emp_isotope_dist)

    def calculate_empirical_isotope_dist(self, undeut, n_isotopes=None):
        isotope_dist=undeut['major_species_integrated_intensities'] / max(undeut['major_species_integrated_intensities'])
        isotope_dist= isotope_dist / max(isotope_dist)
        self.isotope_dist = isotope_dist if n_isotopes == None else isotope_dist[0:n_isotopes]
        return self.isotope_dist

    def calculate_theo_isotope_dist(self, n_isotopes=None):
        f = molmass.Formula(self.sequence)
        isotope_dist=np.array([x[1] for x in f.spectrum().values()])
        isotope_dist= isotope_dist / max(isotope_dist)
        self.isotope_dist = isotope_dist if n_isotopes == None else isotope_dist[0:n_isotopes]
        return self.isotope_dist

    def calculate_backexchange(self, fulldeut): 
        out= fmin_powell(lambda bx: rmse(fulldeut['major_species_integrated_intensities'] / max(fulldeut['major_species_integrated_intensities']),
                                               isotope_poibin_dist(self.isotope_dist,
                                                                   1e9,
                                                                   self.rates,
                                                                   self.blank_fes,
                                                                   nbins = len(fulldeut['major_species_integrated_intensities']),
                                                                   backexchange=expit(bx))),
                                2, disp=False)
        self.backexchange=expit(out)
        return self.backexchange
    
    def calculate_backexchange_time(self, fulldeut): 
        out= fmin_powell(lambda bx: rmse(fulldeut['major_species_integrated_intensities'] / max(fulldeut['major_species_integrated_intensities']),
                                               isotope_poibin_dist(self.isotope_dist,
                                                                   1e9,
                                                                   self.rates,
                                                                   self.blank_fes,
                                                                   nbins = len(fulldeut['major_species_integrated_intensities']),
                                                                   backexchange=np.exp(bx))),
                                2, disp=False)
        self.backexchange=np.exp(out)
        return self.backexchange

    def data_rmse(self, timepoints, data, fes):
        if len(np.ravel(timepoints)) == 1:
            timepoints, data = [timepoints], [data]
        concat_data=[]
        concat_model=[]
        for timepoint, dataset in zip(timepoints, data):
            concat_data  += list(dataset['major_species_integrated_intensities'] / max(dataset['major_species_integrated_intensities']))
            concat_model += list(isotope_poibin_dist(self.isotope_dist,
                                               timepoint,
                                               self.rates,
                                               fes,
                                               nbins = len(dataset['major_species_integrated_intensities']),
                                               backexchange=self.backexchange))
        concat_data=np.array(concat_data)
        concat_model=np.array(concat_model)
        return rmse(concat_data[concat_data > 0], concat_model[concat_data > 0])

    def data_rate_rmse(self, timepoints, rates, concat_data, nbins, backexchange):
        concat_model=isotope_poibin_dist_rates_at_timepoints(timepoints, rates, nbins, self.isotope_dist, backexchange)
        return mean_squared_error(concat_data[concat_data > 0], concat_model[concat_data > 0])
        
    
    def fit_data_scalar(self, timepoints, data, init_fes):
        scale= fmin_powell(lambda scale: self.data_rmse(timepoints,
                                                      data,
                                                      np.where(init_fes == -1, -1, init_fes * scale[0])),
                                1, disp=False)
        return np.where(init_fes == -1, -1, init_fes * scale)
    
    def fit_data(self, timepoints, data, init_fes, scalar_init = False, corr_values=None, corr_weight=None):
        if scalar_init: init_fes = self.fit_data_scalar(timepoints, data, init_fes)
        
        if corr_values == None or corr_weight == 0:
            opt= fmin_powell(lambda fes: self.data_rmse(timepoints,
                                                          data,
                                                          np.where(self.allowed == 1, np.maximum(fes, -4), -4)) + max(fes)/5000.0,
                             init_fes, disp=False)
        else:
            #opt= fmin_powell(lambda fes: (self.data_rmse(timepoints,
            #                                              data,
            #                                              np.where(self.allowed == 1, np.maximum(fes, -4), -4)) +
            #                                  (corr_weight * (1.0-(np.corrcoef(fes, corr_values)[0][1])))),
            #                 init_fes, disp=False)

            opt= fmin_powell(lambda fes: (self.data_rmse(timepoints,
                                                          data,
                                                          np.where(self.allowed == 1, np.maximum(fes, -4), -4)) +
                                              (corr_weight *        (np.average( ((fes[0:-1] - fes[1:]) * np.maximum(corr_values[0:-1] * corr_values[1:], 0.2))**2.0)**0.5) )+
                                              max(fes)/5000.0 ),

                             init_fes, disp=False)#,maxfun=10)

            
        return np.where(self.allowed == 1, np.maximum(opt, -4), -4)
        
    def fit_rates(self, timepoints, init_rates=[], fixed_rates=[], niter=200, T=0.00003, stepsize=0.02, data=None, concat_data=None, nbins=None):
        if len(np.ravel(timepoints)) == 1:
            timepoints, data = [timepoints], [data]
        
        if data != None:
            concat_data=get_concat_data(timepoints, data)
            nbins = len(data[0]['major_species_integrated_intensities'])
        
        if len(init_rates) == 0:
            init_rates = np.linspace(1,-12,len([x for x in self.rates if x != 0]))

        fitbx = np.array(self.backexchange)
        if len(fitbx) != len(timepoints):
            fitbx = np.reshape(fitbx, (len(timepoints),))

        out=basinhopping(lambda rates: self.data_rate_rmse(timepoints, np.concatenate((rates, fixed_rates)), concat_data, nbins, fitbx),
                            init_rates, 
                            niter=niter, #ideally 1000
                            T=T,
                            stepsize=stepsize,
                            minimizer_kwargs={'options': {'maxiter': 1}})
        return out

    def fit_rates_and_backexchange(self, timepoints, init_rates = [], init_bx_rates=[], fixed_rates=[], niter=200, T=0.00003, stepsize=0.02, data=None, concat_data=None, nbins=None):
        if len(np.ravel(timepoints)) == 1:
            timepoints, data = [timepoints], [data]
        
        if data != None:
            concat_data=get_concat_data(timepoints, data)
            nbins = len(data[0]['major_species_integrated_intensities'])
        
        if len(init_bx_rates) == 0 and len(init_rates) == 0:
            init_bx_rates = np.concatenate(([logit(self.backexchange)], np.linspace(1,-12,len([x for x in self.rates if x != 0]))))
        elif init_rates != 0:
            init_bx_rates = np.concatenate(([logit(self.backexchange)], init_rates))

        out=basinhopping(lambda bx_rates: self.data_rate_rmse(timepoints, np.concatenate((bx_rates[1:], fixed_rates)), concat_data, nbins, expit(bx_rates[0])),
                            init_bx_rates, 
                            niter=niter, #ideally 1000
                            T=T,
                            stepsize=stepsize,
                            minimizer_kwargs={'options': {'maxiter': 1}})
        return out
    
    def fit_rates_withfast(self, timepoints, init_rates, niter=200, T=0.00003, stepsize=0.02, data=None, concat_data=None, nbins=None):
        if len(np.ravel(timepoints)) == 1:
            timepoints, data = [timepoints], [data]
        
        if data != None:
            concat_data=get_concat_data(timepoints, data)
            nbins = len(data[0]['major_species_integrated_intensities'])

        
        fast_rates = np.array([5 for x in range(len(self.Nterm + self.Cterm)-2)])

        out=basinhopping(lambda rates: self.data_rate_rmse(timepoints, np.concatenate((rates, fast_rates)), concat_data, nbins, self.backexchange),
                            init_rates,
                            niter=niter, #ideally 1000
                            T=T,
                            stepsize=stepsize,
                            minimizer_kwargs={'options': {'maxiter': 1}})
        return out, np.array(sorted(np.concatenate((out.x, fast_rates))))
        
def make_plots(data, p, fit, title, fn):
    timelabels=['5s','13s','20s','60s','2.5m','7.5m','21m','50m','1hr55m','3hr51m','6hr25m','8hr30m','1hr30m 95C']
    plt.figure(figsize=(8.5, 11))
    
    gs = gridspec.GridSpec(14, 2)
    ax = plt.subplot(gs[0, 1])
    plt.bar(range(50), data['integrated_intensities'][0]/max(data['integrated_intensities'][0]),
        alpha=0.2, color='blue')
    plt.text(0.99,0.98,'UN',transform=ax.transAxes,horizontalalignment='right',verticalalignment='top',fontsize=8)
    plt.yticks([])
    plt.ylim(0, 1.1)
    plt.xticks([])
    
    allmodel=[]
    for i in range(2, 15):
        ax = plt.subplot(gs[i-1, 1])
        plt.bar(range(50), data['major_species_integrated_intensities'][i]/max(data['major_species_integrated_intensities'][i]),
                alpha=0.3, color='blue')
        
        mytp = float(timepoints[i-1]) if i < 14 else 1e12
        fithist = isotope_poibin_dist(p.isotope_dist,
                                               mytp,
                                               p.rates,
                                               fit,
                                               nbins = len(data['integrated_intensities'][i]),
                                               backexchange=p.backexchange)
        
        
        
        if i != 14: allmodel.append(fithist)
        
        plt.bar(range(50), fithist,alpha=0.4,color='red',fill='red',edgecolor='red')
        plt.yticks([])
        plt.ylim(0, 1.1)
        plt.xticks([])
        plt.text(0.01,0.98,timelabels[i-2],transform=ax.transAxes,horizontalalignment='left',verticalalignment='top',fontsize=8)
    plt.xticks([0,10,20,30,40,50])
        
    sns.despine()
    
    ax = plt.subplot(gs[0:3, 0])
    
    temp_fit=np.copy(fit)
    temp_fit[temp_fit != temp_fit] = 100
    plot_decay(p, xs, fes=None, color='black', backexchange=p.backexchange,label='')
    plot_decay(p, xs, fes=temp_fit, color='red', backexchange=p.backexchange,label='')
    
    deuts=[x - data['major_species_centroid'][0] for x in data['major_species_centroid'][1:]]
    plt.scatter(timepoints, deuts,40,alpha=0.5)
    plt.ylim(0,40)
    plt.xlabel('Seconds')
    plt.ylabel('Exchanged amides')
    plt.text(0.02,0.98,'Maximum labeling: %.1f%%' % (100*p.backexchange),transform=ax.transAxes,
             horizontalalignment='left',verticalalignment='top',fontsize=10)
    
    ax = plt.subplot(gs[4:9, 0])
    
    
    plt.scatter(np.arange(len(fit))[p.allowed] - len(p.Nterm) + 1,
                fit[p.allowed],
                facecolor=np.where(p.hbonds[p.allowed] == 1, 'blue', 'none'), linewidth=1.0)
    

    for k in range(-4,int(max(fit[p.allowed == True]))+1):
        plt.plot([0,50],[k,k],linewidth=0.5,color='grey',alpha=0.5,zorder=-2)
    
    plt.xlim(0,len(fit) - len(p.Nterm)+1)
    ylim=ax.get_ylim()
    for k in [5,10,15,20,25,30,35,40,45]:
        plt.plot([k,k],ylim,linewidth=0.5,color='grey',alpha=0.5,zorder=-2)
    plt.yticks(range(-2,int(max(fit[p.allowed == True]))+2))
    plt.ylim(-2.1,int(max(fit[p.allowed == True]))+1)
    
    #plt.ylim(0,40)
    plt.xlabel('Position')
    plt.ylabel('deltaG deprotect (kcal/mol)')
    
    
    
    gau_kde=gaussian_kde(fit[p.allowed == True],bw_method=0.1)
    gau_range=np.linspace(-4,max(fit[p.allowed == True])+1,5000)
    extrema=gau_range[argrelextrema(gau_kde.evaluate(gau_range), np.greater)]
    extrema_values=[(gau_kde.evaluate(x), x) for x in extrema]
    extrema_values.sort()
    extrema_values.reverse()
    extrema=sorted([x[1] for x in extrema_values[0: min(3, len(extrema_values))]])
    extrema.reverse()
    with open('%s.data' % fn,'w') as exfile:
        exfile.write(' '.join(['%.2f' % x for x in extrema]))
    
    
    plt.title('Energy levels:  ' + '  '.join(['%.2f' % x for x in extrema]) + ' (kcal/mol)')
    
    
    
    
    ax = plt.subplot(gs[10:, 0])
    raveldata=np.ravel([data['major_species_integrated_intensities'][i] / max(data['major_species_integrated_intensities'][i]) for i in range(2,14)])
    
    
    
    ravelmodel=np.ravel(allmodel)
    
    sns.regplot(raveldata[raveldata > 0],
                ravelmodel[raveldata >0])
    
    plt.xlabel('Data (Relative peak height)')
    plt.ylabel('Model (Relative peak height)')
    
    plt.text(0.02,0.98,'R^2 %.2f\nRMSE %.2f' % (np.corrcoef(raveldata[raveldata>0], ravelmodel[raveldata > 0])[0][1]**2.0,
                                                rmse(raveldata[raveldata > 0], ravelmodel[raveldata > 0]) ) ,transform=ax.transAxes,
             horizontalalignment='left',verticalalignment='top',fontsize=10)
    
    
    
    plt.subplots_adjust(left=0.1,right=0.95,top=0.95,bottom=0.07)
    plt.suptitle(title)
    plt.savefig(fn)
    plt.close()
    

def prep_sa_input(p, rect,rates=[],active_rates=[]):
    too_slow = len(p.rates[p.rates == 0])
    too_fast = len([x for x in p.Nterm if x != 'P']) + len([x for x in p.Cterm if x != 'P'])
    
    if len(p.Nterm) > 0:
        if p.Nterm[0] != 'P': too_fast -= 1
    
    if active_rates==[]:
        active_rates=np.exp(sorted(rates)[too_slow: (len(rates) if too_fast == 0 else -too_fast)])
    else:
        active_rates=np.exp(sorted(active_rates))
    
    
    n=len(active_rates)
    
    fe_grid = p.calc_fe_from_rate_grid(active_rates)

    pair_energies=np.zeros((n,n,n,n))#r1, r2, r1_energy, r2_energy

    new_pairlist=[]
    for aa1, aa2, d in p.pairlist:
        if p.hbonds[aa1] == 1 and p.hbonds[aa2] == 1 and p.allowed[aa1] and p.allowed[aa2]:
            new_pairlist.append((aa1, aa2, d))
    new_pairlist=np.array(new_pairlist)


    active_aa1s = np.array([p.active_index[x] for x in new_pairlist[:,0].astype(int)])
    active_aa2s = np.array([p.active_index[x] for x in new_pairlist[:,1].astype(int)])


    for aa1, aa2, d in zip(active_aa1s, active_aa2s, new_pairlist[:,2]):
        for rate_idx1 in range(n):
            for rate_idx2 in range(n):
                pair_energies[aa1, aa2, rate_idx1, rate_idx2] = -np.log(rect(d,
                                                                     abs(fe_grid[aa1, rate_idx1] -
                                                                         fe_grid[aa2, rate_idx2]  )))
                
    return {'protein': p,
            'active_rates': active_rates,
            'fe_grid': fe_grid,
            'pair_energies': pair_energies,
            'active_aa1s': active_aa1s,
            'active_aa2s': active_aa2s,
            'n': n}

class hx_mapping(Annealer):
    """Test annealer with a travelling salesman problem."""

    def move(self):
        """Swaps two cities in the route."""
        a = random.randint(0, len(self.state) - 1)
        b = random.randint(0, len(self.state) - 1)
        self.state[a], self.state[b] = self.state[b], self.state[a]

    def energy(self):
        """Calculates the length of the route."""

        return (self.weights['pair_e'] * self.pair_e() + 
                (self.weights['full_burial_e'] * self.full_burial_e()) +    #75
                (self.weights['hbond_burial_e'] * self.hbond_burial_e()) +   #50
                (self.weights['hbond_rank_e'] * self.hbond_rank_e()) +    #100
                (self.weights['top_stdev'] * self.top_stdev())  +
                (self.weights['distance_to_ss_e'] * self.distance_to_ss_e()) +
                (self.weights['distance_to_nonpolar_e'] * self.distance_to_nonpolar_e()))       #80
        #return 1000 * self.hbond_rank_e()
    
    def load_input(self, prepped):
        self.protein = prepped['protein']
        self.active_rates = prepped['active_rates']
        self.fe_grid = prepped['fe_grid']
        self.pair_energies = prepped['pair_energies']
        self.active_aa1s = prepped['active_aa1s']
        self.active_aa2s = prepped['active_aa2s']
        self.n = prepped['n']
        if 'weights' in prepped:
            self.weights = prepped['weights']
        else:
            self.weights={'pair_e': 50,
               'full_burial_e': 120,
               'hbond_burial_e': 14,
               'hbond_rank_e': 60,
               'distance_to_ss_e': 60,
               'distance_to_nonpolar_e': 45,
               'top_stdev': 10}
        
        self.hbond_ranks_best = sum(np.arange(self.n - sum(self.protein.hbonds), self.n))
        self.hbond_ranks_worst = sum(np.arange(sum(self.protein.hbonds)))
        
    
    def pair_e(self):
        e = []
        for aa1, aa2 in zip(self.active_aa1s, self.active_aa2s):
            e.append(self.pair_energies[aa1, aa2, self.state[aa1], self.state[aa2]])
        return np.average(e)
    
    def full_burial_e(self):
        curr_fes=np.array([self.fe_grid[aa][i] for aa,i in enumerate(self.state)])
        if min(curr_fes) == max(curr_fes): return 0.0
        return -np.corrcoef(np.maximum(curr_fes, 0), np.maximum(60, self.protein.burial[self.protein.idx_to_full]))[0][1]
    
    def hbond_burial_e(self):
        curr_fes=np.array([self.fe_grid[aa][i] for aa,i in enumerate(self.state)])
        hbond_fes=curr_fes[self.protein.hbonds[self.protein.idx_to_full] == 1]
        if min(hbond_fes) == max(hbond_fes): return 0.0
        return -np.corrcoef(np.maximum(hbond_fes, 0),
                np.maximum(60, self.protein.burial[(self.protein.hbonds == 1) & (self.protein.allowed) == True]))[0][1]
    
    def distance_to_nonpolar_e(self):
        curr_fes=np.array([self.fe_grid[aa][i] for aa,i in enumerate(self.state)])
        if min(curr_fes) == max(curr_fes): return 0.0
        return -np.corrcoef(np.maximum(curr_fes, 0),
                self.protein.distance_to_nonpolar[self.protein.idx_to_full])[0][1]
    
    def distance_to_ss_e(self):
        curr_fes=np.array([self.fe_grid[aa][i] for aa,i in enumerate(self.state)])
        if min(curr_fes) == max(curr_fes): return 0.0
        return -np.corrcoef(np.maximum(curr_fes, 0),
                self.protein.distance_to_ss[self.protein.idx_to_full])[0][1]
    
    def hbond_rank_e(self):
        curr_fes=np.array([self.fe_grid[aa][i] for aa,i in enumerate(self.state)])
        argsort = np.argsort(curr_fes)
        
        ranks = np.empty(len(curr_fes), int)
        ranks[argsort] = np.arange(len(curr_fes))
        hbond_ranks = ranks[self.protein.hbonds[self.protein.idx_to_full] == 1]
        
        return -  (float(sum(hbond_ranks)) - self.hbond_ranks_worst) / (self.hbond_ranks_best - self.hbond_ranks_worst)
    
    def top_stdev(self):
        curr_fes=np.array([self.fe_grid[aa][i] for aa,i in enumerate(self.state)])
        return np.std(sorted(curr_fes)[-3:]) #/ np.median(curr_fes[curr_fes > 0])
    
    def full_fes(self, sa_out=None):
        if sa_out==None: sa_out=self.state
        
        curr_fes=np.array([self.fe_grid[aa][i] for aa,i in enumerate(sa_out)])
        
        out=np.zeros((len(self.protein.sequence)))
        out[self.protein.idx_to_full] = curr_fes
        out[self.protein.rates == 0] = np.nan
        if len(self.protein.Cterm) > 0:
            out[-len(self.protein.Cterm):] = -2
        if len(self.protein.Nterm) > 2:
            out[1:len(self.protein.Nterm)] = -2

        return out
#xs = np.logspace(-1,5,5000)
#timepoints=[0.25, 5, 13, 20, 60, 150, 450, 1260, 3000, 6900, 13860, 23100, 30600, xs[-1], xs[-1]]
#p = hxprot('%s' % datafilename.split('_z')[0].replace('hx_','').replace('_NOSUMO',''), Nterm=Nterm)

