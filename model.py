import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt

iptg = 0.
mrna = 2.
srna = 0.
actmrna = 0.
gfp = 0.

x0 = [iptg, mrna, srna, actmrna, gfp]
t = np.linspace(0., 700., 1000)

par = {'IPTG_decay':0.015,
       'mRNA_decay_1':0.05,
       'sRNA_kM':.15,
       'sRNA_decay':0.03,
       'mRNA_decay_2':0.03,
       'activation':0.02,
       'translation':1.,
       'GFP_decay':0.004,
       'n':4,
       'leaky':.4,
       'kinf':1.,
       'lambda':1.1
       }

def K_inf():
    return par['kinf']

def eff_IPTG(iptg):
    return iptg*par['IPTG_decay']

def transk_IPTG(iptg):
    return iptg**par['n']/(iptg**par['n']+par['sRNA_kM']**par['n'])

def eff_mRNA(mrna):
    return mrna*par['mRNA_decay_1']

def akt_sRNAmRNA(mrna,srna):
    return mrna*par['activation']* srna

def eff_sRNA(srna):
    return srna*par['sRNA_decay']

def eff_actmRNA(actmrna):
    return actmrna*par['mRNA_decay_2']

def transl_actmRNA(actmrna):
    return actmrna*par['translation']

def eff_GFP(gfp):
    return gfp*par['GFP_decay']

def DGL(t, x0):
    IPTG=x0[0]
    mRNA=x0[1]
    sRNA=x0[2]
    actmRNA=x0[3]
    GFP=x0[4]

    dIPTG = K_inf()*np.exp(-1*par['lambda']*t) - eff_IPTG(IPTG)
    dmRNA = transk_IPTG(IPTG) - eff_mRNA(mRNA) - akt_sRNAmRNA(mRNA,sRNA) + par['leaky']
    dsRNA = transk_IPTG(IPTG) - eff_sRNA(sRNA) - akt_sRNAmRNA(mRNA,sRNA) + par['leaky']
    dactmRNA = akt_sRNAmRNA(mRNA,sRNA) - eff_actmRNA(actmRNA)
    dGFP = transl_actmRNA(actmRNA) - eff_GFP(GFP)

    return [dIPTG, dmRNA, dsRNA, dactmRNA, dGFP]


integrator = scipy.integrate.ode(DGL).set_integrator('lsoda').set_initial_value(x0, 0)

array = [x0]

cnt = 1
while cnt < len(t):
    array.append(integrator.integrate(t[cnt]))

    cnt += 1

array = np.array(array)
plt.plot(t, array[:,0])
plt.plot(t, array[:,1])
plt.plot(t, array[:,2])
plt.plot(t, array[:,3])
#plt.plot(t, array[:,4])
plt.show()

ind_array = [0,2,4,8,16]

for i in ind_array:
    par['kinf']=i
    integrator = scipy.integrate.ode(DGL).set_integrator('lsoda').set_initial_value(x0, 0)

    array = [x0]

    cnt = 1
    while cnt < len(t):
        array.append(integrator.integrate(t[cnt]))

        cnt += 1
    array = np.array(array)
    plt.plot(t,array[:,4])

plt.show()


for i in ind_array:
    par['kinf']=i
    integrator = scipy.integrate.ode(DGL).set_integrator('lsoda').set_initial_value(x0, 0)

    array = [x0]

    cnt = 1
    while cnt < len(t):
        array.append(integrator.integrate(t[cnt]))

        cnt += 1
    array = np.array(array)
    plt.plot(t,array[:,0])

plt.show()"# synbio" 
# synbio
