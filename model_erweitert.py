import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt

iptg = 0.
mrna = 0.
srna = 0.
actmrna = 0.
gfp = 0.
cells = 0.01
substrate = 2.5

x0 = [iptg, mrna, srna, actmrna, gfp, cells, substrate]
t = np.linspace(0., 700., 1000)

par = {'IPTG_decay':0.015,
       'mRNA_decay_1':0.1,
       'sRNA_kM':100,
       'sRNA_decay':0.03,
       'mRNA_decay_2':0.009,
       'activation':0.02,
       'translation':1.,
       'GFP_decay':0.01,
       'n':2,
       'leaky':.4,
       'kinf':0.,
       'y':0.5,    #0.5/0.68
       'Cells_decay':0.002, #0.002/0.003
       'muemax':0.06,    #0.06/0.08
       'kmue':6    #6/8
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
    return mrna*par['activation']*srna

def eff_sRNA(srna):
    return srna*par['sRNA_decay']

def eff_actmRNA(actmrna):
    return actmrna*par['mRNA_decay_2']

def transl_actmRNA(actmrna):
    return actmrna*par['translation']

def eff_GFP(gfp):
    return gfp*par['GFP_decay']

def mue(substrate):
    return par['muemax']*(substrate/(substrate+par['kmue']))

def eff_Cells(cells):
    return cells*par['Cells_decay']

def DGL(t, x0):
    IPTG=x0[0]
    mRNA=x0[1]
    sRNA=x0[2]
    actmRNA=x0[3]
    GFP=x0[4]
    Cells=x0[5]
    Substrate=x0[6]

    dIPTG = (K_inf()*Cells) - eff_IPTG(IPTG)
    dmRNA = (transk_IPTG(IPTG)*Cells) - eff_mRNA(mRNA) - akt_sRNAmRNA(mRNA,sRNA) + (par['leaky']*Cells)
    dsRNA = (transk_IPTG(IPTG)*Cells) - eff_sRNA(sRNA) - akt_sRNAmRNA(mRNA,sRNA) + (par['leaky']*Cells)
    dactmRNA = akt_sRNAmRNA(mRNA,sRNA) - eff_actmRNA(actmRNA)
    dGFP = (transl_actmRNA(actmRNA)*Cells) - eff_GFP(GFP)
    dCells = mue(Substrate)*Cells- eff_Cells(Cells)
    dSubstrate = (-1/par['y'])*mue(Substrate)*Cells


    return [dIPTG, dmRNA, dsRNA, dactmRNA, dGFP,dCells,dSubstrate]


integrator = scipy.integrate.ode(DGL).set_integrator('lsoda').set_initial_value(x0, 0)

array = [x0]


fig, ax1 = plt.subplots()
ax2=ax1.twinx()
cnt =1
while cnt < len(t):
    array.append(integrator.integrate(t[cnt]))

    cnt += 1

array = np.array(array)
ax1.plot(t, array[:,5],label="Cells")
ax1.legend()


ind_array = [0.001,1,2,4,8,16]

for i in ind_array:
    par['kinf']=i
    integrator = scipy.integrate.ode(DGL).set_integrator('lsoda').set_initial_value(x0, 0)

    array1 = [x0]

    cnt = 1
    while cnt < len(t):
        array1.append(integrator.integrate(t[cnt]))

        cnt += 1
    array1 = np.array(array1)
    ax2.plot(t,array1[:,4],label=str(i) + "nMIPTG")
    ax2.legend()
A=(1600,3000,3400,4300,4700,5300)
B=(605,605,605,605,605,605)
ax2.scatter(B,A)
plt.show()

