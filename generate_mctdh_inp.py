from tabulate import tabulate
import numpy as np
#import matplotlib.pyplot as plt
import os, re
import math

###############
##     settings


file_LVC = 'LVC.template.txt'

myformat = "presto"

#
#allmodes=range(7,220)
#allmodes=range(102,122)
#threshold 0.4
#allmodes= range(47,72)
# threshold 0.3 
allmodes = [7 ,8 ,9 ,11 ,12 ,14 ,15 ,16 ,17 ,18 ,19 ,20 ,21 ,22 ,23 ,24 ,25 ,26 ,27 ,28 ,31 ,35 ,36 ,38 ,40 ,41 ,46 ,66 ,71 ,72 ,190 ,191]

selectedmodes=[]


# if number of SPFs should be set to different values for subsets of freqs:
modenumbers=allmodes
TDH_modes=[]

numberOfSpf=[1]
gridpoints=[10]

for i in allmodes:
    if i in selectedmodes:
        numberOfSpf.append(1)
        gridpoints.append(10)
    if i not in selectedmodes:
        TDH_modes.append(i)
        numberOfSpf.append(1)
        gridpoints.append(10)

# if only selected modes should be considered further but no additional TDH
# modenumbers=selectedmodes
# numberOfSpf=len(selectedmodes)*[1] #[9,0,4,2,10,2,4,0,0]
# gridpoints=len(selectedmodes)*[10] #[12,74,6,6,46,6,20,12,10]

with_soc=True 

ns= 10  #9
nt= 14   #14
Singlets=list(range(1,ns+1))
#Singlets.remove(2)
#Singlets.remove(3)
Triplets=range(1,nt+1)
initialState=0 # initial state 
init_M=1       # multiplicity of initial stats



#thresholds in eV
minLambda=0.001
minKappa=0.001
minDelta=0.001 # only relevant if with_soc=True

# setup for relaxation run
# runtype = "relaxation"
# tfinal="400"
# tout= "20.0"
# with_ini=False
# with_dm=False
# with_cosom=False
# title="cpmp-BSE"

# setup for propagation run, path should be set to relaxation directory in variable title
# runtype = "propagation"
# tfinal="2000"
# tout= "5.0"
# with_ini=False
# path_to_prev="./Cobalt-R"
# with_dm=True
# with_cosom=True
# title="Cobalt-3"
######################################

# # setup for propagation run, path should be set to relaxation directory in variable title
# runtype = "propagation"
# tfinal="1000"
# tout= "5.0"
# with_ini=False
# with_dm=True
# with_cosom=False
# title="cpmp-BSE-Spectrum"
# ######################################

# setup for propagation run, path should be set to relaxation directory in variable title
runtype = "propagation"
tfinal="1000"
tout= "5.0"
with_ini=False
with_dm=True
with_cosom=True
title="mctdh"
######################################


# setup for spectra run, path should be set to relaxation directory in variable title
# runtype = "propagation"
# tfinal="1500"
# tout= "0.1"
# with_ini=False
# path_to_prev="./Cobalt-R"
# with_dm=True
# with_cosom=False
# title="Cobalt"
######################################

direction ="Y" # choose XYZ or N=nothing

omega="2.34" # < carrier frequency in ev
E0 ="0.0015"  # < field amplitude>  , au  
T0="100.0"  # center Gaussian pulse in  fs 
Width="24.0" # Gaussian pulse width  (24=40fs FWHM of E^2)

jobname=title


#################

#if direction not in ["X","Y","Z"]:
#   print("you need to put x, y or z")

eV=27.211399
class parameter:
    def __init__(self, M, n, value):
        self.M=M
        self.n=n
        self.value=value*eV
        self.unit="eV"

class eps(parameter):
    def __init__(self, M, n, value):
        super().__init__(M, n, value)
        if M==1:
            self.name="eS"+str(n)
        elif M==3:
            self.name="eT"+str(n)
        else:
            print('no M')
    def setvalue(self, new):
        self.value=new

class ka(parameter):
    def __init__(self, M, n, s, value):
        self.s=s
        super().__init__(M,n,value)
        if M==1:
            self.name="kS"+str(n)+"_"+str(s)
        elif M==3:
            self.name="kT"+str(n)+"_"+str(s)
        else:
            print('no M')

class la(parameter):
    def __init__(self, M, n, m, s, value):
        self.s=s
        super().__init__(M,n,value)
        if M==1:
            self.name="lS"+str(n)+"_"+str(m)+"_"+str(s)
        elif M==3:
            self.name="lT"+str(n)+"_"+str(m)+"_"+str(s)
        else:
            print('no M')
        self.m = m

class fre():
    def __init__(self, s, value):
        self.s=s
        self.value=value*eV
        self.name="f"+str(s)
        self.unit="eV"


class spinOrbitCoupling():
    def __init__(self, V, W, n, m, re, im):
        self.n=n
        self.m=m
        self.V=V
        self.W=W
        label=["S","m","0","p"]
        self.interaction=label[W-1]+label[V-1]
        self.name_re="r" + self.interaction + str(self.n)+"_"+str(self.m)
        self.name_im="i" + self.interaction + str(self.n)+"_"+str(self.m)
        self.re=re
        self.im=im
        self.absolute=(re**2+im**2)**0.5
        self.unit="eV"

class dipoleMoment():       #these are the only parameters, that dont get converted to eV
    def __init__(self, xyz, row, col, value, M):
        self.M=M
        if M==1:
            self.n=row+1
            self.m=col+1
        else:
            self.n=row-14
            self.m=col-14
        x=["S", None, "T"]
        self.name = "d" + xyz + x[M-1] + str(self.n) + "_" + str(self.m)
        self.value=value
        self.unit="AU"

frequencies=[]
epsilons=[]
kappas=[]
lambdas=[]
plot_kap=[]
plot_lam=[]
soc_r=np.array([])
soc_i=np.array([])
dm=[]

with open(file_LVC,"r") as LVC:
    for content in LVC:         # unit transformations of objects are in the class definition
        if 'epsilon' in content:
            for i in range(int(LVC.readline())):
                h=LVC.readline().split()                
                epsilons.append(eps(int(h[0]), int(h[1]), float(h[2]))) # '-1' makes doublet to singlet
        if 'kappa' in content:
            for i in range(int(LVC.readline())):
                h=LVC.readline().split()
                if float(h[3]) != 0:
                    kappas.append(ka(int(h[0]), int(h[1]), int(h[2]), float(h[3])))  # '-1' makes doublet to singlet
                    plot_kap.append(float(h[3]))
        if 'lambda' in content:
            for i in range(int(LVC.readline())):
                h=LVC.readline().split()
                if float(h[4]) != 0:
                    lambdas.append(la(int(h[0]), int(h[1]), int(h[2]), int(h[3]), float(h[4]))) # '-1' makes doublet to singlet
                    plot_lam.append(float(h[4]))
        if 'SOC R' in content:
            h=[]
            for i in range(ns+3*nt):
                h.append(LVC.readline().split())
            soc_r=np.array(h, dtype=(float))*eV
        if 'SOC I' in content:
            h=[]
            for i in range(ns+3*nt):
                h.append(LVC.readline().split())
            soc_i=np.array(h, dtype=(float))*eV
        if 'DMX R' in content and direction=="X":
            h=[]
            for i in range(ns+3*nt):
                h.append(LVC.readline().split())
            dm=np.array(h, dtype=(float))
        if 'DMY R' in content and direction=="Y":
            h=[]
            for i in range(ns+3*nt):
                h.append(LVC.readline().split())
            dm=np.array(h, dtype=(float))
        if 'DMZ R' in content and direction=="Z":
            h=[]
            for i in range(ns+3*nt):
                h.append(LVC.readline().split())
            dm=np.array(h, dtype=(float))



## correct the values of kappa to the value of the ground state for every mode

#for mode in modenumbers:
#    for kap in kappas:
#        if kap.M == 3 and kap.n ==1 and kap.s == mode:
#            ground_mode = kap.value
#            kap.value = 0.0
#        if kap.s == mode and kap.M == 3 and kap.n != 1:
#            kap.value = kap.value - ground_mode
#        if kap.s == mode and kap.M == 1:
#            mys = kap.s
#            for kap1 in kappas:
#                if kap1.M == 3 and kap1.n == 1 and kap1.s == mys:
#                    ground_mode = kap1.value
#            kap.value = kap.value - ground_mode


#print(kappas)
#print(lambdas)
#plt.hist(plot_kap)
#plt.savefig('SetB_kappas.png')

#plt.hist(plot_lam)
#plt.savefig('SetB_lambdas.png')



# frequencies are copied from coord.file, first 6 zeros are removed
# frequencies are copied from coord.file, first 6 zeros are removed
f0=[0.0001679465  , 0.0001700425  , 0.0002294115  , 0.0002405745  , 0.0002500517  , 0.0002830851  , 
0.0003238644  , 0.0003404039  , 0.0003987705  , 0.0004723098  , 0.0005131345  , 0.0005541416  , 
0.0006069495  , 0.0006529229  , 0.0006629013  , 0.0006881434  , 0.0007439130  , 0.0007737570  , 
0.0008485720  , 0.0008959579  , 0.0009176916  , 0.0010406671  , 0.0010544273  , 0.0010754775  , 
0.0010836334  , 0.0012100717  , 0.0012504864  , 0.0013021552  , 0.0013602485  , 0.0013617521  , 
0.0014143322  , 0.0014177039  , 0.0014639963  , 0.0014797612  , 0.0015162575  , 0.0015521158  , 
0.0016274776  , 0.0016356335  , 0.0016428325  , 0.0017098106  , 0.0017510455  , 0.0019131143  , 
0.0019249152  , 0.0019656489  , 0.0019887951  , 0.0020582792  , 0.0020697156  , 0.0021151423  , 
0.0021283556  , 0.0021428448  , 0.0021725976  , 0.0021949237  , 0.0022685996  , 0.0024157237  , 
0.0024209180  , 0.0024537691  , 0.0024632008  , 0.0026833629  , 0.0026900607  , 0.0028807889  , 
0.0028905395  , 0.0029544649  , 0.0029550572  , 0.0029976134  , 0.0030014407  , 0.0030089131  , 
0.0030132872  , 0.0031358526  , 0.0031385864  , 0.0033104970  , 0.0033225257  , 0.0033633049  , 
0.0033652641  , 0.0034848224  , 0.0035009063  , 0.0035276064  , 0.0035341219  , 0.0035843783  , 
0.0035860642  , 0.0036266611  , 0.0036286204  , 0.0036700830  , 0.0036707665  , 0.0037341907  , 
0.0037411618  , 0.0039060556  , 0.0039088806  , 0.0040716785  , 0.0040731365  , 0.0041594790  , 
0.0041597980  , 0.0042556633  , 0.0042592628  , 0.0043183585  , 0.0043191786  , 0.0045196118  , 
0.0045296358  , 0.0045318684  , 0.0045366070  , 0.0045847674  , 0.0045864988  , 0.0046782179  , 
0.0046783546  , 0.0046968533  , 0.0046970356  , 0.0047281553  , 0.0047306158  , 0.0047307980  , 
0.0047408219  , 0.0047463351  , 0.0047542176  , 0.0047620089  , 0.0047757235  , 0.0049526916  , 
0.0049559266  , 0.0049700512  , 0.0049716004  , 0.0050994967  , 0.0051011370  , 0.0051429186  , 
0.0051442399  , 0.0051667027  , 0.0051673405  , 0.0051957721  , 0.0051963188  , 0.0052851674  , 
0.0052869444  , 0.0053790279  , 0.0053820807  , 0.0053950207  , 0.0053957497  , 0.0054476008  , 
0.0054490132  , 0.0054804064  , 0.0054819100  , 0.0055821949  , 0.0055854755  , 0.0058418150  , 
0.0058503353  , 0.0058997715  , 0.0059005917  , 0.0059192271  , 0.0059226899  , 0.0059564068  , 
0.0059668408  , 0.0060339101  , 0.0060378285  , 0.0060799746  , 0.0060801113  , 0.0061403916  , 
0.0061414852  , 0.0062892927  , 0.0062950792  , 0.0066145240  , 0.0066243656  , 0.0066575358  , 
0.0066576269  , 0.0067053773  , 0.0067124396  , 0.0067184540  , 0.0067265643  , 0.0067724466  , 
0.0067774130  , 0.0067804657  , 0.0067832451  , 0.0068774245  , 0.0068806595  , 0.0069108225  , 
0.0069142853  , 0.0069407576  , 0.0069479111  , 0.0073275905  , 0.0073310078  , 0.0073541995  , 
0.0073572978  , 0.0073864128  , 0.0073914248  , 0.0074508850  , 0.0074541655  , 0.0074707050  , 
0.0074733933  , 0.0075351772  , 0.0075395513  , 0.0080267603  , 0.0080300864  , 0.0137923477  , 
0.0137925300  , 0.0141740320  , 0.0141740775  , 0.0143290841  , 0.0143291297  , 0.0145760375  , 
0.0145769032  , 0.0145954019  , 0.0145958576  , 0.0145976801  , 0.0145979990  , 0.0146123515  , 
0.0146133083  , 0.0146471619  , 0.0146476631  , 0.0146503513  , 0.0146525384  , 0.0146678477  , 
0.0146686678  , 0.0146952768  , 0.0146954135  , 0.0147057108  , 0.0147060298  , 0.0147202455  , 
0.0147212479  , 0.0147460344  , 0.0147472190]

for i in range(len(f0)):
    frequencies.append(fre(i+7, f0[i]))

#setting S_0 as ground state
#of=epsilons[15].value
#for i in epsilons:
#    i.setvalue(i.value-of)

#organizing DM parameters
if with_dm:
    dm_S=[]
    dm_T=[]
    init = initialState 		#Singlets + initial n - 1 (python indexes)
    if init_M==1:
        for j in range(ns):
#            print (init, j, dm[init,j])
            if dm[init, j] != 0:
#               print(dm[init,j])# and j in Singlets:
               dm_S.append(dipoleMoment(direction, init, j, dm[init, j], 1))
    if init_M==3:
        for j in range(ns,ns+nt*3+1):
            print (init+ns, j, dm[init,j])
            if dm[init+ns, j] != 0:# and j in Triplets:
                dm_T.append(dipoleMoment(direction, init+ns, j, dm[init+ns, j], 3))

    dms=[dm_T,dm_S]

if with_soc:
    soc=[]  
    couSOC=0
    label=[[1,range(0,ns),1],[2,range(ns,ns+nt),-(ns-1)],[3,range(ns+nt,ns+2*nt),-(ns+nt-1)],[4,range(ns+2*nt,ns+3*nt),-(ns+2*nt-1)]]
    for i in range(ns+nt*3):
        for j in range(ns+nt*3):
            if i != j:
                V=0
                W=0
                n=0
                m=0
                for k in label:
                    if i in k[1]:
                        W=k[0]
                        n = i + k[2]
                    if j in k[1]:
                        V=k[0]
                        m = j + k[2]
                if ((V==1 and m in Singlets) or (V>1 and m in Triplets)) and ((W==1 and n in Singlets) or (W>1 and n in Triplets)):
                    if n!=m:
                        e1=0
                        e2=0
                        for epsilon in epsilons:
                            if V==1:
                                if epsilon.M==1 and epsilon.n==m:
                                    e1 = epsilon.value
                            else:
                                if epsilon.M==3 and epsilon.n==m:
                                    e1 = epsilon.value
                            if W==1:
                                if epsilon.M==1 and epsilon.n==n:
                                    e2 = epsilon.value
                            else:
                                if epsilon.M==3 and epsilon.n==n:
                                    e2 = epsilon.value
                        de = e1 - e2
                        normV_squ=(soc_r[i,j])**2 + (soc_i[i,j])**2
                        delta= math.sqrt(normV_squ) #4*normV_squ/(de**2+4*normV_squ)           # <----- delta(SOC)
                    else:
                        delta=0
                    if soc_r[i,j]!=0 and W<V:
                        couSOC+=1
                    if soc_i[i,j]!=0:
                        couSOC+=1
                    if delta >= minDelta:
                        soc.append(spinOrbitCoupling(V, W, n, m, soc_r[i,j], soc_i[i,j]))



sel_eps=[epsilons[0]] if (1 in Triplets or initialState==1 or initialState==0) else []

for i in epsilons:
    if i.value!=0 and (i.M ==1 and i.n in Singlets):
        sel_eps.append(i)

for i in epsilons:
    if i.value!=0 and (i.M == 3 and i.n in Triplets):
        sel_eps.append(i)


#for i in sel_eps:
#    print(i)

sel_ka=[]
counter=0
full=0
for i in kappas:
    if (i.M == 1 and i.n in Singlets) and i.s in modenumbers:
        full+=1
        if abs(i.value) >= minKappa:
           sel_ka.append(i)
           counter+=1
for i in kappas:
    if (i.M == 3 and i.n in Triplets) and i.s in modenumbers:
        full+=1
        if abs(i.value) >= minKappa:
           sel_ka.append(i)
           counter+=1
print("Kappa terms: ", counter, "/", full)

counter=0
sel_la=[]
full=0
for i in lambdas:
    if i.M ==1 and i.n in Singlets and i.m in Singlets and i.s in modenumbers:
       full+=1
       if abs(i.value) >= minLambda:
          sel_la.append(i)
          counter+=1
for i in lambdas:
    if i.M == 3 and i.n in Triplets and i.m in Triplets and i.s in modenumbers:
       full+=1
       if abs(i.value) >= minLambda:
          sel_la.append(i)
          counter+=1
print("Lambda terms: ", counter, "/", full)

newModes=[]
for i in sel_la:
    if i.s not in newModes:
        newModes.append(i.s)
for i in sel_ka:
    if i.s not in newModes:
        newModes.append(i.s)
discardedModes=list(set(modenumbers) - set(newModes))
modenumbers=sorted(newModes)

minims=str(minKappa)+"/"+str(minLambda)+"/"+str(minDelta)
f="OP_DEFINE-SECTION\ntitle\n" + title + """
end-title
end-op_define-section

# modes """+str(modenumbers)+""" are used
# minimal Kappa/Lambda/Delta="""+minims

if bool(discardedModes):
    f+="\n# "+str(discardedModes)+" are modes where kappas and lambdas fell under the threshholds and got discarded"
f+="\nPARAMETER-SECTION\n"


sel_fre=[]
for i in frequencies:
    if i.value!=0 and i.s in modenumbers:
        sel_fre.append(i)

def writeParameters(l):
    s=""
    for i in l:
        s += i.name + " = " + str(i.value) + ", " + i.unit + "\n"
    return s

f += writeParameters(sel_fre)
f += writeParameters(sel_eps)
f += writeParameters(sel_ka)
f += writeParameters(sel_la)

if with_soc:
#SOC Parameters
    counter=0
#real     this will be asumed to be symmetrical
    for i in soc:
        if i.re!=0 and i.W<i.V:
            f += i.name_re + " = " + str(i.re) + ", " + i.unit + "\n"
            counter+=1
#imaginary
        if i.im!=0:
            f += i.name_im + " = " + str(i.im) + ", " + i.unit + "\n"
            counter+=1
    print("SOC terms: ", counter, "/", couSOC)

# DM Parameters
if with_dm:
    if init_M==3:
        f += writeParameters(dm_T)
    if init_M==1:
        f += writeParameters(dm_S)

f += "end-parameter-section\n\n"

if with_cosom:
    f+="""LABELS-SECTION
cosom=cos[omega]
GA=gauss[SIG,T0]
end-labels-section\n\n"""
f+="""HAMILTONIAN-SECTION
----------------------------------------
"""

f=f.replace("e-0", "d-0")

######################### HAMILTON table

table = [["modes","el"]]
for i in modenumbers:
    table[0].append("Q"+str(i))
if with_cosom:
    table[0].append("Time")

print (table[0])

for i in sel_fre:
    h1=["1"]
    h2=["1"]
    for m in modenumbers:
        h1+=["q^2"] if i.s==m else ["1"]
        h2+=["KE"] if i.s==m else ["1"]
    if with_cosom:
        h1+=["1"]
        h2+=["1"]
    table+=[["0.5*"+i.name]+h1]
    table+=[[i.name]+h2]

# epsilon_n
COL=len(modenumbers)+2 if with_cosom else len(modenumbers)+1

e_table=[]

for i in sel_eps:
    if i.M==1:
        S=str(i.n)
        e_table+=[[i.name, "S"+S+"&"+S]]#+ ["1"]*COL]

for u in range(3):
    for i in sel_eps:
        if i.M==3:
            S=str(i.n+u*len(Triplets)+len(Singlets))
            e_table+=[[i.name, "S"+S+"&"+S]]#+ ["1"]*COL]

# kappa_n_mode

for i in sel_ka:
    if i.M==1:
        S=str(i.n)
        h=["S"+S+"&"+S]
        for m in modenumbers:
            h+=["q"] if i.s==m else ["1"]
        if with_cosom: h+=["1"]
        table+=[[i.name]+h]

for u in range(3):
    for i in sel_ka:
        if i.M==3:
            S=str(i.n+u*len(Triplets)+len(Singlets))
            h=["S"+S+"&"+S]
            for m in modenumbers:
                h+=["q"] if i.s==m else ["1"]
            if with_cosom: h+=["1"]
            table+=[[i.name]+h]

#lambda
for i in sel_la:
    if i.M==1:
        h=["S"+str(i.n)+"&"+str(i.m)]
        for m in modenumbers:
            h+=["q"] if i.s==m else ["1"]
        if with_cosom: h+=["1"]
        table+=[[i.name]+h]


for u in range(3):
    for i in sel_la:
        if i.M==3:
            h=["S"+str(i.n+u*len(Triplets))+"&"+str(i.m+u*len(Triplets))]
            for m in modenumbers:
                h+=["q"] if i.s==m else ["1"]
            if with_cosom: h+=["1"]
            table+=[[i.name]+h]

#SOC
soc_table=[]

if with_soc:
    for i in soc:
    #real being symmetric
        # V = 1 means that index n (left index) is the triplet. 
        # W = 1 means that index m (right index) is the triplet, so we check which one is not equal to 1 
        # and set a or b accordingly to shift the label by the number of singlet states plus a multiple of the number of triplet states
        # this multiple is determined by the index V or W which can be 2,3 or 4 and means m_s = -1, 0 or 1
        if init_M == 1:
            a= len(Singlets) + (i.V-2)*len(Triplets) if i.V>1 else 0 
            b= len(Singlets) + (i.W-2)*len(Triplets) if i.W>1 else 0
        elif init_M == 3:
            a=(i.V-2)*len(Triplets) if i.V>1 else 3*len(Triplets)
            b=(i.W-2)*len(Triplets) if i.W>1 else 3*len(Triplets)
        else:
            print('Variable init_M is not set to 1 or 3!')
            exit()
        if i.re!=0 and i.W<i.V:
            surf=["S"+str(i.n+b)+"&"+str(i.m+a)]
            soc_table += [[i.name_re] + surf] #+ ["1"]*(COL-1)]
#       table += [[i.name_re] + surf]
    #imaginary being antisymmetric
        if i.im!=0:
            if i.n+a==i.m+b:
                surf=["S"+str(i.n+b)+"&"+str(i.m+a)]
            else:
                surf=["Z"+str(i.n+b)+"&"+str(i.m+a)]
            soc_table += [["I*"+i.name_im] + surf] #+ ["1"]*(COL-1)]
#       table += [["I*"+i.name_im] + surf]
#    print(table)

if with_cosom:

    dm_tab = []
#    for i in modenumbers:
#        dm_tab[0].append("Q"+str(i))
#    if with_cosom:
#        dm_tab[0].append("Time")
    if init_M==1:

        for i in dm_S[1:]:
            surf=["1 "+"S"+str(i.n)+"&"+str(i.m)]
#            dm_tab += [[i.name+"*E0"] + surf + ["1"]*(len(modenumbers))]
            dm_tab += [[i.name+"*E0"] + surf]
            if with_cosom: dm_tab[-1].append(str(len(modenumbers)+2)+" cosom*GA")


f+=tabulate(table, tablefmt=myformat, headers='firstrow')
f+="\n"
f+=tabulate(e_table, tablefmt=myformat)

if with_soc:
    f+="\n"
    f+=tabulate(soc_table, tablefmt=myformat)

if with_cosom:
    f+="\n"
    f+=tabulate(dm_tab, tablefmt=myformat) #, headers='firstrow')                    
f+="""
----------------------------------------
end-hamiltonian-section

"""
if with_dm:
    f+="""
HAMILTONIAN-SECTION_DIPOLE
----------------------------------------
"""

    dm_tab = [["modes","el"]]
#    for i in modenumbers:
#        dm_tab[0].append("Q"+str(i))
#    if with_cosom:
#        dm_tab[0].append("Time")

    if init_M==3:
        for u in range(3):
           for i in dm_T:
              surf=["S"+str(len(Singlets)+i.n+u*len(Triplets))+"&"+str(len(Singlets)+i.m+u*len(Triplets))]
              dm_tab += [[i.name] + surf + ["1"]*(len(modenumbers))]

#              if with_cosom: dm_tab[-1].append(str(len(modenumbers)+1)+" cosom*GA")

    if init_M==1:
        for i in dm_S[1:]:
            surf=["1 "+"S"+str(i.n)+"&"+str(i.m)]
#            dm_tab += [[i.name] + surf]
            dm_tab += [[i.name] + surf]
#            if with_cosom: dm_tab[-1].append(str(len(modenumbers)+1)+" cosom*GA")

    f+=tabulate(dm_tab, tablefmt=myformat, headers='firstrow')
    f+="""\n----------------------------------------
end-hamiltonian-section
"""
f+="end-operator"

with open("tmp.op","w+") as operator:
    operator.write(f)
    operator.seek(0)
    #print(file.read())

#generate input-file
with open(title + ".inp","w+") as input:
	input.write("""#######################################################################
###                  Janus complex
###           LVC coupling model based on (alpha,omega)-LC-BLYP             
#######################################################################

ALLOC-SECTION
maxpar    =     60000
maxfac    =      5000
maxsub    =       183
maxkoe    =    100000
maxhop    =     50000
maxhtm    =     50000
maxsub    =       250
maxnhtmshift =     20
end-alloc-section


RUN-SECTION
"""+runtype+"""
name = """+title+"""
tfinal ="""+tfinal+""".0 tout ="""+tout+"""
psi gridpop 
end-run-section

OPERATOR-SECTION
opname = """+title+"""
alter-parameters
#field_parameters (gaussian pulse E(t)=E0*cos(omega*t)*exp(-(t-T0)^2/(2*C^2);
# FWHM=2*sqrt(2*ln(2))*C; notation SIG:=1/2*C^2)
omega = """+omega+""", ev # <carrier frequency> 
E0  = """+E0+""", au # <field amplitude> , au  # amplitude
T0  = """+T0+""" , fs # <time>  , fs 
C   = """+Width+""", fs #<width> (24=40fs FWHM of E^2)
        SIG = 0.5/C/C
    end-alter-parameters
end-operator-section

sbasis-section
no-redundancy-check
""")
	for m, spf in zip(modenumbers, numberOfSpf):
		input.write("    Q"+str(m)+" = "+str(spf)+"\n")
	input.write("""end-sbasis-section
pbasis-section
""")
	for m, gp in zip(modenumbers, gridpoints):
		input.write("    Q"+str(m)+"     HO     "+str(gp)+"   0.00       1.0    1.00\n")
	input.write("""    el     el     """+str(3*len(Triplets)+len(Singlets))+"""
end-pbasis-section

INTEGRATOR-SECTION
 VMF
 RK5 = 1.0d-7
 proj-h
#CMF/var =  0.02 ,   1.0E-5
#BS/spf  =     7 ,   1.0E-5 ,   2.5E-04
#SIL/A   =     5 ,   1.0E-5
end-integrator-section

""")
    
	if with_ini:
		input.write("INIT_WF-SECTION"+"\n")
		input.write("file="+str(path_to_prev)+"/restart \n")

	if with_ini==False:
		input.write("""INIT_WF-SECTION
build
init_state = """+str(initialState+1)+"\n")
		for m in modenumbers:
			input.write("    Q"+str(m)+"     HO    0.00       0.00      1.0\n")
		input.write("""
end-build
""")
	if with_dm and not with_cosom:
		input.write("operate = DIPOLE")
	input.write("""
end-init_wf-section

end-input
""")

table=tabulate(table, tablefmt=myformat, headers='firstrow')
print(table,file=open("Ham_only.txt", "a"))

position=184
n=100
leng=15

# print(len(modenumbers), " is modenubers in tmp.op")
# print(len(modenumbers), " is modenubers in tmp.op")

newline=[]
fold=open("tmp.op", "r")
fnew=open(title + ".op", "w")

def reformat_input_line(input_line):
    values = [v.strip() for v in input_line.split('|') if v.strip()]
    first_non_one_index = None
    for i, value in enumerate(values):
        if value != '1':
            first_non_one_index = i
            break
    second_non_one_index = None
    for i in range(first_non_one_index + 1, len(values)):
        if values[i] != '1':
            second_non_one_index = i
            break
    third_non_one_index = None
    for i in range(second_non_one_index + 1, len(values)):
        if values[i] != '1':
            third_non_one_index = i
            break
    if third_non_one_index == None:
        count_ones = second_non_one_index - first_non_one_index
        output_line = f"{values[0]} |{count_ones} {values[second_non_one_index]}"
    else:
        count_ones = third_non_one_index - second_non_one_index+1
        if count_ones == 2:
            output_line = f"{values[0]} |1 {values[1]} |2 {values[third_non_one_index]}"
        else:
            output_line = f"{values[0]} |1 {values[1]} |{count_ones} {values[third_non_one_index]}"
    return output_line

for line in fold:
    if len(line)<n and 'modes' in line and '[' not in line and "--" not in line and "discarded" not in line:
        print(line.strip(), file=fnew)     
    elif len(line)>n and 'modes' in line and '[' not in line and "--" not in line and "discarded" not in line:
        elems=line.split('|')
        for i in range(1, len(modenumbers)+3, leng):
            chunk='modes  | ' + ' | '.join(elems[i:i+leng])
            print(chunk.rstrip(), file=fnew)
    elif 'modes' not in line and '|' in line:# and 'q' in line or 'K' in line:
        if 'q' in line or 'KE' in line and 'd' not in line:
            new_line=reformat_input_line(line)
            print(new_line.strip(), file=fnew)
        elif 'dX' in line or 'dY' in line or 'dZ' in line:
            new_line=re.sub(r'\|\s+','|',line)
            print(new_line.strip(), file=fnew)
        else:
            print(line.strip(), file=fnew)
    else:
        print(line.strip(), file=fnew)
fold.close()
fnew.close()

#os.remove("tmp.op")
#os.remove("Ham_only.txt")

# #generate slurm
# with open(title+".slurm","+w") as slurm:
# 	slurm.write("""#!/bin/bash
# #SBATCH --tmp=10GB
# #SBATCH --partition=Infiniband --account=infinib
# #SBATCH --exclude=compute-1-1
# #SBATCH --output=%j.out
# #SBATCH --error=%j.err
# #SBATCH --nodes=1
# #SBATCH --cpus-per-task=2      # number of requested CPUs
# #SBATCH --job-name="""+jobname+"""          # name of the job, choose a unique name!!!
# #SBATCH --mem=10GB
# #SBATCH --time=4-00:00:00


# export SCRATCH_DIR=/scratch/$USER-$SLURM_JOB_ID    # location of the scratch directory on the computing node         (loaded from programs module)
# export CALC_DIR=/scratch/$USER-$SLURM_JOB_ID       # location of the output directory on the computing node

# module load Mctdh/85.14

# export INPUT_DIR=$(pwd)

# echo "Job was started on Node: $SLURM_NODELIST"
# echo "JobID was: $SLURM_JOB_ID"
# echo "User id used is: $USER"
# echo "Submit Directory is $SLURM_SUBMIT_DIR"

# # perform calculations
# cp """+title+""".inp $CALC_DIR
# cp """+title+""".op $CALC_DIR
# cd $CALC_DIR
# mctdh85 -mnd """+title+"""
# cp $CALC_DIR/* $INPUT_DIR
# cp -r $SCRATCH_DIR/* $INPUT_DIR
# #end-run-section""")

# #################



