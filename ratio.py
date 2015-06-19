from pylab import *
from scipy import *

n_rs = []
a_rs = []
n_ths = []
a_ths = []
n_fs = []
a_fs = []
a_proper_time = []
n_proper_time = []
r_ratio = []
th_ratio = []
f_ratio = []


for analytical,numerical in zip(open("numerical.dat"),open("analytical.dat")):
    n_time = float(numerical.split()[0].replace(r'{','').replace(r'}',''))
    a_time = float(analytical.split()[0].replace(r'{','').replace(r'}',''))
    analytical_r = float(analytical.split()[1].replace(r'{','').replace(r'}',''))
    numerical_r = float(numerical.split()[1].replace(r'{','').replace(r'}',''))
    analytical_th = float(analytical.split()[2].replace(r'{','').replace(r'}',''))
    numerical_th = float(numerical.split()[2].replace(r'{','').replace(r'}',''))
    analytical_f = float(analytical.split()[3].replace(r'{','').replace(r'}',''))
    numerical_f = float(numerical.split()[3].replace(r'{','').replace(r'}',''))
    a_proper_time.append(a_time)
    n_proper_time.append(n_time)
    n_rs.append(numerical_f)
    a_rs.append(analytical_f)
    n_ths.append(numerical_f)
    a_ths.append(analytical_f)
    n_fs.append(numerical_f)
    a_fs.append(analytical_f)
    r_ratio.append(1-numerical_r/analytical_r)
    if analytical_th != 0:
        th_ratio.append(1-numerical_th/analytical_th)
    if analytical_f != 0:
        f_ratio.append(1-numerical_f/analytical_f)

proper_time = a_proper_time

subplot(311)
plot(proper_time,r_ratio)
xlim(0,proper_time[-1])
subplot(312)
plot(proper_time,th_ratio)
xlim(0,proper_time[-1])
subplot(313)
plot(proper_time[1:],f_ratio)
xlim(0,proper_time[-1])

savefig('ratio.png')
