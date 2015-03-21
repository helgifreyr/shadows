from math import *

file = open('harmonic-x^2.dat','w')

for i in range(-10,11):
    file.write(str(i)+'\t'+'%10.6f'%(i*i)+'\n')
