import glob
import numpy as np
import os

juhan = glob.glob('*f90')
gareth = glob.glob('../../gareth_patch/Horizon5-master/*f90')

old = []
for j in juhan:
    segs = j.split('.')
    if(len(segs)==3 and not 'light_cone' in j):
        comb = segs[0]+'.'+segs[2]
        old.append('../'+comb)
    else:
        old.append('')

gstring=''
for g in gareth:
    gstring+= g + ' , '

gar = []
nj = len(juhan)

for j in range(nj):
    if(len(old[j])>0):
        active = old[j].split('/')[-1]
        if(active in gstring):
            gar.append('../../gareth_patch/Horizon5-master/'+active)
        else:
            gar.append('')
    else:
        active = juhan[j].split('.')[0]+'.'+juhan[j].split('.')[-1]
        if(active in gstring):
            gar.append('../../gareth_patch/Horizon5-master/'+active)
        else:
            gar.append('')

lookup = np.array([juhan,old,gar])


touched = []
# New files
ind=(lookup[2]=='')
ni = ind.sum()
print ""
print "New files"
print "========="
for i in xrange(ni):
    print(lookup[0,ind][i])
    touched.append(lookup[0,ind][i])
print "========="    
print(nj,len(touched))
print "=========" 


print ""
print "Modified openmp but same basic with main branch"
print "========="
for j in range(nj):
    if(len(lookup[1,j])!=0 and len(lookup[2,j])!=0):
       a=lookup[1,j]
       b=lookup[2,j]
       cmd = 'diff {0} {1} > mydiff'.format(a,b)
       os.system(cmd)
       f=open('mydiff','r')
       dat=f.read()
       f.close()
       if(len(dat)==0):
          print(lookup[0,j])
          touched.append(lookup[0,j])
print "========="  
print(nj,len(touched))
print "========="

print ""
print "Unmodified openmp and same basic with main branch"
print "========="
for j in range(nj):
    if(len(lookup[1,j])==0 and len(lookup[2,j])!=0):
       a=lookup[0,j]
       b=lookup[2,j]
       cmd = 'diff {0} {1} > mydiff'.format(a,b)
       #print(cmd)
       os.system(cmd)
       f=open('mydiff','r')
       dat=f.read()
       f.close()
       if(len(dat)==0):
          print(lookup[0,j])
          touched.append(lookup[0,j])
print "========="  
print(nj,len(touched))
print "========="
        
print ""
print "Needs converging"
print "========="
for j in range(nj):
    if(not lookup[0,j] in touched):
        print(lookup[0,j])


print ""
print "Needs converging : - meld"
print "========="
for j in range(nj):
    if(not lookup[0,j] in touched):
        print('meld {0} {1}'.format(lookup[0,j],lookup[2,j]))



        

