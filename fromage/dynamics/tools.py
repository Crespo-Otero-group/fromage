## tools for PyRAIMD
## Jingbai Li Feb 13 2020
## Jingbai Li May 15 2020 add readinitcond

import time,datetime,json
from fromage.dynamics.periodic_table import Element
import numpy as np

def whatistime():
    ## This function return current time

    return datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S')

def howlong(start,end):
    ## This function calculate time between start and end

    walltime=end-start
    walltime='%5d days %5d hours %5d minutes %5d seconds' % (int(walltime/86400),int((walltime%86400)/3600),int(((walltime%86400)%3600)/60),int(((walltime%86400)%3600)%60))
    return walltime

def NACpairs(ci):
    ## This function generate a dictionary for non-adiabatic coupling pairs

    pairs={}
    n=0
    for i in range(ci):
        for j in range(i+1,ci):
            n+=1
            pairs[n]=[i+1,j+1]
            pairs[str([i+1,j+1])]=n
       	    pairs[str([j+1,i+1])]=n
    return pairs

def S2F(M):
    ## This function convert 1D string (e,x,y,z) list to 2D float array

    M=[[float(x) for x in row.split()[1:4]] for row in M]
    return M

def C2S(M):
    ## This function convert 2D complex array to 2D string array

    M=[[str(x) for x in row] for row in M]
    return M

def S2C(M):
    ## This function convert 2D string array back to 2D complex array

    M=[[complex(x) for x in row] for row in M]
    return M

def Readcoord(title):
    ## This function read coordinates from a xyz file
    ## This function return coordinates in a numpy array
    ## The elements are presented by the nuclear number
    ## This function also return a list of atomic mass in amu
    ## 1 g/mol = 1822.8852 amu

    file=open('%s.xyz'% (title)).read().splitlines()
    natom=int(file[0])
    xyz=[]
    mass=np.zeros((natom,1))
    for i,line in enumerate(file[2:2+natom]):
        e,x,y,z=line.split()
        xyz.append([e,x,y,z])
        m=Element(e).getMass()
        mass[i,0:1]=m*1822.8852

    return xyz,mass

def Readinitcond(trvm):
    ## This function read coordinates from sampled initial condition
    ## This function return coordinates in a numpy array
    ## The elements are presented by the nuclear number
    ## This function also return a list of atomic mass in amu
    ## 1 g/mol = 1822.8852 amu
    natom=len(trvm)
    xyz=[]
    velo=np.zeros((natom,3))
    mass=np.zeros((natom,1))
    for i,line in enumerate(trvm):
        e,x,y,z,vx,vy,vz,m,chrg=line
        xyz.append([e,x,y,z])
        m=Element(e).getMass()
        velo[i,0:3]=float(vx),float(vy),float(vz)
        mass[i,0:1]=float(m)*1822.8852
    
    return xyz,mass,velo    

def Printcoord(xyz):
    ## This function convert a numpy array of coordinates to a formatted string

    coord=''
    for line in xyz:
        e,x,y,z=line
        coord+='%-5s%16.8f%16.8f%16.8f\n' % (e,float(x),float(y),float(z))

    return coord

def Markatom(xyz,marks,prog):
    new_xyz=[]

    for n,line in enumerate(xyz):
        e,x,y,z=line
        e = marks[n].split()[0]
        new_xyz.append([e,x,y,z])

    return new_xyz

def GetInvR(R):
    ## This functoin convert coordinates to inverse R matrix
    ## This function is only used for interfacing with PyRAIMD

    invr=[]
    q=R[1:]
    for atom1 in R:
        for atom2 in q:
            d=np.sum((atom1-atom2)**2)**0.5
            invr.append(1/d)
        q=q[1:]

    invr=np.array(invr)
    return invr

def Checkpoint(traj):
    ## obsolete

    ## This function print current information
    ## This function append output to .log, .md.energies and .md.xyz

    Chk=traj.copy()        
    title=Chk['title']     ## title
    dir=Chk['dir']         ## output directory
    temp=Chk['temp']       ## temperature
    t=Chk['t']             ## time step size
    ci=Chk['ci']           ## ci dimension
    old_state=Chk['old']   ## the previous state or the current state before surface hopping
    state=Chk['state']     ## the current state or the new state after surface hopping
    iter=Chk['iter']       ## the current iteration
    T=Chk['T']             ## atom list
    R=Chk['R']             ## coordiantes
    V=Chk['V']             ## velocity
    Ekin=Chk['Ekin']       ## kinetic energy
    E=Chk['E']             ## potential energy
    G=Chk['G']             ## gradient
    N=Chk['N']             ## non-adiabatic coupling
    At=Chk['At']           ## population (complex array)
    hoped=Chk['hoped']     ## surface hopping detector
    natom=len(T)           ## number of atoms
    ERR=Chk['ERR']         ## error of e, g, and nac in nn adaptive sampling
    verbose=Chk['verbose'] ## print level

    ## prepare a comment line for xyz file
    cmmt='%s coord %d state %d'	% (title,iter+1,old_state)

    ## prepare the surface hopping detection section according to Molcas output format
    if   hoped == 0:
        hop_info=' A surface hopping is not allowed\n  **\n At state: %3d\n' % (state)
    elif hoped == 1:
       	hop_info=' A surface hopping event happened\n  **\n From state: %3d to state: %3d\n' % (old_state,state)
        cmmt+=' to %d CI' % (state)
    elif hoped == 2:
        hop_info=' A surface hopping is frustrated\n  **\n At state: %3d\n' % (state)

    ## prepare population and potential energy info
    pop=''.join(['%28.16f' % (x) for x in np.real(np.diag(At))])
    pot=''.join(['%28.16f' % (x) for x in E])

    ## prepare non-adiabatic coupling pairs
    pairs=NACpairs(ci)

    ## start to output
    log_info=' Iter: %8d  Ekin = %28.16f au T = %8.2f K dt = %10d CI: %3d\n Root chosen for geometry opt %3d\n' % (iter,Ekin,temp,t,ci,old_state)
    log_info+='\n Gnuplot: %s%28.16f%s\n  **\n  **\n  **\n%s\n' % (pop,E[old_state-1],pot,hop_info)

    if verbose >= 1:
        xyz=np.concatenate((T,R),axis=1)
        log_info+="""
  &coordinates in Angstrom
-------------------------------------------------------
%s-------------------------------------------------------
""" % (Printcoord(xyz))
        velo=np.concatenate((T,V),axis=1)
        log_info+="""
  &velocities in Bohr/au
-------------------------------------------------------
%s-------------------------------------------------------
""" % (Printcoord(velo))
        for n,g in enumerate(G):
            grad=np.concatenate((T,g),axis=1)
            log_info+="""
  &gradient %3d in Eh/Bohr
-------------------------------------------------------
%s-------------------------------------------------------
""" % (n+1,Printcoord(grad))
        for m,n in enumerate(N):
            nac=np.concatenate((T,n),axis=1)
       	    log_info+="""
  &non-adiabatic coupling %3d - %3d in 1/Bohr
-------------------------------------------------------
%s-------------------------------------------------------
""" % (pairs[m+1][0],pairs[m+1][1],Printcoord(nac))

        if ERR != [None,None,None]:
            log_info+="""
  &error iter %-10s
-------------------------------------------------------
  Energy   error:             %-10.4f
  Gradient error:             %-10.4f
  Nac      error:             %-10.4f
-------------------------------------------------------

""" % (iter,ERR[0],ERR[1],ERR[2])


    if dir == None:
        logpath=os.getcwd()
    else:
        logpaht=dir

    print(log_info)
    mdlog=open('%s/%s.log' % (logpath,title),'a')
    mdlog.write(log_info)
    mdlog.close() 

    energy_info='%8.2f%28.16f%28.16f%28.16f%s\n' % (iter*t,E[old_state-1],Ekin,E[old_state-1]+Ekin,pot)
    mdenergy=open('%s/%s.md.energies' % (logpath,title),'a')
    mdenergy.write(energy_info)
    mdenergy.close()

    xyz_info='%d\n%s\n%s' % (natom,cmmt,Printcoord(np.concatenate((T,R),axis=1)))
    mdxyz=open('%s/%s.md.xyz' % (logpath,title),'a')
    mdxyz.write(xyz_info)
    mdxyz.close()


#    Chk['A']=C2S(Chk['A'])
#    Chk['H']=C2S(Chk['H'])
#    Chk['D']=C2S(Chk['D'])
#    Chk['At']=C2S(Chk['At'])
#    Chk['Ht']=C2S(Chk['Ht'])
#    Chk['Dt']=C2S(Chk['Dt'])

#    with open('%s/%s.chk.json' % (logpath,title),'w') as chk_file:
#        json.dump(Chk,chk_file)

