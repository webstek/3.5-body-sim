# -*- coding: utf-8 -*-
"""
@author: Kyle Webster
"""

import numpy as np
import matplotlib.pyplot as plt

#Launch Window Testing Simulation
#Patched Conics vs. 3.5-Body System (3.5-Body System Only for Now Though)

#4-Body (Sun-Earth-Mars-Spaceship)
#Treated as Point masses

###############################################################################
###############################################################################

##########################  Initial Conditions  ###############################

#Use for referencing following vars and such in functions
index={1:'S', 2:'E', 3:'M'}


#Some Constants
#arcseconds per hour to radians conversion factor
arc_rad = np.pi/648000

#Gravitational Parameters (GM)
mu_S = 1.3271244004193938*10**11 #km^3/s^2
mu_E = 3.98600435436*10**5       #km^3/s^2
mu_M = 4.2828375214*10**4        #km^3/s^2


#Import Data for Initial Conditions
#Data starts on line 35 with titles on line 32, don't read last 84 lines
fname_S = "horizons_results (4).csv"
fname_E = "horizons_results_Earth.csv"
fname_M = "horizons_results_Mars.csv"

#Required format of data to be read from csv
dataformat = np.dtype([('RA','f8'),('DEC','f8'),('dRA*cosD','f8'),('dDEC','f8')
            ,('r','f8'), ('rdot','f8')])

#load HORIZONS results into python data
Data_S = np.loadtxt(fname_S, comments='*', delimiter=',', dtype=dataformat,
                    usecols=[0,1,2,3,4,5], skiprows=1)
Data_E = np.loadtxt(fname_E, comments='*', delimiter=',', dtype=dataformat,
                    usecols=[0,1,2,3,4,5], skiprows=1)
Data_M = np.loadtxt(fname_M, comments='*', delimiter=',', dtype=dataformat,
                    usecols=[0,1,2,3,4,5], skiprows=1)


#Initial Position Vectors above Velocity Vectors in Spherical Coordinates
R_S = [[Data_S[0][4],Data_S[0][0],Data_S[0][1]],
       [Data_S[0][5],Data_S[0][2],Data_S[0][3]]]
R_E = [[Data_E[0][4],Data_E[0][0],Data_E[0][1]],
       [Data_E[0][5],Data_E[0][2],Data_E[0][3]]]
R_M = [[Data_M[0][4],Data_M[0][0],Data_M[0][1]],
       [Data_M[0][5],Data_M[0][2],Data_M[0][3]]]


#Which have xyz position place holders of
#Note the X direction is the direction of the vernal equinox of Earth
X_S = [0,0,0]
X_E = [0,0,0]
X_M = [0,0,0]


#Define Transformation from [R,theta,phi] to [x,y,z]
def XofR_Theta_Phi(R,Theta,Phi,i):
    x = R*np.cos(Phi)*np.cos(Theta)
    y = R*np.cos(Phi)*np.sin(Theta)
    z = R*np.sin(Phi)

    globals()['X_{}'.format(index[i])] = [x, y, z]
    

#Initial Positions of Sun, Earth, and Mars in xyz
XofR_Theta_Phi(R_S[0][0],R_S[0][1],R_S[0][2], 1) #Sun
XofR_Theta_Phi(R_E[0][0],R_E[0][1],R_E[0][2], 2) #Earth
XofR_Theta_Phi(R_M[0][0],R_M[0][1],R_M[0][2], 3) #Mars


#And velocity place holders of
Xdot_S = [0,0,0]
Xdot_E = [0,0,0]
Xdot_M = [0,0,0]

#Define Transformation from [RDot, ThetaDot*cosPhi, PhiDot] to [xdot,ydot,zdot]
def XdotofSphereDot(R,Rdot,T,TdotCosP,P,Pdot, i):
    xdot = ((Rdot*np.cos(np.radians(P))-R*np.sin(np.radians(P))*Pdot*arc_rad)*
            np.cos(np.radians(T))-(R)*TdotCosP*arc_rad*np.sin(np.radians(T)))
    ydot = ((Rdot*np.cos(np.radians(P))-R*np.sin(np.radians(P))*Pdot*arc_rad)*
            np.sin(np.radians(T))+(R)*TdotCosP*arc_rad*np.cos(np.radians(T)))
    zdot = Rdot*np.sin(np.radians(P))+R*np.cos(np.radians(P))*Pdot*arc_rad
    
    globals()['Xdot_{}'.format(index[i])] = [xdot, ydot, zdot]

#Initial Velocities of Sun, Earth, and Mars in xyz
XdotofSphereDot(R_S[0][0],R_S[1][0],R_S[0][1],R_S[1][1],R_S[0][2],R_S[1][2], 1)
XdotofSphereDot(R_E[0][0],R_E[1][0],R_E[0][1],R_E[1][1],R_E[0][2],R_E[1][2], 2)
XdotofSphereDot(R_M[0][0],R_M[1][0],R_M[0][1],R_M[1][1],R_M[0][2],R_M[1][2], 3)


#IMPOTANT MATRIX THAT NEEDS EXPANDING FOR A VESSEL
#Before doing accelerations, we set up the kinematics matrix
K=np.array([[X_S[0],Xdot_S[0],0,X_S[1],Xdot_S[1],0,X_S[2],Xdot_S[2],0]
           ,[X_E[0],Xdot_E[0],0,X_E[1],Xdot_E[1],0,X_E[2],Xdot_E[2],0]
           ,[X_M[0],Xdot_M[0],0,X_M[1],Xdot_M[1],0,X_M[2],Xdot_M[2],0]])


    
#The acceleration vectors will be
Xddot_S = [0,0,0]
Xddot_E = [0,0,0]
Xddot_M = [0,0,0]

#Define a function to compute the acceleration vector
#def accX1toX2(i,j): #acceleration of index 1 to index 2
#    X1 = globals()['X_{}'.format(index[i])]
#    X2 = globals()['X_{}'.format(index[j])]
#    mu = globals()['mu_{}'.format(index[j])]
#    R = np.subtract(X2, X1)
#    
#    xddot = -mu/(np.linalg.norm(R)**3)*R[0]
#    yddot = -mu/(np.linalg.norm(R)**3)*R[1]
#    zddot = -mu/(np.linalg.norm(R)**3)*R[2]
#    
#    globals()['Xddot_{0}_{1}'.format(index[i],index[j])]=(
#                                                [xddot, yddot, zddot])

def accX1toX2(i,j): #acceleration of index 1 to index 2
    X1 = [K[i-1][0],K[i-1][3],K[i-1][6]]
    X2 = [K[j-1][0],K[j-1][3],K[j-1][6]]
    mu = globals()['mu_{}'.format(index[j])]
    R = np.subtract(X2, X1)
    
    xddot = -mu/(np.linalg.norm(R)**3)*R[0]
    yddot = -mu/(np.linalg.norm(R)**3)*R[1]
    zddot = -mu/(np.linalg.norm(R)**3)*R[2]
    
    globals()['Xddot_{0}_{1}'.format(index[i],index[j])]=(
                                                [xddot, yddot, zddot])
    
    
#Define another function that adds the accelerations
def totalaccel(i): #total acceleration of object i
    for k in range(1,len(index)+1):
        if k != i:
            for u in range(0,3):
                globals()['Xddot_{}'.format(index[i])][u] = (
                      globals()['Xddot_{}'.format(index[i])][u]+
                      globals()['Xddot_{0}_{1}'.format(index[i],index[k])][u])


#Now calculate all the accellerations between bodies
def calcaccels():
    for i in range(1,len(index)+1):
        for j in range(1,len(index)+1):
            if j != i:
                accX1toX2(i,j)

        totalaccel(i)

#Calculate the Accelerations of the Objects
calcaccels()

#Now we can put these into the kinematic matrix by a function
def newaccels():
    K[0][2]=Xddot_S[0]
    K[1][2]=Xddot_S[1]
    K[2][2]=Xddot_S[2]
    K[0][5]=Xddot_E[0]
    K[1][5]=Xddot_E[1]
    K[2][5]=Xddot_E[2]
    K[0][8]=Xddot_M[0]
    K[1][8]=Xddot_M[1]
    K[2][8]=Xddot_M[2]
    
#Pout the New accelerations into the kinematics matrix
newaccels()

#Print the Kinematics Matrix to Check Values
print(K)

x=K[:,0]
y=K[:,3]

#Coloring Points
color=['orange', 'b', 'r']

#Actual Plot
plt.scatter(x, y, c=color)
plt.grid(b='true', which='both')
plt.show()
###############################################################################
###############################################################################





########################  Iterative Function  #################################
#We want to take the initial conditions and move the xyz by t times xdot, ydot,
# zdot, then change xdot, ydot, zdot, by t times xddot, yddot, zddot, and
#finally compute new xddot, yddot, zddot, for the new positions of objects

n = 1                    #Number of Iterations to run
t = 3600 * 24              #Interval Between steps - in seconds
K_prime = np.zeros((3,9))  #Place Holder Matrix

#Define the transformation matrix that will give iterated x and xdot
T = np.identity(9)
T[1][0]=t
T[2][1]=t
T[4][3]=t
T[5][4]=t
T[7][6]=t
T[8][7]=t


#define a function that produces a new kinematics matrix after a time interval
def iterate():
    global K
    for q in range(0,n):
        K = np.matmul(K, T) #adds t t*xdot to x, and t*xddot to xdot
    
        calcaccels()  #new accelerations
    
        newaccels()   #put the new accels into K
        
        print(K)
        
        x=K[:,0]
        y=K[:,3]

#Coloring Points
        color=['orange', 'b', 'r']

#Actual Plot
        plt.scatter(x, y, c=color)
        plt.grid(b='true', which='both')
        plt.show()

#Now Call the Iterate Function which runs n times in steps of t seconds
iterate()


#########################  Plot of Positions  #################################

##Note that z position is not being taken into account
#x=K[:,0]
#y=K[:,3]
#
##Coloring Points
#color=['orange', 'b', 'r']
#
##Actual Plot
#plt.scatter(x, y, c=color)
#plt.grid(b='true', which='both')
#plt.show()