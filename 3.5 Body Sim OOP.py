# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 17:52:26 2020

@author: kylej
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pandas as pd

#3.5 Body Simulation using OOP
#Bodies treated as point-masses

#%%--------------------------------------------------------------------------#
#Constants and Important Variables for the Simulation
G = 6.6743015e-20 #km^3 kg^-1 s^-2


dt = 5 #86,400s (1 day) Iterative interval in seconds
dt_hours = dt/3600 #Iteration time in hours

T = 365 #total time in days to simulate --> Earth year: 365, Mars year: 687
D = int(T * (24/dt_hours)) #Number of Iterations

dV = 11.78
dR = 6371 + 200

#Number of points and data points to display
check_rate = dt * 2000 #Number of iterations between checks
data_rate = check_rate #number of iterations bettween data outputs
plot_rate = check_rate #Number of iterations between graphs
n = int(D//check_rate)
print(n) #total number of rows in the data out table

#Animation Settings
# vid_length = 60 #Desired video length in seconds
# fps = n/vid_length #required fps to have the video vid_length seconds long

ani_name = str(f'{dt_hours}hrdt{T}d')
modifier = str('')
filename = str(f'{ani_name}{modifier}.mp4')
#%%--------------------------------------------------------------------------#
#Importing Data from HORIZONS
#function to pass row indices to the loadtxt method
def gen_rows(filePath, row_indices=[]):
    with open(filePath) as f:
        # using enumerate to track line no.
        for i, line in enumerate(f):
            #if line no. is in the row index list, then return that line
            if i in row_indices:
               yield line

#Reading txt files into python
sun_txt = gen_rows("./horizons_results_sun.txt",row_indices=[23])
earth_txt = gen_rows("./horizons_results_earth.txt",row_indices=[23])
mars_txt = gen_rows("./horizons_results_mars.txt",row_indices=[23])


#The Data in rectangular coordinates
#Format is X, p
sun_data = np.loadtxt(sun_txt, delimiter=",",usecols=[2,3,4,5,6,7])
earth_data = np.loadtxt(earth_txt, delimiter=",",usecols=[2,3,4,5,6,7])
mars_data = np.loadtxt(mars_txt, delimiter=",",usecols=[2,3,4,5,6,7])

#reshaping data
sun_data.shape = (2,3)
earth_data.shape = (2,3)
mars_data.shape = (2,3)

#setting a dictionary for accessing the data
main_csv = {'sun':sun_data, 'earth':earth_data, 'mars':mars_data, 'ship':[]}

#index dictionary
ind_dict = {'sun':0, 'earth':1, 'mars':2, 'ship':3}

#Output data array
data_set = np.array([])


#%%--------------------------------------------------------------------------#
#Physics things
#masses in kg
masses = {'sun':1.989e+30, 'earth':5.9722e+24, 'mars':6.417e+23}
# radii earth: 6371km


#Needed physics equations
#F*dt = dp
#F = -GMm~r/r^3 for radial vector ~r --> Fx = -GMm x / (x^2+y^2+z^2)^3/2

#%%--------------------------------------------------------------------------#
#Physics Functions
def _calc_Rvec(ref,body):
    R = [0,0,0]
    for i in range(3):
        R[i] = body.X[i]-ref.X[i]
    return R

def _calc_R(ref, body):
    R = _calc_Rvec(ref,body)
    return np.sqrt((R[0])**2+(R[1])**2+(R[2])**2)

def _calc_Vvec(ref, body):
    m = body.mass
    M = ref.mass
    V = [0,0,0]
    for i in range(3):
        V[i] = body.p[i]/m - ref.p[i]/M
    return V

def _calc_V(ref,body):
    V = _calc_Vvec(ref,body)
    return np.sqrt((V[0])**2+(V[1])**2+(V[2])**2)

def _V_polar(ref,body):
    V = _calc_Vvec(ref,body)
    Vtheta = np.arctan2(V[1],V[0])
    R = _calc_Rvec(ref,body)
    r = _calc_R(ref,body)
    xy_proj = np.sqrt((R[0])**2+(R[1])**2)
    Vphi = np.arccos(xy_proj/r)
    return [Vtheta, Vphi]
    
def _circ_V(ref,body):
    R = _calc_R(ref,body)
    circ_V = np.sqrt(ref.mu/R)
    print('Circular Velocity: ', circ_V)
    
def _calc_Energy(ref,body):
    r = _calc_R(ref,body)
    v = _calc_V(ref,body)
    return (v**2)/(2) - (ref.mu)/(r)
    # E = E_spfc * body.mass
    
def _calc_hvec(ref,body):
    R = _calc_Rvec(ref,body)
    V = _calc_Vvec(ref,body)
    return np.cross(R,V)
    
def _calc_eccvec(ref, body):
    V = _calc_Vvec(ref,body)
    r = _calc_R(ref,body)
    R = _calc_Rvec(ref,body)
    H = _calc_hvec(ref,body)
    VCH = np.cross(V,H)
    
    e = [0,0,0]
    for i in [0,1,2]:
        e[i] = (VCH[i])/(ref.mu) - (R[i])/(r)    
    return e

def _calc_ecc(ref, body):
    e_vec = _calc_eccvec(ref, body)
    return np.sqrt((e_vec[0])**2+(e_vec[1])**2+(e_vec[2])**2)

def _R_polar(ref,body):
    r = _calc_R(ref,body)
    R = _calc_Rvec(ref,body)
    rtheta = np.arctan2(R[1],R[0])
    xy_proj = np.sqrt((R[0])**2+(R[1])**2)
    if r!=0:
        rphi = np.arccos(xy_proj/r)
    else:
        rphi = 0
    return [r,rtheta,rphi]
      
def deltaV_prograde(ref,body,deltaV):
    V_polar = _V_polar(ref,body)
    Vtheta = V_polar[0]
    Vphi = V_polar[1]
    
    V = _calc_Vvec(ref,body)
    V[0] += deltaV * np.cos(Vtheta) * np.cos(Vphi)
    V[1] += deltaV * np.sin(Vtheta) * np.cos(Vphi)
    V[2] += deltaV * np.sin(Vphi)
    
    for i in range(3):
        body.p[i] = V[i] * body.mass
        
def deltaR(ref,body,deltaR):
    R_polar = _R_polar(ref,body)
    rtheta = R_polar[1]
    rphi = R_polar[2]
    
    X = _calc_Rvec(ref,body)
    
    dX = [0,0,0] 
    dX[0] = deltaR * np.cos(rtheta) * np.cos(rphi)
    dX[1] = deltaR * np.sin(rtheta) * np.cos(rphi)
    dX[2] = deltaR * np.sin(rphi)
    
    for i in range(3):
        body.X[i] = dX[i] + ref.X[i] + X[i]


def calc_physics():
    for body in bodies:
        F = [0,0,0]
        for source in bodies:
            if body != source:
                r = _calc_R(source,body)
                
                for i in range(3):
                    F[i] += -source.mu*body.mass*(body.X[i]-source.X[i])/(r**3)
                    
                p0 = body.p
                X0 = body.X
        for i in range(3):
            body.p[i] = p0[i] + F[i]*dt
            body.X[i] = X0[i] + (p0[i]*dt)/(body.mass) \
                              + (F[i]*dt*dt)/(2*body.mass)
                                      
def load_data():
    for body in bodies:
        body.load_data(main_csv[str(body)])
        
def export_data():
    counter = 0
    for body in bodies:
        clmS, clmE = counter*3, (counter+1)*3
        body.data = out_data[:,clmS:clmE]
        
        counter+=1
    out_csv = pd.DataFrame(data=out_data)
    return out_csv


#%%--------------------------------------------------------------------------#
#Animation/Plotting Functions
def get_artists():
    artists=list()
    for body in bodies:
        new_artist = body.get_artist()
        artists.append(new_artist)
    return artists
        
def init():
    plt.grid(b='true', which='both',zorder=1)
    plt.ylim(-psl,psl)
    plt.xlim(-psl,psl)
    plt.yticks(ticks)
    plt.xticks(ticks)
    
    artists = get_artists()
    return artists

def animate(i):
    pts = list()

    ax.set_xlim(-psl,psl)
    ax.set_ylim(-psl,psl)
    
    # x, y = earth.data[i-1,0],earth.data[i-1,1]
    # ax.set_xlim(x - 1e5, x + 1e5)
    # ax.set_ylim(y - 1e5, y + 1e5)
    for body in bodies:
        
        thisx = body.data[i-1,0]
        thisy = body.data[i-1,1]
        j = ind_dict[str(body)]
        
        pt_plot[j].set_data(thisx,thisy)
        pts.append(pt_plot[j])
    return pts
    

#%%--------------------------------------------------------------------------#
#Iterative function
def iterate():
    counter = 0
    out_data = data_set
    for i in range(D):
        calc_physics()
        if i%data_rate==0:
            counter+=1
            req_data = [sun.X[0],sun.X[1],sun.X[2],
                        earth.X[0],earth.X[1],earth.X[2],
                        mars.X[0],mars.X[1],mars.X[2],
                        ship.X[0],ship.X[1],ship.X[2]]
            
            out_data = np.append(out_data,req_data)
            
        for j in [10/9,5/4,10/7,5/3,2,5/2,10/3,5,10]:
            if i == D//j:
                completion = round(float(1/j)*100,0)
                print(f'{completion}%')
        if i==D-1:
            print('Physics complete.')
    
    out_data.shape = (counter, len(req_data))
    return out_data
            

#%%--------------------------------------------------------------------------#
#Body Class
class Body:
    def __init__(self, name, mass, c='k',s=1,zo=2):
        self.name = name    
        self.mass = mass
        self.mu = G * self.mass #km^3 s^-2
        self.c = c
        self.s = s
        self.zo = zo
        self.data = np.array([])
        
    def __str__(self):
        return self.name

    @property
    def X(self):
        return self._X
    @X.setter
    def X(self, new_X):
        self._X = new_X

    @property
    def p(self):
        return self._p
    @p.setter
    def p(self, new_p):
        self._p = new_p

    def load_data(self, data): #load the appropriate csv for the body
        X = data[0]
        p = self.mass*data[1]
        
        self._X = X
        self._p = p
        
    def get_artist(self):
        dot, = plt.plot(self.X[0],self.X[1],c=self.c,zorder=self.zo,  \
                        marker='.',ms=self.s)
        return dot
    
    def plot(self):
        ax.plot([],[],c=self.c,zorder=self.zo,\
                 marker='.',ms=self.s)

class Ship(Body):
    def load_data(self, data):
        #Setting spacecraft elements
        
        self.X = [earth.X[0],earth.X[1],earth.X[2]]
        self.p = [earth.p[0]/earth.mass,earth.p[1]/earth.mass, \
                  earth.p[2]/earth.mass]
        

#%%--------------------------------------------------------------------------#
#Creating the Instance Objects
sun = Body('sun', masses['sun'],'orange',3)
earth = Body('earth', masses['earth'],'b')
mars = Body('mars', masses['mars'],'r')
ship = Ship('ship',1,'g',2,3)

bodies = [sun, earth, mars, ship]

#Giving instance objects X, p values from HORIZONS
load_data()

deltaV_prograde(sun, ship, dV)
deltaR(sun,ship,dR)


#%%--------------------------------------------------------------------------#
#Alternatively the physics can be run but no image will be plotted.
#Plotting area for the zoomed plot
out_data = iterate()
csv_out = export_data()

csv_out.to_csv('./3.5BodySimResults.csv',header=False,index=False)
sys.exit('csv file created.')

#%%--------------------------------------------------------------------------#
#Creating animation by reading the out_data at certain points
#Figure Variables and Reading data to plot
psl = 3e8
ticks = np.linspace(-psl,psl,9)

fig = plt.figure(figsize=(4, 4))
ax = fig.add_subplot(111, autoscale_on=False,
                     xlim=(-psl, psl), ylim=(-psl, psl))
ax.set_aspect('equal')
ax.grid()

pt_plot = get_artists()

ani = FuncAnimation(fig, animate, frames=n, blit=True)
ani.save(filename,fps=60,dpi=300)
plt.show()











