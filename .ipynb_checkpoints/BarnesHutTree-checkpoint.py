import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import datetime
from mpl_toolkits.mplot3d import Axes3D
import copy
import pandas as pd
import vispy
import pyqtgraph

class particle():

    def __init__(self, name, x, y, z, mass, vx, vy, vz):

        self.name = str(name)
        self.x = x
        self.y = y
        self.z = z
        self.mass = mass
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.radius = (self.mass/((4/3)*np.pi*1392))**(1/3)

    def acc(self, mass, d, x, y, z, dt):

        G = 6.67408*10**(-11)
        self.vx += G*mass/d**(1.5)*(x-self.x)*dt
        self.vy += G*mass/d**(1.5)*(y-self.y)*dt
        self.vz += G*mass/d**(1.5)*(z-self.z)*dt

    def move(self, dt):

        self.x += self.vx * dt
        self.y += self.vy * dt
        self.z += self.vz * dt
        
    def trait(self):
        return {'x':self.x, 'y':self.y, 'z':self.z, 'vx':self.vx, 'vy':self.vy, 'vz':self.vz, 'mass':self.mass, 'radius':self.radius}

class quadtreeNode():

    def __init__(self, smallx, smally, smallz, largex, largey, largez):

        self.xs = smallx
        self.xl = largex
        self.xc = (largex+smallx)/2

        self.ys = smally
        self.yl = largey
        self.yc = (largey+smally)/2

        self.zs = smallz
        self.zl = largez
        self.zc = (largez+smallz)/2

        self.nsss = None
        self.nsls = None
        self.nlss = None
        self.nlls = None
        self.nssl = None
        self.nsll = None
        self.nlsl = None
        self.nlll = None
        self.p = None

        self.cmxsum = 0.0
        self.cmysum = 0.0
        self.cmzsum = 0.0
        self.msum=0.0

    def addParticle(self, p):
        if self.nsss ==  None:
            if self.p == None:
                self.cmxsum = p.x*p.mass
                self.cmysum = p.y*p.mass
                self.cmzsum = p.z*p.mass
                self.msum = p.mass
                self.p = p
                return
            else:
                self.nsss = quadtreeNode(self.xs, self.ys, self.zs, self.xc, self.yc, self.zc)
                self.nsls = quadtreeNode(self.xs, self.yc, self.zs, self.xc, self.yl, self.zc)
                self.nlss = quadtreeNode(self.xc, self.ys, self.zs, self.xl, self.yc, self.zc)
                self.nlls = quadtreeNode(self.xc, self.yc, self.zs, self.xl, self.yl, self.zc)

                self.nssl = quadtreeNode(self.xs, self.ys, self.zc, self.xc, self.yc, self.zl)
                self.nsll = quadtreeNode(self.xs, self.yc, self.zc, self.xc, self.yl, self.zl)
                self.nlsl = quadtreeNode(self.xc, self.ys, self.zc, self.xl, self.yc, self.zl)
                self.nlll = quadtreeNode(self.xc, self.yc, self.zc, self.xl, self.yl, self.zl)

                tmp = self.p
                if (abs(tmp.x-p.x)<0.001 and abs(tmp.y-p.y)<0.001 and abs(tmp.z-p.z)<0.001):
                    return
                self.cmxsum=0.0
                self.cmysum=0.0
                self.cmzsum=0.0
                self.msum=0.0
                self.p = None
                self.addParticle(tmp)
                self.addParticle(p)
        else:
            self.msum += p.mass
            self.cmxsum += p.x*p.mass
            self.cmysum += p.y*p.mass
            self.cmzsum += p.z*p.mass


            if p.x < self.xc and p.y < self.yc and p.z < self.zc:
                self.nsss.addParticle(p)
            elif p.x >= self.xc and p.y < self.yc and p.z < self.zc:
                self.nlss.addParticle(p)
            elif p.x < self.xc and p.y >= self.yc and p.z < self.zc:
                self.nsls.addParticle(p)
            elif p.x >= self.xc and p.y >= self.yc and p.z < self.zc:
                self.nlls.addParticle(p)

            elif p.x < self.xc and p.y < self.yc and p.z >= self.zc:
                self.nssl.addParticle(p)
            elif p.x >= self.xc and p.y < self.yc and p.z >= self.zc:
                self.nlsl.addParticle(p)
            elif p.x < self.xc and p.y >= self.yc and p.z >= self.zc:
                self.nsll.addParticle(p)
            elif p.x >= self.xc and p.y >= self.yc and p.z >= self.zc:
                self.nlll.addParticle(p)
            return

    def gravForce(self, p, dt, assumptionValue=0.5):

        if self.cmxsum == None:
            return
        if self.msum == 0.0:
            return 0.0
        x = self.cmxsum/self.msum
        y = self.cmysum/self.msum
        z = self.cmzsum/self.msum
        s = ((self.xl-self.xs)+(self.yl-self.ys)+(self.zl-self.zs))/2
        d = np.sqrt((x-p.x)**2+(y-p.y)**2+(z-p.z)**2)

        if d != 0.0:
            if (s/d) < assumptionValue:
                p.acc(mass=self.msum, d=d, x=x, y=y, z=z, dt=dt)
            else:
                if self.nsss == None:
                    p.acc(mass = self.msum, d=d, x=x, y=y, z=z, dt=dt)
                else:
                    self.nsss.gravForce(p, assumptionValue)
                    self.nsls.gravForce(p, assumptionValue)
                    self.nlss.gravForce(p, assumptionValue)
                    self.nlls.gravForce(p, assumptionValue)
                    self.nssl.gravForce(p, assumptionValue)
                    self.nsll.gravForce(p, assumptionValue)
                    self.nlsl.gravForce(p, assumptionValue)
                    self.nlll.gravForce(p, assumptionValue)

class quadtree ():

    def __init__(self, particles, xs=-10000, xl=10000, ys=-10000, yl=10000 , zs=-10000, zl=10000):

        self.xs = xs
        self.xl = xl
        self.ys = ys
        self.yl = yl
        self.zs = zs
        self.zl = zl
        self.particles = particles
        for n in range(len(particles)):
            p = particles[n]
            if p.x < self.xs:
                self.xs=p.x
            if p.x > self.xl:
                self.xl=p.x
            if p.y < self.ys:
                self.ys = p.y
            if p.y > self.yl:
                self.yl = p.y
            if p.z < self.zs:
                self.zs = p.z
            if p.z > self.zl:
                self.zl = p.z
        self.xs -= 0.5
        self.ys -= 0.5
        self.zs -= 0.5
        self.xl += 0.5
        self.yl += 0.5
        self.zl += 0.5
        
                
        self.timestep = 0
        pnames = []
        for p in self.particles:
            pnames.append(p.name)
        self.log = pd.DataFrame(columns = pnames)
        
        row_list = {}
        for p in self.particles:
            
            row_list.update({p.name: p.trait()})
        
        
        self.log.loc[self.timestep] = row_list

        
        self.q = quadtreeNode(self.xs, self.ys, self.zs, self.xl, self.yl, self.zl)
        for n in range(len(particles)):
            self.q.addParticle(particles[n])
            
    

    def timeStep (self, dt, assumptionValue=0.5):
        for p in self.particles:
            self.q.gravForce(p, dt=dt, assumptionValue=assumptionValue)
        for p in self.particles:
            p.move(dt)
            
        self.timestep += dt
        row_list = {}
        for p in self.particles:    
            row_list.update({p.name: p.trait()})
        
        self.log.loc[self.timestep] = row_list
        
    def plot(self):
        fig = plt.figure()
        ax = fig.add_subplot(111,projection='3d')
        phi = np.linspace(0,2*np.pi,15)
        theta = np.linspace(0,np.pi, 15)

        for p in self.particles:
            xm = p.radius * np.outer(np.cos(phi), np.sin(theta)) + p.x
            ym = p.radius * np.outer(np.sin(phi), np.sin(theta)) + p.y
            zm = p.radius * np.outer(np.ones(np.size(phi)), np.cos(theta)) + p.z
            ax.plot_surface(xm, ym, zm, rstride=1, cstride=1)
        plt.show()

    def copy(self):
        newq = copy.deepcopy(self)
        return newq
    
    def saveLog(self, name):
        self.log.to_csv(name)


def init():
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    
    phi = np.linspace(0,2*np.pi,15)
    theta = np.linspace(0,np.pi, 15)

    
    



def animate(i):
    print(i)
    ax.clear()
    q = simulation[i]

    phi = np.linspace(0,2*np.pi,10)
    theta = np.linspace(0,np.pi, 10)

    
    

    for p in q.particles:
        xm = p.radius * np.outer(np.cos(phi), np.sin(theta)) + p.x
        ym = p.radius * np.outer(np.sin(phi), np.sin(theta)) + p.y
        zm = p.radius * np.outer(np.ones(np.size(phi)), np.cos(theta)) + p.z
        ax.plot_surface(xm, ym, zm, rstride=1, cstride=1)
    
    #ax.set_xlim((10000*(q.q.cmxsum//10000))-10000,(10000*(q.q.cmxsum//10000))+10000)
    #ax.set_ylim((10000*(q.q.cmysum//10000))-10000,(10000*(q.q.cmysum//10000))+10000)
    #ax.set_zlim((10000*(q.q.cmzsum//10000))-10000,(10000*(q.q.cmzsum//10000))+10000)
    
    #ax.draw()
    #ax.show()
    ax.view_init(elev=20., azim=i/16)


# Set up formatting for the movie files
print("Starting simulation")
startTime = time.time()
#Writer = animation.writers['ffmpeg']
#writer = Writer(fps=20, metadata=dict(artist='Me'), bitrate=18000)
#mywriter = animation.FFMpegWriter()

duration = 100
dt = 100000
v=0.0001
particles = []
n=1000
print("Simulation of {0} particles, {1} second timestep, {2} m/s initial velocity, and {3} timesteps".format(n, dt, v, duration))
      
print("Creating particles")
for i in range(n):
    particles.append(particle(name = '{0}'.format(i),
                              x=np.random.randint(-5000,5000),
                              y = np.random.randint(-5000,5000),
                              z = np.random.randint(-5000,5000),
                              mass = np.random.randint(1000, 10000),
                              vx = v*np.random.random()-v*np.random.random(),
                              vy = v*np.random.random()-v*np.random.random(),
                              vz = v*np.random.random()-v*np.random.random()
                              ))
#q = quadtree(particles=particles)

#simulation=[]
Q = quadtree(particles=particles)
#simulation.append(Q)
print("Advancing simulation")
for i in range(duration):
    #tempq=simulation[i].copy()
    #tempq.timeStep(dt)
    Q.timeStep(dt)
    #simulation.append(tempq)
    print(i)

Q.saveLog('TestFile.csv')
'''
fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')

print("Animating")
ani=animation.FuncAnimation(fig, animate, duration, interval=1, init_func = init, blit=False)
print("Saving video")
Timestamp = ('{:%Y-%b-%d_%H:%M:%S}'.format(datetime.datetime.now()))
ani.save("GravEx%iParticles_%s.mp4"%(n, Timestamp), fps=30, bitrate=18000)

#ani.save("GravitationExperiment.mp4")
#plt.show()
'''
finishTime = time.time()
print("Total Calculation Time is %0.4f seconds."%(finishTime-startTime))
