import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import datetime
from mpl_toolkits.mplot3d import Axes3D


class particle():

    def __init__(self, x, y, z, mass, vx, vy, vz):

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
        self.q = quadtreeNode(self.xs, self.ys, self.zs, self.xl, self.yl, self.zl)
        for n in range(len(particles)):
            self.q.addParticle(particles[n])

    def timeStep (self, dt, assumptionValue=0.5):
        for p in self.particles:
            self.q.gravForce(p, dt=dt, assumptionValue=assumptionValue)
        for p in self.particles:
            p.move(dt)
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

particles = []

for n in range(50):
    particles.append(particle(x=np.random.randint(-10000,10000),
                                y = np.random.randint(-10000,10000),
                                z = np.random.randint(-10000,10000),
                                mass = np.random.randint(10, 100),
                                vx = np.random.randint(0,10),
                                vy = np.random.randint(0,10),
                                vz = np.random.randint(0,10)
                                ))
q = quadtree(particles=particles)
q.plot()
q.timeStep(dt=100)
q.plot()

def init():
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    phi = np.linspace(0,2*np.pi,15)
    theta = np.linspace(0,np.pi, 15)

    for p in self.particles:
        xm = p.radius * np.outer(np.cos(phi), np.sin(theta)) + p.x
        ym = p.radius * np.outer(np.sin(phi), np.sin(theta)) + p.y
        zm = p.radius * np.outer(np.ones(np.size(phi)), np.cos(theta)) + p.z
        ax.plot_surface(xm, ym, zm, rstride=1, cstride=1)



def animate(i):
    ax.clear()
    Cloud = simulation[i]
    X, Y, Z, Mass, Radius = Cloud.status()
    CMX, CMY, CMZ = Cloud.cm()
    ax.set_xlim3d(CMX-space,CMX+space)
    ax.set_ylim3d(CMY-space,CMY+space)
    ax.set_zlim3d(CMZ-space,CMZ+space)
    phi = np.linspace(0,2*np.pi,10)
    theta = np.linspace(0,np.pi, 10)

    #newCloud = newCloud.copy()
    #newCloud.timestep(dt)
    #X, Y, Z, Mass, Radius = newCloud.status()

    for j in range(len(X)):
        xm = Radius[j] * np.outer(np.cos(phi), np.sin(theta)) + X[j]
        ym = Radius[j] * np.outer(np.sin(phi), np.sin(theta)) + Y[j]
        zm = Radius[j] * np.outer(np.ones(np.size(phi)), np.cos(theta)) + Z[j]
        ax.plot_surface(xm, ym, zm, rstride=1, cstride=1, color=Cloud.clumps[j].color)

    #ax.draw()
    #ax.show()
    ax.view_init(elev=20., azim=i)


# Set up formatting for the movie files
startTime = time.time()
Writer = animation.writers['ffmpeg']
writer = Writer(fps=20, metadata=dict(artist='Me'), bitrate=18000)
#mywriter = animation.FFMpegWriter()

duration = 5000
dt = 1000
v=0.001
simulation=[]
newCloud = Cloud(N=n, V=v, Space = space, Mass = mass)
simulation.append(newCloud)
for i in range(duration):
    tempCloud=simulation[i].copy()
    tempCloud.timeStep(dt)

    simulation.append(tempCloud)


fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')


ani=animation.FuncAnimation(fig, animate, duration, interval=1, init_func = init, blit=False)
Timestamp = ('{:%Y-%b-%d_%H:%M:%S}'.format(datetime.datetime.now()))
ani.save("GravEx%iParticles_%iSpace_%s.mp4"%(n,space, Timestamp), writer=writer)
finishTime = time.time()
#ani.save("GravitationExperiment.mp4")
#plt.show()

print("Total Calculation Time is %0.4f seconds."%(finishTime-startTime))
