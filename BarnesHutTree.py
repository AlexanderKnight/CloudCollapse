class particle():

    def __init__(x, y, z, mass, vx, vy, vz):

        self.x = x
        self.y = y
        self.z = z
        self.mass = mass
        self.vx = vx
        self.vy = vy
        self.vz = vz

class quadtreeNode():

    def __init__(smallx, smally, largex, largey):

        self.xs = smallx
        self.xl = largex
        self.xc = (largex+smallx)/2

        self.ys = smally
        self.yl = largey
        self.yc = (largey+smally)/2

        #self.zs = smallz
        #self.zl = largez
        #self.zc = (largez+smallz)/2

        self.nss = None
        self.nsl = None
        self.nls = None
        self.nll = None
        self.p = None

        self.cmxsum = 0.0
        self.cmysum = 0.0
        #self.cmzsum = 0.0

    def addParticle(p):
        if self.nss ==  None:
            if self.p == None:
                self.cmxsum = p.x*p.mass
                self.cmysum = p.y*p.mass
                #self.cmzsum = p.z*p.mass
                self.p = p
                return
            else:
                self.nss = quadtreeNode(self.xs, self.ys, self.xc, self.yc)
                self.nsl = quadtreeNode(self.xs, self.yc, self.xc, self.yl)
                self.nls = quadtreeNode(self.xc, self.ys, self.xl, self.cy)
                self.nll = quadtreeNode(self.xc, self.yc, self.xl, self.yl)
                tmp = self.p
                if (abs(tmp.x-p.x)<0.001 and abs(tmp.y-p.y)<0.001):
                    return
                self.cmxsum=0.0
                self.cmysum=0.0
                self.msum=0.0
                self.p = None
                self.addParticle(tmp)
                self.addParticle(p)
        else:
            self.cmxsum += p.x*p.mass
            self.cmysum += p.y*p.mass
            self.msum += p.m

            if p.x < self.xc and p.y < self.yc:
                self.nss.addParticle(p)
            elif p.x >= self.xc and p.y < self.yc:
                self.nls.addParticle(p)
            elif p.x < self.xc and p.y >= self.yc:
                self.nsl.addParticle(p)
            elif p.x >= self.xc and p.y >= self.yc:
                self.nll.addParticle(p)
            return

class quadtree ():

    def __init__(particles):

        xs = 10000
        xl = -10000
        ys = 10000
        yl = -10000
        for n in range(len(particles)):
            p = particles[n]
            if p.x < xs:
                xs=p.x
            if p.x > xl:
                xl=p.x
            if p.y < ys:
                ys = p.y
            if p.y > yl:
                yl = p.y
        xs -= 0.5
        ys -= 0.5
        xl += 0.5
        yl += 0.5
        q = quadtreeNode(xs, ys, xl, yl)
        for n in range(len(particles)):
            q.addParticle(particles[n])
