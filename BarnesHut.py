class Quad():

    def __init__(particles=None, parent=None, children=None, com = None,
                mass=None,limits=None):

        self.particles = particles
        self.parent = parent
        self.children = children

        if self.particles != None:
            self.mass = 0.0
            for particle in self.particles:
                self.mass += particle.mass
            self.com = [0,0,0]
            for particle in self.particles:
                self.com[0] += (particle.loc[0]*particle.mass)/self.mass
                self.com[1] += (particle.loc[1]*particle.mass)/self.mass
                self.com[2] += (particle.loc[2]*particle.mass)/self.mass
        else:
            self.com = com
            self.mass = mass

        self.limits = limits

    def spawnTree ():
        particlesLoc = {'swf':[], 'swb'=[], 'sef'=[], 'seb'=[],
                        'nwf':[], 'nwb'=[], 'nef'=[], 'neb'=[]}
        for particle in self.particles:
            if  self.limits[0][0] <= particle.loc[0] < self.limits[0][1]:
                particlesLoc['swf'].append(particle)
            elif self.limits[0][1] <= particle.loc[0] < self.limits[0][2]:
                particlesLoc['sef'].append(particle)
            elif self.limits[]
        self.children = {'NW':Quad(parent=self, )}

class particle():

    def __init__(loc=[0,0,0], velocity=[0,0,0], mass = None):

        self.loc=loc
        self.velocity = velocity
        self.mass = mass
