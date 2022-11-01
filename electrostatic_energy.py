# ELEMENTARY MOLECULAR ENERGY CALCULATION USING MOLECULAR MECHANICS
# Reference A. Lancaster and G.Webster, Python for the Life Sciences

import numpy as np
import matplotlib.pyplot as plt


#1. CONSTRUCT A VECTOR CLASS
class Vector3D(object):
    """ A 3-D vector with x y and z coordinate"""
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        self.np_vector = np.array([self.x, self.y, self.z])
        
    def magnitude(self):
        """Return the magnitude of self"""
        return np.linalg.norm(self.np_vector)
    
    def unitVector(self):
        """Return the unit vector of self"""
        
        return self.np_vector/self.magnitude()
    
    def transform(self, other, degrees):
        """Create a transformation matrix with self and
        return the tranformation of the vector, other"""
    
        theta = np.radians(degrees)
        u = self.unitVector()
        sint = np.sin(theta)
        cost = np.cos(theta)
        a = cost + u[0]**2 * (1 - cost)
        b = u[0] * u[1] * (1 - cost) - u[2] * sint
        c = u[0] * u[2] * (1 - cost) + u[1] * sint
        d = u[1] * u[0] * (1 - cost) + u[2] * sint
        e = cost + u[1]**2 * (1 - cost)
        f = u[1] * u[2] * (1 - cost) - u[0] * sint
        g = u[2] * u[0] * (1 - cost) - u[1] * sint
        h = u[2] * u[1] * (1 - cost) + u[0] * sint
        i = cost + u[2]**2 * (1 - cost)
        
        rotmat = np.array([ [a, b, c], 
                         [d, e, f],
                         [g, h,i] ])
    
        tx = np.matmul(rotmat[0,:], other.np_vector)
        ty = np.matmul(rotmat[1,:], other.np_vector)
        tz = np.matmul(rotmat[2,:], other.np_vector)
        
        return Vector3D(tx, ty, tz)
    
#2. CONSTRUCT AN ATOM CLASS  
class Atom:
    """An atom with a symbol, charge and coordinates"""
    def __init__(self, name, x, y, z, q):
        self.name = name
        self.x = x
        self.y = y
        self.z = z
        self.q = q
        
    def distance(self, other):
        """Calculate the distance between self and another atom, other"""
        xd = self.x - other.x
        yd = self.y - other.y
        zd = self.z - other.z
        
        return np.sqrt(xd**2 + yd**2 + zd**2)
    
    def electrostatic(self, other):
        """Calculate the electrostatic energy between self and another atom"""
        r = self.distance(other) * 1.0e-10
        q1 = self.q * 1.6e-19
        q2 = other.q * 1.6e-19
        
        return 0.000239 * (9.0e9 * 6.02e23 * q1 * q2)/(4.0*r)
    


#3. CREATING VARIABLES FOR SIMULATION
serCa = Atom('C', -44.104, 2.133, -16.495, 0.07)  #Serine alpha carbon
serCb = Atom('C', -45.239, 1.307, -17.044, 0.05)  #Serine beta carbon
serOx = Atom('O', -44.722, 0.368, -18.048, -0.66) #Serine terminal oxygen
argNH1 = Atom('N', -45.692, 1.823, -20.906, -0.8) #Arginine Nitrogen

#Rotation axis -> serine C-C bond
rotation_axis = Vector3D(serCa.x - serCb.x, serCa.y - serCb.y,
                        serCa.z - serCb.z)
#Vector to be rotated -> dfference vector between the serine beta carbon and terminal oxygen
CO_vector = Vector3D(serOx.x -serCb.x, serOx.y -serCb.y,
                             serOx.z -serCb.z)



#4. COMPUTING ELECTRROSTATIC INTERACTION
angleData =[]
distData = []
eData = []

def electrostatic():
    """Calculate electrostatic interaction between serine oxygen and a
    arginine NH1"""
    for angle in range(0, 360, 10):
        CO_vector_prime = rotation_axis.transform(CO_vector, angle)
        serOx_new = Atom('Ox', CO_vector_prime.x + serCb.x, CO_vector_prime.y 
                      + serCb.y, CO_vector_prime.z + serCb.z, -0.66)
        dist = serOx_new.distance(argNH1)
        elec = serOx_new.electrostatic(argNH1)
        angleData.append(angle)
        distData.append(dist)
        eData.append(elec)
        # print(angle, dist, elec)

electrostatic()

#5. PLOTTING
plt.title('Distance/Electrostatic Energy vs. Bond Rotation')
plt.plot(angleData, eData, 'red', angleData, distData, 'blue')
plt.axis([0, 360, 0, 20])
plt.xlabel('Rotation Angle (Degree)')
plt.ylabel('Distance (A)   E(kCal/mol)')
plt.show()




















    

