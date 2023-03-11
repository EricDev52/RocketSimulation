from math import pi, sqrt, e

# The gravitational constant
G = 6.6743*(10**-11)

class Atmosphere:
    def __init__(self, seaLevelDensity, scaleHeight) -> None:
        self.seaLevelDensity = seaLevelDensity # The density of the atmosphere at altitude zero in kg/m³
        self.scaleHeight = scaleHeight # The scale height of the atmosphere in meters. Set it to zero for a constant density
        # Scale height: The height over which the density decreases by the factor of e. More infos: https://en.wikipedia.org/wiki/Scale_height
    
    def density(self, altitude): # The approximate density of the atmosphere at the given altitude
        if self.scaleHeight == 0: return self.seaLevelDensity
        return self.seaLevelDensity*e**(-altitude/self.scaleHeight)

class AstronomicalObject:
    def __init__(self, mass : float, radius : float, atmosphere : Atmosphere = None) -> None:
        self.mass = mass # The mass of the object
        self.radius = radius # The radius of the object
        self.atmosphere = atmosphere # The atmosphere of the object
        if mass == 0 or radius == 0: self.mass = 0; self.radius = 0 # If mass or radius is zero there is nothing and we are in the void
    
    def gravitationalForce(self, otherObject, distance): # The gravitational force between these two objects
        if distance == 0: return 0
        return G * ((self.mass * otherObject.mass) / (distance**2))
    
    def gravity(self, altitude): # Acceleration of gravity at the given altitude
        if self.mass == 0: return 0
        return G * (self.mass / ((self.radius + altitude)**2))

    def orbitVelocity(self, altitude): # The velocity that is needed to orbit around the object at the given altitude
        if self.mass == 0: return None
        return sqrt((G * self.mass) / (self.radius + altitude))

    def escapeVelocity(self, altitude): # The velocity that is needed to escape the gravitational field of the object at the given altitude
        if self.mass == 0: return 0
        return sqrt(( 2 * G * self.mass) / (self.radius + altitude))
    
    def orbitTime(self, altitude): # The time one orbit period will take at the given altitude
        if self.mass == 0: return None
        return 2 * pi * sqrt(((self.radius + altitude)**3) / (G * self.mass))
    
    def geostationaryAltitude(self, rotationTime): # The altitude where a satellite needs to be for a geostationary orbit
        if self.mass == 0 or rotationTime == 0: return None
        return pow((G * self.mass * rotationTime**2) / (4 * pi**2), 1/3) - self.radius

sun = AstronomicalObject(1.989*(10**30), 696342000)
earth = AstronomicalObject(5.9722*(10**24), 6371000, Atmosphere(1.2, 8500))
earthVacuum = AstronomicalObject(5.9722*(10**24), 6371000, Atmosphere(0, 0))
earthWater = AstronomicalObject(5.9722*(10**24), 6371000, Atmosphere(997, 0))
mars = AstronomicalObject(6.417*(10**23), 3389500, Atmosphere(0.02, 11100))
void = AstronomicalObject(0, 0)