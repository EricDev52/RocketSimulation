'''
This script is about simulating the physics of a rocket launch
It is assumed that the fuel consumption rate and the thrust of the rocket stays constant during launch until its out of fuel
It is assumed that the rocket never rotates and only accelerates straight up and also only falls straight down
Values like force, acceleration and velocity are vectors, but since the rocket can only go up or down here, the vector is represented as a number
A negative number represents a down direction, a positive number represents an up direction, and the absolute value is the magnitude of the vector
This simulation can also be used for just simulating free fall for example. Just set fuelMass, burnTime and thrustForce to zero and set a initial velocity and altitude
Always keep in mind that when looking at the simulation plot, you only see the data for these data points. For anything between you theoretically dont know what happened there

Default units used:
- mass in kg
- time in seconds
- force in newton
- acceleration in m/s²
- velocity in m/s
- altitude in m
- volume in m³
'''

from astronomicalObject import *
import pandas as pd

class DragModel:
    def __init__(self, areaTop : float, coefficientTop : float, areaBottom : float, coefficientBottom : float) -> None:
        self.areaTop = areaTop # The reference area size where drag is applied when launching in m²
        self.coefficientTop = coefficientTop # The drag coefficient of the area at the top ### Infos about the coefficient: https://en.wikipedia.org/wiki/Drag_coefficient
        self.areaBottom = areaBottom # The reference area size where drag is applied when falling down in m². This is also where you would set a parachute for example
        self.coefficientBottom = coefficientBottom # The drag coefficient of the area at the bottom

class Rocket:
    def __init__(self, astronomicalObject : AstronomicalObject, mass : float, fuelMass : float, burnTime : float, thrustForce : float, dragModel : DragModel | None = None, velocity : float = 0, altitude : float = 0, volume : float = 0) -> None:
        # Constants of the rocket:
        self.astronomicalObject = astronomicalObject # The astronomical object from which the rocket is launched
        self.dragModel = dragModel # The drag model of the rocket which declares for example its cross-sectional area
        self.mass = mass # The mass of the complete rocket including fuel before launch
        self.volume = volume # The displaced volume in m³, assuming full immersion in the surrounding medium. The volume is only used for calculating buoyancy and is therefore optional.
        self.fuelMass = fuelMass # The mass that gets removed during launch like fuel and liquid oxygen
        self.burnTime = burnTime # The time the rocket fires before running out of fuel
        self.thrustForce = thrustForce # The force that is applied by thrust of the rocket. Keep in mind that this value doesnt depend on the fuelConsumption and will be the same even when you burn less fuel

        if burnTime != 0:
            if (self.fuelMass >= self.mass): raise ValueError("Fuel mass cant be greater or equal than total mass.")
            self.fuelConsumption = fuelMass / burnTime # The mass that is consumed per second during burn
        else:
            self.fuelConsumption = 0
            self.thrustForce = 0
            if (self.fuelMass != 0):
                self.fuelMass = 0
                print("Burn time is zero. Fuel mass wont be burned.")

        # Changing variables for the simulation:
        self.velocity : float = velocity # The velocity of the rocket
        self.altitude : float = altitude # The altitude of the rocket

    def currentMass(self, t : float) -> float: # The mass of the rocket at time t
        if t >= self.burnTime: return self.mass - self.fuelMass
        return self.mass - self.fuelConsumption * t
    
    def gravityForce(self, t : float): # The force of gravity on the rocket at time t
        return self.currentMass(t) * -self.astronomicalObject.gravity(self.altitude)

    def dragForce(self): # The drag force on the rocket which is based on its velocity, altitude, DragModel and the atmosphere
        if self.dragModel == None or self.astronomicalObject.atmosphere == None or self.velocity == 0: return 0
        if self.velocity > 0: return -0.5 * self.astronomicalObject.atmosphere.density(self.altitude) * (self.velocity**2) * self.dragModel.areaTop * self.dragModel.coefficientTop
        else: return 0.5 * self.astronomicalObject.atmosphere.density(self.altitude) * (self.velocity**2) * self.dragModel.areaBottom * self.dragModel.coefficientBottom

    def buoyancyForce(self): # The upward force from the displaced surrounding medium
        # Buoyancy is often negligible in air, but plays a major role in water, for example.
        if self.volume == 0 or self.astronomicalObject.atmosphere == None: return 0
        return self.astronomicalObject.atmosphere.density(self.altitude) * self.volume * self.astronomicalObject.gravity(self.altitude)

    def force(self, t : float): # The total force on the rocket at time t
        f = self.gravityForce(t) + self.dragForce() + self.buoyancyForce()
        if t < self.burnTime: f += self.thrustForce
        return f

    def acceleration(self, t : float): # The acceleration of the rocket at time t
        return self.force(t) / self.currentMass(t)

    def thrustWeightRatio(self, t : float): # The thrust weight ratio of the rocket at time t
        if t >= self.burnTime or self.gravityForce(t) == 0: return 0
        return self.thrustForce / abs(self.gravityForce(t))
    
    def simulation(self, seconds : float, dataPoints : int = 100, dt : float = 0.01) -> pd.DataFrame: # Simulates the rocket for the given amount of seconds and returns data of it
        '''
        Simuliert die Rakete über 'seconds'.
        - Integration erfolgt mit Schrittweite höchstens dt, bei hoher Beschleunigung kleiner.
        - Daten werden nur 'dataPoints' mal gespeichert.
        - Gibt ein pandas DataFrame zurück.
        '''
        if dt <= 0: raise ValueError("dt must be greater than zero.")
        if seconds < 0: raise ValueError("seconds must not be negative.")
        if dataPoints < 1: raise ValueError("dataPoints must be at least 1.")

        # ursprünglichen Zustand sichern (weil die Physikfunktionen darauf zugreifen)
        v0, h0 = self.velocity, self.altitude
        sampleTimes : list[float]
        if dataPoints <= 1: sampleTimes = [0.0]
        else:
            step = seconds / (dataPoints - 1)
            sampleTimes = [i * step for i in range(dataPoints)] # Die zeiten zu den ein Datenpunkt gespeichert werden soll
            sampleTimes[-1] = seconds
        
        columns = ["t","mass","thrust","gravity","drag","force","acc","vel","alt","g","density","buoyancy"]
        rows : list[list[float]] = [] # Hier werden die Daten gespeichert

        def recordRow(t : float): # Hilfsfunktion zum erstellen eines Datenpunkts (Zeile in der späteren Tabelle)
            mass = self.currentMass(t)
            thrust = self.thrustForce if t < self.burnTime else 0.0
            gravity = self.gravityForce(t)
            drag = self.dragForce()
            force = self.force(t)
            buoyancy = self.buoyancyForce()
            acceleration = self.acceleration(t)
            g = self.astronomicalObject.gravity(self.altitude)
            density = self.astronomicalObject.atmosphere.density(self.altitude) if self.astronomicalObject.atmosphere != None else 0.0     
            rows.append([t, mass, thrust, gravity, drag, force, acceleration, self.velocity, self.altitude, g, density, buoyancy])
        
        t = 0
        recordRow(t) # Save initial sample
        nextSample = 1
        
        while t < seconds and nextSample < dataPoints:
            targetTime = sampleTimes[nextSample]
            while t < targetTime:
                step = min(dt, targetTime - t) # Step = deltaTime except when time left to targetTime is smaller than deltaTime
                a = self.acceleration(t)
                # Smaller steps so that the speed does not change abruptly.
                if a != 0: step = min(step, 0.01 * max(abs(self.velocity), 1.0) / abs(a))
                if t < self.burnTime: step = min(step, self.burnTime - t)
                # Adjust the velocity and altitude: Is updated several times per second even when you get fewer data points
                # This is particularly important when dealing with high forces, because otherwise it would be assumed that this force persists for one second
                self.altitude += self.velocity * step + 0.5 * a * step**2
                self.velocity += a * step
                if self.altitude <= 0 and self.astronomicalObject.radius > 0: self.velocity = 0; self.altitude = 0 # Hit the ground
                t += step
            recordRow(targetTime) # Data record row to target time after possibly several delta time steps
            nextSample += 1
        
        # Set these values back to its initial values after the simulation
        self.velocity, self.altitude = v0, h0
        return pd.DataFrame(rows, columns=columns)