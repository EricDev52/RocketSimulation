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
'''

from astronomicalObject import *
import matplotlib.pyplot as plt
import numpy as np
from math import ceil

class DragModel:
    def __init__(self, areaTop, coefficientTop, areaBotton, coefficientBottom) -> None:
        self.areaTop = areaTop # The reference area size where drag is applied when launching in m²
        self.coefficientTop = coefficientTop # The drag coefficient of the area at the top ### Infos about the coefficient: https://en.wikipedia.org/wiki/Drag_coefficient
        self.areaBotton = areaBotton # The reference area size where drag is applied when falling down in m². This is also where you would set a parachute for example
        self.coefficientBottom = coefficientBottom # The drag coefficient of the area at the bottom

class Rocket:
    def __init__(self, astronomicalObject : AstronomicalObject, mass : float, fuelMass : float, burnTime : float, thrustForce : float, dragModel : DragModel = None, velocity : float = 0, altitude : float = 0) -> None:
        # Constants of the rocket:
        self.astronomicalObject = astronomicalObject # The astronomical object from which the rocket is launched
        self.dragModel = dragModel # The drag model of the rocket which declares for example its cross-sectional area
        self.mass = mass # The mass of the complete rocket including fuel before launch
        self.fuelMass = fuelMass # The mass that gets removed during launch like fuel and liquid oxygen
        self.burnTime = burnTime # The time the rocket fires before running out of fuel
        self.thrustForce = thrustForce # The force that is applied by thrust of the rocket. Keep in mind that this value doesnt depend on the fuelConsumption and will be the same even when you burn less fuel

        if burnTime != 0:
            if (self.fuelMass >= self.mass): raise ValueError("Fuel mass cant be greater or equal than total mass.")
            self.fuelConsumption = fuelMass / burnTime # The mass that is consumed per second during burn
        else:
            self.fuelConsumption = 0; thrustForce = 0;
            if (self.fuelMass != 0): self.fuelMass = 0; print("Burn time is zero. Fuel mass wont be burned.")

        # Changing variables for the simulation:
        self.velocity = velocity # The velocity of the rocket
        self.altitude = altitude # The altitude of the rocket

    def currentMass(self, t): # The mass of the rocket at time t
        if t >= self.burnTime: return self.mass - self.fuelMass
        return self.mass - self.fuelConsumption * t
    
    def gravityForce(self, t): # The force of gravity on the rocket at time t
        return self.currentMass(t) * -self.astronomicalObject.gravity(self.altitude)

    def dragForce(self): # The drag force on the rocket which is based on its velocity, altitude, DragModel and the atmosphere
        if self.dragModel == None or self.astronomicalObject.atmosphere == None or self.velocity == 0: return 0
        if self.velocity > 0: return -0.5 * self.astronomicalObject.atmosphere.density(self.altitude) * (self.velocity**2) * self.dragModel.areaTop * self.dragModel.coefficientTop
        else: return 0.5 * self.astronomicalObject.atmosphere.density(self.altitude) * (self.velocity**2) * self.dragModel.areaBotton * self.dragModel.coefficientBottom

    def force(self, t): # The total force on the rocket at time t
        F = self.gravityForce(t) + self.dragForce()
        if t < self.burnTime: F += self.thrustForce
        return F

    def acceleration(self, t): # The acceleration of the rocket at time t
        return self.force(t) / self.currentMass(t)

    def thrustWeightRatio(self, t): # The thrust weight ratio of the rocket at time t
        if t >= self.burnTime or self.gravityForce(t) == 0: return 0
        return self.thrustForce / abs(self.gravityForce(t))
    
    def simulation(self, seconds : int, dataPoints : int = 100): # Simulates the rocket for the given amount of seconds and returns data of it
        initialVelocity = self.velocity
        initialAltitude = self.altitude
        x, massY, thrustY, gravityY, dragY, forceY, accelerationY, velocityY, altitudeY, gY, densityY = ([] for i in range(11)) # The data lists that we will return
        updatesPerDataPoint = ceil(100 / (dataPoints / seconds)) # Target update frame rate is 100 FPS
        updateFPS = (updatesPerDataPoint * dataPoints) / seconds
        for i in range(dataPoints):
            t = (i/dataPoints) * seconds
            # Fill the data lists
            x.append(t)
            massY.append(self.currentMass(t));
            if t < self.burnTime: thrustY.append(self.thrustForce)
            else: thrustY.append(0)
            gravityY.append(self.gravityForce(t))
            dragY.append(self.dragForce())
            forceY.append(self.force(t))
            accelerationY.append(self.acceleration(t))
            velocityY.append(self.velocity)
            altitudeY.append(self.altitude)
            gY.append(self.astronomicalObject.gravity(self.altitude))
            if self.astronomicalObject.atmosphere != None: densityY.append(self.astronomicalObject.atmosphere.density(self.altitude))
            else: densityY.append(0)
            for u in range(updatesPerDataPoint):
                # Adjust the velocity and altitude: Is updated several times per second even when you get fewer data points
                # This is particularly important when dealing with high forces, because otherwise it would be assumed that this force persists for one second
                self.velocity += self.acceleration(t + (u / updatesPerDataPoint) * (seconds / dataPoints)) / updateFPS
                self.altitude += self.velocity / updateFPS
                if self.altitude <= 0 and self.astronomicalObject.radius > 0: self.velocity = 0; self.altitude = 0 # Hit the ground

        self.velocity = initialVelocity # Set these values back to its initial values after the simulation
        self.altitude = initialAltitude
        return x, massY, thrustY, gravityY, dragY, forceY, accelerationY, velocityY, altitudeY, gY, densityY
    
    def plotSimulation(self, seconds : int, dataPoints : int = 100): # Simulates the rocket and plots the data of it afterwards using matplotlib
        x, massY, thrustY, gravityY, dragY, forceY, accelerationY, velocityY, altitudeY, gY, densityY = self.simulation(seconds, dataPoints)
        velocityY = np.array(velocityY) * 3.6 # For the plot it should show km/h and km
        altitudeY = np.array(altitudeY) / 1000

        figure, axis = plt.subplots(3, 2, sharex="all", gridspec_kw=dict(height_ratios=[4, 4, 1], width_ratios=[1, 1], hspace=0, wspace=0, left=0.1, right=0.85, top=0.92))
        axis[0, 1].yaxis.tick_right(); axis[1, 1].yaxis.tick_right(); axis[2, 1].yaxis.tick_right()
        axis[2, 0].set_xlabel("Time (s)", color="white"); axis[2, 1].set_xlabel("Time (s)", color="white")
        axis[2, 0].tick_params(axis='x', labelcolor="white"); axis[2, 1].tick_params(axis='x', labelcolor="white");
        figure.patch.set_facecolor("#1a1f33")
        for r in range(3):
            for c in range(2): axis[r, c].set_facecolor("#232633")
        positiveDragY, negativeDragY = [], []
        for i in range(len(x)):
            if dragY[i] > 0: positiveDragY.append(dragY[i]); negativeDragY.append(0)
            elif dragY[i] < 0: positiveDragY.append(0); negativeDragY.append(abs(dragY[i]))
            else: positiveDragY.append(0); negativeDragY.append(0)
        forces = np.vstack([abs(np.array(gravityY)), abs(np.array(thrustY)), negativeDragY, positiveDragY])
        thrustSign = "+"
        if self.thrustForce < 0: thrustSign = "-"
        axis[0, 0].set_title("Absolute Forces Stacked", color="white", y=0.9)
        axis[0, 0].stackplot(x, forces, labels=["Gravity (-N)", "Thrust ("+thrustSign+"N)", "Drag (-N)", "Drag (+N)"], colors=["tab:blue", "#d63d1e", "#048787", "#878704"])
        axis[0, 0].legend(fancybox=True, shadow=True); axis[0, 0].tick_params(axis='y', labelcolor="tab:blue");

        axis[0, 1].set_title("F=ma", color="white", y=0.9)
        plot2_1 = axis[0, 1].twinx()
        plot2_2 = axis[0, 1].twinx()
        axis[0, 1].tick_params(axis='y', labelcolor="tab:blue")
        plot2_1.tick_params(axis='y', labelcolor="#bf3444")
        plot2_2.tick_params(axis='y', labelcolor="orange")
        p1, = axis[0, 1].plot(x, forceY, label="Force (N)", color="tab:blue")
        p2, = plot2_1.plot(x, accelerationY, label="Acceleration (m/s²)", color="#bf3444")
        p3, = plot2_2.plot(x, massY, label="Mass (kg)", color="orange")
        axis[0, 1].legend(handles=[p1, p2, p3], fancybox=True, shadow=True)
        axis[0, 1].yaxis.tick_right()
        plot2_1.spines['right'].set_position(('outward', 50))
        plot2_2.spines['right'].set_position(('outward', 100))

        axis[1, 0].plot(x, velocityY, marker=".", color="orange", label="Velocity (km/h)"); axis[1, 0].tick_params(axis='y', labelcolor="orange"); axis[1, 0].legend(fancybox=True, shadow=True); axis[1, 0].fill_between(x, 0, velocityY, facecolor="orange", alpha=0.5)
        axis[1, 1].plot(x, altitudeY, marker=".", color="green", label="Altitude (km)"); axis[1, 1].tick_params(axis='y', labelcolor="green"); axis[1, 1].legend(fancybox=True, shadow=True); axis[1, 1].fill_between(x, 0, altitudeY, facecolor="green", alpha=0.5)
        axis[2, 0].plot(x, gY, color="tab:blue", label="g (m/s²)"); axis[2, 0].tick_params(axis='y', labelcolor="tab:blue"); axis[2, 0].legend(fancybox=True, shadow=True)
        axis[2, 1].plot(x, densityY, color="tab:blue", label="ρ (kg/m³)"); axis[2, 1].tick_params(axis='y', labelcolor="tab:blue"); axis[2, 1].legend(fancybox=True, shadow=True)
        plt.show()

if __name__ == "__main__":
    # Real stats of the SpaceX Falcon 9 first stage rocket
    # Burn time and thrust of the first stage. The fuel mass is only the mass which is used by the first stage during launch and is roughly estimated
    # But the rocket mass is the mass of everything, as that is what the first stage must lift
    # It launches from the earth, but you could also set other planets for example. These variables are found in the astronomicalObject class
    # The reason the real Falcon 9 is a bit faster at stage seperation than this simulation says is probably because in the real world the rocket does a gravity turn of course which also helps to get faster
    F9 = Rocket(earth, 549054, 180000, 162, 7607000, DragModel(10.752, 0.4, 10.752, 0.4))
    print("Thrust weight ratio at start: " + str(round(F9.thrustWeightRatio(0), 3)))
    F9.plotSimulation(600, 600)