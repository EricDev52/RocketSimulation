from rocket import *
from visualization import plotSimulation

# Real stats of the SpaceX Falcon 9 first stage rocket
# Burn time and thrust of the first stage. The fuel mass is only the mass which is used by the first stage during launch and is roughly estimated
# But the rocket mass is the mass of everything, as that is what the first stage must lift
# It launches from the earth, but you could also set other planets for example. These variables are found in the astronomicalObject class
# The reason the real Falcon 9 is a bit faster at stage seperation than this simulation says is probably because in the real world the rocket does a gravity turn of course which also helps to get faster

F9 = Rocket(earth, 549054, 180000, 162, 7607000, DragModel(10.752, 0.4, 10.752, 0.4))
print("Thrust weight ratio at start: " + str(round(F9.thrustWeightRatio(0), 3)))

data = F9.simulation(1000, 1000)

#data.to_csv("output.csv", index=False)
plotSimulation(data)