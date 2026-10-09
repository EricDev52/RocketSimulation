import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

thrustForce = 1 # Is normally rocket.thrustForce. But its only neccesarry to know if thrustSign should be + or -

def plotSimulation(dataset : pd.DataFrame): # Plots the data of a rocket simulation using matplotlib
    df = dataset.copy()
    x : list[float] = df["t"]
    massY : list[float] = df["mass"]
    thrustY : list[float] = df["thrust"]
    gravityY : list[float] = df["gravity"]
    dragY : list[float] = df["drag"]
    forceY : list[float] = df["force"]
    accelerationY : list[float] = df["acc"]
    velocityY : list[float] = df["vel"]
    altitudeY : list[float] = df["alt"]
    gY : list[float] = df["g"]
    densityY : list[float] = df["density"]

    velocityY = np.array(velocityY) * 3.6 # For the plot it should show km/h and km
    altitudeY = np.array(altitudeY) / 1000

    figure, axis = plt.subplots(3, 2, sharex="all", gridspec_kw=dict(height_ratios=[4, 4, 1], width_ratios=[1, 1], hspace=0, wspace=0, left=0.1, right=0.85, top=0.92))
    axis[0, 1].yaxis.tick_right(); axis[1, 1].yaxis.tick_right(); axis[2, 1].yaxis.tick_right()
    axis[2, 0].set_xlabel("Time (s)", color="white"); axis[2, 1].set_xlabel("Time (s)", color="white")
    axis[2, 0].tick_params(axis='x', labelcolor="white"); axis[2, 1].tick_params(axis='x', labelcolor="white");
    figure.patch.set_facecolor("#1a1f33")
    for r in range(3):
        for c in range(2): axis[r, c].set_facecolor("#232633")
    positiveDragY : list[float] = []
    negativeDragY : list[float] = []
    for i in range(len(x)):
        if dragY[i] > 0: positiveDragY.append(dragY[i]); negativeDragY.append(0)
        elif dragY[i] < 0: positiveDragY.append(0); negativeDragY.append(abs(dragY[i]))
        else: positiveDragY.append(0); negativeDragY.append(0)
    forces = np.vstack([abs(np.array(gravityY)), abs(np.array(thrustY)), negativeDragY, positiveDragY])
    thrustSign = "+"
    if thrustForce < 0: thrustSign = "-"
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