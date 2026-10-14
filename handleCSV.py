import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# Load simulation data
csv_path = Path("Trajectory.csv")
data = pd.read_csv(csv_path)

time = data["time"]

x = data["x"]
y = data["y"]
z = data["z"]

vx = data["vx"]
vy = data["vy"]
vz = data["vz"]


# Create one figure with a 2x2 grid
fig = plt.figure(figsize=(15, 10))

fig.suptitle("Flight Simulation", fontsize=18)


# --------------------------------
# Plot 1) Altitude vs Time
# --------------------------------

ax1 = fig.add_subplot(2, 2, 1)

ax1.plot(time, z, linewidth=2)

ax1.set_xlabel("Time (s)")
ax1.set_ylabel("Altitude (m)")
ax1.set_title("Altitude vs Time")
ax1.grid(True)


# --------------------------------
# Plot 2) Velocity vs Time
# --------------------------------

ax2 = fig.add_subplot(2, 2, 2)

ax2.plot(time, vx, label="X Velocity")
ax2.plot(time, vy, label="Y Velocity")
ax2.plot(time, vz, label="Z Velocity")

ax2.set_xlabel("Time (s)")
ax2.set_ylabel("Velocity (m/s)")
ax2.set_title("Velocity vs Time")

ax2.legend()
ax2.grid(True)


# --------------------------------
# Plot 3) 3D Trajectory
# --------------------------------

ax3 = fig.add_subplot(2, 2, 3, projection="3d")

ax3.plot(x, y, z, linewidth=2)

ax3.scatter(
    x.iloc[0], y.iloc[0], z.iloc[0],
    label="Start", s=50
)

ax3.scatter(
    x.iloc[-1], y.iloc[-1], z.iloc[-1],
    label="Landing", s=50
)

ax3.set_xlabel("X Position (m)")
ax3.set_ylabel("Y Position (m)")
ax3.set_zlabel("Altitude (m)")

ax3.set_title("3D Trajectory")
ax3.legend()


# --------------------------------
# Plot 4) Vertical Velocity vs Altitude
# --------------------------------

ax4 = fig.add_subplot(2, 2, 4)

ax4.plot(vz, z, linewidth=2)

ax4.set_xlabel("Vertical Velocity (m/s)")
ax4.set_ylabel("Altitude (m)")
ax4.set_title("Vertical Velocity vs Altitude")

ax4.grid(True)


# --------------------------------
# Display dashboard
# --------------------------------

plt.tight_layout()

try:
    plt.show(block=True)
finally:
    csv_path.unlink(missing_ok=True)