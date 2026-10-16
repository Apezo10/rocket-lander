from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


# Rocket dimensions and animation settings
folder = Path(__file__).resolve().parent
csv_file = folder / "Trajectory.csv"
gif_file = folder / "RocketOrientation.gif"

length = 124.4       # meters
radius = 9.0 / 2.0    # meters
fps = 20
playback_speed = 20  # simulation seconds per GIF second
max_frames = 400
#---------------------------------------------------------------------------------------------------


# Read flight data
# Positions and velocities in the CSV are now in ECI. Altitude is exported
# separately by C++, since ECI Z is not height above the ground.
data = np.genfromtxt(csv_file, delimiter=",", names=True)
data = np.atleast_1d(data)


columns = ["time", "x", "y", "z", "altitude", "vx", "vy", "vz", "qw", "qx", "qy", "qz"]


if data.dtype.names is None or not all(name in data.dtype.names for name in columns):
    raise ValueError("Missing columns in Trajectory.csv. Re-export with altitude from C++.")

if len(data) < 2:
    raise ValueError("Trajectory.csv needs at least two data rows.")

time = data["time"]

position = np.column_stack((data["x"], data["y"], data["z"]))
quaternion = np.column_stack((data["qw"], data["qx"], data["qy"], data["qz"]))


if not np.all(np.isfinite(time)) or not np.all(np.isfinite(position)) or not np.all(np.isfinite(quaternion)):
    raise ValueError("The CSV contains invalid numbers.")

if not np.all(np.isfinite(data["altitude"])):
    raise ValueError("The CSV contains invalid altitude values.")

if np.any(np.diff(time) <= 0):
    raise ValueError("CSV times must increase.")

q_length = np.linalg.norm(quaternion, axis=1)
if np.any(q_length < 1e-12):
    raise ValueError("A quaternion in the CSV is zero.")
quaternion = quaternion / q_length[:, None]


# Coordinate conversions
# -----------------------------------------------------------------------------
# C++ integrates in ECI, but the plots are easier to read relative to Earth.
# Use the same simplified Earth rotation convention as the C++ code.
EARTH_OMEGA = 7.2921159e-5  # rad/s
WGS84_A = 6378137.0
WGS84_F = 1.0 / 298.257223563
WGS84_E2 = WGS84_F * (2.0 - WGS84_F)


def eci_to_ecef(coords, times):
    # Rotate vector components into Earth-fixed axes. This is not a
    # rotating-frame velocity correction; that is handled separately below.
    theta = EARTH_OMEGA * np.asarray(times)
    c, s = np.cos(theta), np.sin(theta)
    return np.column_stack((
        c * coords[:, 0] + s * coords[:, 1],
        -s * coords[:, 0] + c * coords[:, 1],
        coords[:, 2],
    ))


def geodetic_angles(ecef):
    # Iterate latitude because Earth is modeled as a WGS84 ellipsoid.
    x, y, z = ecef.T
    longitude = np.arctan2(y, x)
    p = np.hypot(x, y)
    latitude = np.arctan2(z, p * (1.0 - WGS84_E2))
    for _ in range(10):
        sinlat = np.sin(latitude)
        N = WGS84_A / np.sqrt(1.0 - WGS84_E2 * sinlat**2)
        latitude = np.arctan2(z + WGS84_E2 * N * sinlat, p)
    return latitude, longitude


def enu_basis(lat, lon):
    # Local east, north, and up unit vectors, expressed in ECEF.
    slat, clat = np.sin(lat), np.cos(lat)
    slon, clon = np.sin(lon), np.cos(lon)
    return np.array([
        [-slon, clon, 0.0],
        [-slat * clon, -slat * slon, clat],
        [clat * clon, clat * slon, slat],
    ])


def eci_to_ecef_rotation(t):
    theta = EARTH_OMEGA * t
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]])


# Convert each recorded position back to Earth-fixed coordinates.
ecef_position = eci_to_ecef(position, time)
latitudes, longitudes = geodetic_angles(ecef_position)
reference_basis = enu_basis(latitudes[0], longitudes[0])
reference_ecef = ecef_position[0]
local_enu = (ecef_position - reference_ecef) @ reference_basis.T

# Use launch-site east/north displacement for the trajectory view.
# Keep the actual WGS84 altitude for Z rather than tangent-plane height.
altitude_km = data["altitude"] / 1000.0
trajectory = np.column_stack((local_enu[:, 0] / 1000.0,
                              local_enu[:, 1] / 1000.0,
                              altitude_km))

# Remove the velocity due to Earth rotation before plotting airspeed.
# This matches the current C++ aerodynamic model (no weather winds yet).
velocity_eci = np.column_stack((data["vx"], data["vy"], data["vz"]))
rotation_velocity = np.column_stack((-EARTH_OMEGA * position[:, 1],
                                     EARTH_OMEGA * position[:, 0],
                                     np.zeros(len(position))))
velocity_air_eci = velocity_eci - rotation_velocity
velocity_air_ecef = eci_to_ecef(velocity_air_eci, time)
air_speed = np.linalg.norm(velocity_air_eci, axis=1)

# Resolve the air-relative velocity into east, north, up at each position.
velocity_enu = np.array([
    enu_basis(lat, lon) @ v for lat, lon, v in
    zip(latitudes, longitudes, velocity_air_ecef)
])

# Quaternion functions
# -----------------------------------------------------------------------------
# The C++ orientation is body-to-ECI, with the rocket pointing along body +X.
def get_quaternion(t):
    i = np.searchsorted(time, t) - 1
    i = int(np.clip(i, 0, len(time) - 2))

    fraction = (t - time[i]) / (time[i + 1] - time[i])
    q0 = quaternion[i]
    q1 = quaternion[i + 1].copy()

    dot = np.dot(q0, q1)

    # Flip if needed so the interpolation takes the shorter rotation
    if dot < 0:
        q1 = -q1
        dot = -dot

    if dot > 0.9995:
        q = (1 - fraction) * q0 + fraction * q1
    else:
        angle = np.arccos(np.clip(dot, -1, 1))
        q = (np.sin((1 - fraction) * angle) * q0 + np.sin(fraction * angle) * q1) / np.sin(angle)

    return q / np.linalg.norm(q)


def get_rotation(q):
    w, x, y, z = q

    return np.array([
        [1 - 2*(y*y + z*z), 2*(x*y - w*z), 2*(x*z + w*y)],
        [2*(x*y + w*z), 1 - 2*(x*x + z*z), 2*(y*z - w*x)],
        [2*(x*z - w*y), 2*(y*z + w*x), 1 - 2*(x*x + y*y)]
    ])


# Build a cylinder along the rocket's body X-axis
# -----------------------------------------------------------------------------
segments = 24
angles = np.linspace(0, 2*np.pi, segments, endpoint=False)

def make_ring(x):
    return np.column_stack((
        np.full(segments, x),
        radius * np.cos(angles),
        radius * np.sin(angles)
    ))


nose_ring = make_ring(length / 2)
tail_ring = make_ring(-length / 2)

side_faces = []
nose_faces = []
tail_faces = []

for i in range(segments):
    j = (i + 1) % segments

    side_faces.append([tail_ring[i], tail_ring[j], nose_ring[j], nose_ring[i]])
    nose_faces.append([[length / 2, 0, 0], nose_ring[i], nose_ring[j]])
    tail_faces.append([[-length / 2, 0, 0], tail_ring[j], tail_ring[i]])

side_faces = np.array(side_faces)
nose_faces = np.array(nose_faces)
tail_faces = np.array(tail_faces)

# A stripe on the side makes roll visible in the close-up view.
stripe = np.array([
    [-length / 2, radius + 0.15, 0],
    [length / 2, radius + 0.15, 0]
])


def add_rocket(ax):
    body = Poly3DCollection(side_faces, facecolor="silver", edgecolor="gray", linewidth=0.3)
    nose = Poly3DCollection(nose_faces, facecolor="red")
    tail = Poly3DCollection(tail_faces, facecolor="dimgray")

    ax.add_collection3d(body)
    ax.add_collection3d(nose)
    ax.add_collection3d(tail)

    line, = ax.plot([], [], [], color="orange", linewidth=3)
    return body, nose, tail, line


def move_rocket(artists, rotation, center, scale=1.0):
    body, nose, tail, line = artists

    body.set_verts((side_faces @ rotation.T) * scale + center)
    nose.set_verts((nose_faces @ rotation.T) * scale + center)
    tail.set_verts((tail_faces @ rotation.T) * scale + center)

    stripe_points = (stripe @ rotation.T) * scale + center
    line.set_data(stripe_points[:, 0], stripe_points[:, 1])
    line.set_3d_properties(stripe_points[:, 2])

#------------------------------------------------------------------------------------------------------------
# Animation time and trajectory
flight_time = time[-1] - time[0]
frame_count = min(max_frames, max(2, int(flight_time * fps / playback_speed) + 1))
frame_times = np.linspace(time[0], time[-1], frame_count)

frame_positions = np.column_stack([
    np.interp(frame_times, time, trajectory[:, i]) for i in range(3)
])


# Set up the two 3D views
# Left: ground-relative flight path. Right: rocket orientation at full scale.
fig = plt.figure(figsize=(13, 6.5))
ax_path = fig.add_subplot(121, projection="3d")
ax_body = fig.add_subplot(122, projection="3d")
fig.subplots_adjust(left=0.02, right=0.98, top=0.87, bottom=0.07, wspace=0.08)


ax_path.plot(trajectory[:, 0], trajectory[:, 1], trajectory[:, 2], color="lightgray")
trail, = ax_path.plot([], [], [], color="blue", linewidth=2)
path_rocket = add_rocket(ax_path)
body_rocket = add_rocket(ax_body)


# Keep the trajectory axes fixed as the rocket moves
span = np.ptp(trajectory, axis=0)
max_span = max(np.max(span), 0.1)
center = (np.max(trajectory, axis=0) + np.min(trajectory, axis=0)) / 2
half_range = np.maximum(span / 2, 0.2 * max_span)


ax_path.set_xlim(center[0] - half_range[0], center[0] + half_range[0])
ax_path.set_ylim(center[1] - half_range[1], center[1] + half_range[1])
ax_path.set_zlim(center[2] - half_range[2], center[2] + half_range[2])
ax_path.set_box_aspect(2 * half_range)
ax_path.view_init(elev=20, azim=-63)
ax_path.set_title("Rocket trajectory (rocket enlarged)")
ax_path.set_xlabel("East from start (km)")
ax_path.set_ylabel("North from start (km)")
ax_path.set_zlabel("Altitude (km)")

# The close-up uses the cylinder's real size and proportions
zoom = 0.6 * length

ax_body.set_xlim(-zoom, zoom)
ax_body.set_ylim(-zoom, zoom)
ax_body.set_zlim(-zoom, zoom)
ax_body.set_box_aspect((1, 1, 1))
ax_body.view_init(elev=23, azim=-55)
ax_body.set_title("Rocket orientation (true dimensions)")
ax_body.set_xlabel("East (m)")
ax_body.set_ylabel("North (m)")
ax_body.set_zlabel("Up (m)")

# Only the trajectory-view rocket is enlarged so it can be seen from far away
path_scale = 0.13 * max_span / length
actual_speed = flight_time * fps / (frame_count - 1)
title = fig.suptitle("")


def update(frame):
    t = frame_times[frame]

    # Take the body orientation out of ECI and into launch-site ENU.
    rotation = reference_basis @ eci_to_ecef_rotation(t) @ get_rotation(get_quaternion(t))
    rocket_pos = frame_positions[frame]

    move_rocket(path_rocket, rotation, rocket_pos, path_scale)
    move_rocket(body_rocket, rotation, np.zeros(3))

    trail.set_data(frame_positions[:frame + 1, 0], frame_positions[:frame + 1, 1])
    trail.set_3d_properties(frame_positions[:frame + 1, 2])

    title.set_text(
        f"Time: {frame_times[frame]:.1f} s   |   "
        f"Altitude: {rocket_pos[2]:.2f} km   |   "
        f"Playback: {actual_speed:.1f}x"
    )
#------------------------------------------------------------------------------------------------------------



#------------------------------------------------------------------------------------------------------------
# Plot flight and attitude data
def plotFlightData():

    needed = ("vx", "vy", "vz", "p", "q", "r")

    if not all(name in data.dtype.names for name in needed):
        raise ValueError("CSV needs vx, vy, vz, p, q, and r for the plots.")

    # Plot local air-relative velocity rather than inertial XYZ components.
    ve, vn, vu = velocity_enu.T
    speed = air_speed
    altitude = altitude_km

    # Compare the body +X direction (nose) with local up.
    # Comparing against ECI Z would give the wrong tilt outside the poles.
    nose_eci = np.array([get_rotation(q)[:, 0] for q in quaternion])

    nose_ecef = eci_to_ecef(nose_eci, time)

    up_ecef = np.array([enu_basis(lat, lon)[2] for lat, lon in zip(latitudes, longitudes)])

    nose_up = np.einsum("ij,ij->i", nose_ecef, up_ecef)

    tilt = np.degrees(np.arccos(np.clip(nose_up, -1.0, 1.0)))

    # The simulation stores body angular rates in radians per second.
    roll_rate = np.degrees(data["p"])
    pitch_rate = np.degrees(data["q"])
    yaw_rate = np.degrees(data["r"])

    fig, axes = plt.subplots(3, 2, figsize=(14, 12))

    axes[0, 0].plot(time, altitude, linewidth=2)
    axes[0, 0].set_title("Altitude vs. Time")
    axes[0, 0].set_ylabel("Altitude (km)")

    axes[0, 1].plot(time, ve, label="East")
    axes[0, 1].plot(time, vn, label="North")
    axes[0, 1].plot(time, vu, label="Up")
    axes[0, 1].set_title("Air-Relative Velocity (Local ENU)")
    axes[0, 1].set_ylabel("Velocity (m/s)")
    axes[0, 1].legend()

    axes[1, 0].plot(time, speed, linewidth=2)
    axes[1, 0].set_title("Air-Relative Speed vs. Time")
    axes[1, 0].set_ylabel("Air-relative speed (m/s)")

    axes[1, 1].plot(time, tilt, linewidth=2, label="Tilt from upward vertical")
    axes[1, 1].axhline(90, linestyle="--", linewidth=1, color="gray")
    axes[1, 1].set_ylim(0, 180)
    axes[1, 1].set_yticks(np.arange(0, 181, 30))
    axes[1, 1].set_title("Rocket Tilt vs. Time")
    axes[1, 1].set_ylabel("Tilt (degrees)")

    axes[2, 0].plot(time, roll_rate, label="p (roll rate)")
    axes[2, 0].plot(time, pitch_rate, label="q (pitch rate)")
    axes[2, 0].plot(time, yaw_rate, label="r (yaw rate)")
    axes[2, 0].set_title("Body Angular Velocity vs. Time")
    axes[2, 0].set_ylabel("Angular rate (deg/s)")
    axes[2, 0].legend()

    axes[2, 1].axis("off")
    axes[2, 1].text(
        0.05, 0.75,
        "Tilt = 0°: nose up\nTilt = 90°: horizontal\nTilt = 180°: nose down\n\n"
        "p, q, r are angular rates, not Euler angles.",
        transform=axes[2, 1].transAxes, fontsize=12, va="top"
    )

    for ax in axes.flat:
        if ax.axison:
            ax.set_xlabel("Time (s)")
            ax.grid(True, alpha=0.4)

    fig.tight_layout()
    fig.savefig(folder / "FlightData.png", dpi=220)
    plt.close(fig)
    print("Saved: FlightData.png")
#------------------------------------------------------------------------------------------------------------

# PLOT FLIGHT DATA
plotFlightData()

# GENERATE GIF
animation = FuncAnimation(fig, update, frames=frame_count, interval=1000 / fps, blit=False)

print(f"Creating {gif_file.name} ({frame_count} frames)...")
animation.save(gif_file, writer=PillowWriter(fps=fps), dpi=85)
plt.close(fig)

print(f"Saved: {gif_file}")