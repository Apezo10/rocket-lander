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
playback_speed = 15  # simulation seconds per GIF second
max_frames = 400


# Read flight data
data = np.genfromtxt(csv_file, delimiter=",", names=True)
data = np.atleast_1d(data)

columns = ["time", "x", "y", "z", "qw", "qx", "qy", "qz"]

if data.dtype.names is None or not all(name in data.dtype.names for name in columns):
    raise ValueError("Trajectory.csv needs time, position, and quaternion columns.")

if len(data) < 2:
    raise ValueError("Trajectory.csv needs at least two data rows.")

time = data["time"]

position = np.column_stack((
    data["x"],
    data["y"],
    data["z"]
))

quaternion = np.column_stack((
    data["qw"],
    data["qx"],
    data["qy"],
    data["qz"]
))

if not np.all(np.isfinite(time)) or not np.all(np.isfinite(position)) or not np.all(np.isfinite(quaternion)):
    raise ValueError("The CSV contains invalid numbers.")

if np.any(np.diff(time) <= 0):
    raise ValueError("CSV times must increase.")

q_length = np.linalg.norm(quaternion, axis=1)

if np.any(q_length < 1e-12):
    raise ValueError("A quaternion in the CSV is zero.")

quaternion = quaternion / q_length[:, None]


# Quaternion functions
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

        q = (
            np.sin((1 - fraction) * angle) * q0
            + np.sin(fraction * angle) * q1
        ) / np.sin(angle)

    return q / np.linalg.norm(q)


def get_rotation(q):

    w, x, y, z = q

    return np.array([
        [1 - 2*(y*y + z*z), 2*(x*y - w*z), 2*(x*z + w*y)],
        [2*(x*y + w*z), 1 - 2*(x*x + z*z), 2*(y*z - w*x)],
        [2*(x*z - w*y), 2*(y*z + w*x), 1 - 2*(x*x + y*y)]
    ])


# Build a cylinder along the rocket's body X-axis
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

    side_faces.append([
        tail_ring[i],
        tail_ring[j],
        nose_ring[j],
        nose_ring[i]
    ])

    nose_faces.append([
        [length / 2, 0, 0],
        nose_ring[i],
        nose_ring[j]
    ])

    tail_faces.append([
        [-length / 2, 0, 0],
        tail_ring[j],
        tail_ring[i]
    ])

side_faces = np.array(side_faces)
nose_faces = np.array(nose_faces)
tail_faces = np.array(tail_faces)


# A stripe makes rolling visible on the cylinder
stripe = np.array([
    [-length / 2, radius + 0.15, 0],
    [length / 2, radius + 0.15, 0]
])


def add_rocket(ax):

    body = Poly3DCollection(
        side_faces,
        facecolor="silver",
        edgecolor="gray",
        linewidth=0.3
    )

    nose = Poly3DCollection(
        nose_faces,
        facecolor="black"
    )

    tail = Poly3DCollection(
        tail_faces,
        facecolor="dimgray"
    )

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


# Animation time and trajectory
flight_time = time[-1] - time[0]

frame_count = min(
    max_frames,
    max(2, int(flight_time * fps / playback_speed) + 1)
)

frame_times = np.linspace(time[0], time[-1], frame_count)

trajectory = position / 1000  # km

frame_positions = np.column_stack([
    np.interp(frame_times, time, position[:, i]) / 1000
    for i in range(3)
])


# Set up the two 3D views
fig = plt.figure(figsize=(13, 6.5))

ax_path = fig.add_subplot(121, projection="3d")
ax_body = fig.add_subplot(122, projection="3d")

fig.subplots_adjust(
    left=0.02,
    right=0.98,
    top=0.87,
    bottom=0.07,
    wspace=0.08
)

ax_path.plot(
    trajectory[:, 0],
    trajectory[:, 1],
    trajectory[:, 2],
    color="lightgray"
)

trail, = ax_path.plot([], [], [], color="blue", linewidth=2)

path_rocket = add_rocket(ax_path)
body_rocket = add_rocket(ax_body)


# Keep the trajectory axes fixed as the rocket moves
span = np.ptp(trajectory, axis=0)
max_span = max(np.max(span), 0.1)

center = (
    np.max(trajectory, axis=0)
    + np.min(trajectory, axis=0)
) / 2

half_range = np.maximum(span / 2, 0.2 * max_span)

ax_path.set_xlim(center[0] - half_range[0], center[0] + half_range[0])
ax_path.set_ylim(center[1] - half_range[1], center[1] + half_range[1])
ax_path.set_zlim(center[2] - half_range[2], center[2] + half_range[2])

ax_path.set_box_aspect(2 * half_range)
ax_path.view_init(elev=20, azim=-63)

ax_path.set_title("Rocket trajectory (rocket enlarged)")

ax_path.set_xlabel("X (km)")
ax_path.set_ylabel("Y (km)")
ax_path.set_zlabel("Altitude (km)")


# The close-up uses the cylinder's real size and proportions
zoom = 0.6 * length

ax_body.set_xlim(-zoom, zoom)
ax_body.set_ylim(-zoom, zoom)
ax_body.set_zlim(-zoom, zoom)

ax_body.set_box_aspect((1, 1, 1))
ax_body.view_init(elev=23, azim=-55)

ax_body.set_title("Rocket orientation (true dimensions)")

ax_body.set_xlabel("X (m)")
ax_body.set_ylabel("Y (m)")
ax_body.set_zlabel("Z (m)")


# Only the trajectory-view rocket is enlarged so it can be seen from far away
path_scale = 0.13 * max_span / length

actual_speed = flight_time * fps / (frame_count - 1)

title = fig.suptitle("")


def update(frame):

    rotation = get_rotation(get_quaternion(frame_times[frame]))
    rocket_pos = frame_positions[frame]

    move_rocket(path_rocket, rotation, rocket_pos, path_scale)
    move_rocket(body_rocket, rotation, np.zeros(3))

    trail.set_data(
        frame_positions[:frame + 1, 0],
        frame_positions[:frame + 1, 1]
    )

    trail.set_3d_properties(
        frame_positions[:frame + 1, 2]
    )

    title.set_text(
        f"Time: {frame_times[frame]:.1f} s   |   "
        f"Altitude: {rocket_pos[2]:.2f} km   |   "
        f"Playback: {actual_speed:.1f}x"
    )


# Generate GIF
animation = FuncAnimation(
    fig,
    update,
    frames=frame_count,
    interval=1000 / fps,
    blit=False
)

print(f"Creating {gif_file.name} ({frame_count} frames)...")

animation.save(
    gif_file,
    writer=PillowWriter(fps=fps),
    dpi=85
)

plt.close(fig)

print(f"Saved: {gif_file}")