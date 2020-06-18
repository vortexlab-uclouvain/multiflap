import numpy as np

frequency = 4

wingframe_position = np.array([0, -0.0, -0.05])
wingframe_position_tail = np.array([0, -0., .3])
tail_length = 0.25
tail_chord = 0.15
tail_opening = 0.

"""Kinematics"""

# =============================================================================
# Shoulder parameters
# =============================================================================
amplitude_shoulder_x = 0.014
offset_shoulder_x = 0.2
phaseangle_shoulder_x = -np.pi/2

amplitude_shoulder_y =  np.pi/12 #np.deg2rad(20)
offset_shoulder_y =  -np.deg2rad(19)
phaseangle_shoulder_y = np.pi/2

amplitude_shoulder_z = np.deg2rad(42)
offset_shoulder_z = 0.
phaseangle_shoulder_z = np.pi

# =============================================================================
# Elbow parameters
# =============================================================================
offset_elbow_y = np.pi/6
amplitude_elbow_y = np.pi/6
phaseangle_elbow_y = -np.pi/2

offset_elbow_x = 0.
amplitude_elbow_x = np.pi/6
phaseangle_elbow_x = -np.pi/2

# =============================================================================
# Wrist parameters
# =============================================================================
offset_wrist_y = -np.pi/6
amplitude_wrist_y = np.pi/6
phaseangle_wrist_y = np.pi/2

offset_wrist_z = 0.
amplitude_wrist_z = 0.
phaseangle_wrist_z = 0.