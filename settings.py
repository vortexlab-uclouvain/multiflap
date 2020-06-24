import numpy as np


class SimulationsSettings():

    def __init__(self):

        self.frequency = 4
        self.wingframe_position = np.array([0, -0.0, -0.05])
        self.wingframe_position_tail = np.array([0, -0., .3])
        self.tail_length = 0.25
        self.tail_chord = 0.15
        self.tail_opening = 0.
