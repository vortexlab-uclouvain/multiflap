import numpy as np
import settings 
rho = 1.225
kin_viscosity = 1.5e-5
Reynolds = float(15*settings.tail_length/kin_viscosity)

def delta_tail(velocity_vector, tail_span):
    alpha = np.arctan(velocity_vector[1]/velocity_vector[2])
    velocity_module = np.linalg.norm(velocity_vector)
    Force = (np.pi/4)*rho*(velocity_module**2)*alpha*(tail_span**2)
    Lift = Force*np.cos(alpha)
    Cf = 1.328/np.sqrt(Reynolds)
    Df = 0.5*rho*(velocity_module**2)*(tail_span*settings.tail_length/2)*Cf
    Drag_ind = 0.5*Lift*alpha
    Drag = Drag_ind + Df + Force*np.sin(alpha)
    Lat_force = 0.
    return Lat_force, Lift, Drag

def tail_geometry(length, **opening):
    tail_opening=opening.get('tail_opening', np.deg2rad(0))
    tail_span = 2*length*np.tan(tail_opening/2)
    AR_tail = 2*(tail_span/length)
    NP_tail = settings.wingframe_position_tail + np.array([0., 0., (2/3)*length])
    return tail_span, AR_tail, NP_tail


if __name__ == "__main__":
    [b, _, _] = tail_geometry(.25, tail_opening=np.deg2rad(25))
    print(b)
