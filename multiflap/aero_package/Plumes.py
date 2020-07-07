import numpy as np

def plumes():
    lp1 = .25;
    lp2 = .225;
    lp3 = .2;
    lp4 = .175;
    lp5 = .15;
    lp6 = .125;
    lp7 = .1;
    pY1 = np.pi/2;
    pY2 = 3*np.pi/8;
    pY3 = np.pi/4;
    pY4 = np.pi/8;
    pY5 = 0;
    pY6 = 0;
    pY7 = 0;
    pX1 = 0;
    pX2 = 0;
    pX3 = 0;
    pX4 = 0;
    pX5 = 0;
    pX6 = 0;
    pX7 = 0;
    prot = np.array([[pY1, pY2, pY3, pY4, pY5, pY6, pY7],
            [pX1, pX2, pX3, pX4, pX5, pX6, pX7]])
    pvec = np.array([[0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0],
            [lp1, lp2, lp3, lp4, lp5, lp6, lp7]])
    return prot, pvec

if __name__ == "__main__":
   plumes()

