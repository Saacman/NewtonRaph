## slider-crank linkage
## Newton - Raphson 
## Isaac Sanchez - UGTO

import numpy as np
from math import cos, sin, pi, sqrt

# Parameters of the mechanism
e  = 0
a2 = 3
a3 = 8
theta_e = 0
theta_2 = (40 * pi) / 180

#Set initial Values
theta3_s = np.array([[(300 * pi) / 180], [5]])

# Calculate the Jacobian of the system
# Measure the angle from the positive x-axis
def createJacobian(t_3):
    J = np.matrix([[ -a3 * sin(t_3), -1],
                    [ a3 * cos(t_3),  0]])
    return J

# Evaluate the position functions, given the values
def posFunc(t_3, s):
    f1 = a2 * cos(theta_2) + a3 * cos(t_3) - s
    f2 = a2 * sin(theta_2) + a3 * sin(t_3) +e
    return np.matrix([[f1], [f2]])


def updateValues(thetas):
    f = posFunc(thetas[0][0], thetas[1][0])
    print(f"Position Functions:\nf_theta3: {f[0]}\nf_theta4: {f[1]}")
    J = createJacobian(thetas[0][0])
    print(f"Jacobian:\n{J}")
    newF = thetas + np.linalg.inv(J) * (-f)
    return newF


if __name__ == '__main__':
    
    for i in range(5):
        print(f"Iteration {i}\n")
        old = theta3_s
        theta3_s = updateValues(theta3_s)
        print(f" Theta_3: {360 - (theta3_s[0][0] * 180)/ pi}\ns: {theta3_s[1][0]}")
        print("")
        if abs(old[0][0]-theta3_s[0][0]) < 0.00001:
            break
