## 4-Bar Linkage
## Newton - Raphson 
## Isaac Sanchez - UGTO

import numpy as np
from math import cos, sin, pi, sqrt

# Parameters of the mechanism
a1 = 1
a2 = 1
a3 = 2
a4 = sqrt(2)
theta1 = 0
theta2 = (90 * pi) / 180

#Set initial Values
theta3_4 = np.array([[(20 * pi) / 180], [(30* pi) / 180]])

# Calculate the Jacobian of the system
# Measure the angle from the positive x-axis
def createJacobian(t_3, t_4):
    J = np.array([[ -a3 * sin(t_3), a4 * sin(t_4)],
                    [ a3 * cos(t_3), -a4 * cos(t_4)]])
    return J

# Evaluate the position functions, given the values
def posFunc(t_3, t_4):
    f1 = a2 * cos(theta2) + a3 * cos(t_3) -  a1 * cos(theta1) - a4 * cos(t_4)
    f2 = a2 * sin(theta2) + a3 * sin(t_3) -  a1 * sin(theta1) - a4 * sin(t_4)
    return np.array([[f1], [f2]])


def updatePosValues(thetas):
    f = posFunc(thetas[0][0], thetas[1][0])
    print(f"Position Functions:\nf_theta3: {f[0]}\nf_theta4: {f[1]}")
    J = createJacobian(thetas[0][0], thetas[1][0])
    print(f"Jacobian:\n{J}")
    newF = thetas + np.dot(np.linalg.inv(J), (-f))
    print(f"Inverse Jacobian:\n{np.linalg.inv(J)}")
    return newF

# Evaluate the velocity functions, given the values
# TODO Fix theta3 , theta4
def velFunc(w_3, w_4):
    f1 = -a2 * sin(theta2) * w2 - a3 * sin(theta3) * w_3 + a4 * sin(theta4) * w_4
    f2 = a2 * cos(theta2) * w2 + a3 * cos(theta3) * w_3 - a4 * cos(theta4) * w_4
    return np.array([[f1], [f2]])

#TODO: FIX omegas3_4!!
def updateVelValues(omegas):
    f = velFunc(omegas[0][0], omegas[1][0])
    print(f"Velocity Functions:\nf_theta3: {f[0]}\nf_theta4: {f[1]}")
    J = createJacobian(theta3, theta4)
    print(f"Jacobian:\n{J}")
    newF = omegas + np.dot(np.linalg.inv(J), (-f))
    print(f"Inverse Jacobian:\n{np.linalg.inv(J)}")
    return newF

if __name__ == '__main__':

    # Solve the Position of the system
    for i in range(3):
        print(f"Iteration {i}\n")
        old = theta3_4
        theta3_4 = updatePosValues(theta3_4)
        print(f" Theta_3_4:\n{(theta3_4 * 180)/pi}")
        print("")
        if abs(old[0][0]-theta3_4[0][0]) < 0.00001:
            break

