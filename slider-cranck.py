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
w2 = 209.43951
alpha_2 = 0
#Set initial Values
theta3_s = np.array([[(300 * pi) / 180], [5]])
w3_s = np.array([[0], [0]])
alpha3_s = np.array([[0], [0]])


# Calculate the Jacobian of the system
# Measure the angle from the positive x-axis
def createJacobian(t_3):
    J = np.array([[ -a3 * sin(t_3), -1],
                    [ a3 * cos(t_3),  0]])
    return J

# Evaluate the position functions, given the values
def posFunc(t_3, s):
    f1 = a2 * cos(theta_2) + a3 * cos(t_3) - s
    f2 = a2 * sin(theta_2) + a3 * sin(t_3) +e
    return np.array([[f1], [f2]])


def updateValues(thetas):
    f = posFunc(thetas[0][0], thetas[1][0])
    print(f"Position Functions:\nf_theta3: {f[0]}\nf_theta4: {f[1]}")
    J = createJacobian(thetas[0][0])
    print(f"Jacobian:\n{J}")
    newF = thetas + np.dot(np.linalg.inv(J), (-f))
    return newF

# Evaluate the velocity functions, given the values
def velFunc(w_3, s):
    f1 = -a2 * sin(theta_2) * w2 - a3 * sin(theta3_s[0][0]) * w_3 - s
    f2 = a2 * cos(theta_2) * w2 + a3 * cos(theta3_s[0][0]) * w_3
    return np.array([[f1], [f2]])

def updateVelValues(omegas):
    f = velFunc(omegas[0][0], omegas[1][0])
    print(f"Velocity Functions:\nf_w3: {f[0]}\nf_w4: {f[1]}")
    J = createJacobian(theta3_s[0][0])
    print(f"Jacobian:\n{J}")
    newF = omegas + np.dot(np.linalg.inv(J), (-f))
    print(f"Inverse Jacobian:\n{np.linalg.inv(J)}")
    return newF

# Evaluate the Acceleration functions, given the values
def accFunc(alpha_3, s):
    f1 = -a2 * sin(theta_2) * alpha_2 - a2 * cos(theta_2) * pow(w2, 2) -a3 * sin(theta3_s[0][0]) * alpha_3 -a3 * cos(theta3_s[0][0]) * pow(w3_s[0][0], 2) - s
    f2 = a2 * cos(theta_2) * alpha_2 - a2 * sin(theta_2) * pow(w2, 2) +a3 * cos(theta3_s[0][0]) * alpha_3 -a3 * sin(theta3_s[0][0]) * pow(w3_s[0][0], 2)
    return np.array([[f1], [f2]])

def updateAccValues(alphas):
    f = accFunc(alphas[0][0], alphas[1][0])
    print(f"Acceleration Functions:\nf_alpha3: {f[0]}\nf_alpha4: {f[1]}")
    J = createJacobian(theta3_s[0][0])
    print(f"Jacobian:\n{J}")
    newF = alphas + np.dot(np.linalg.inv(J), (-f))
    print(f"Inverse Jacobian:\n{np.linalg.inv(J)}")
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


    
    # Solve the velocities of the system
    for i in range(5):
        print(f"Iteration {i}\n")
        old = w3_s
        w3_s = updateVelValues(w3_s)
        print(f" w3_s:\n{w3_s}")
        print("")
    

    # Solve the accelerations of the system
    for i in range(10):
        print(f"Iteration {i}\n")
        old = alpha3_s
        alpha3_s = updateAccValues(alpha3_s)
        print(f" Alphas_3_4:\n{alpha3_s}") #rad/s2 // in/s2
        print("")
       
