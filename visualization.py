from subprocess import Popen, PIPE
import sys
import filecmp
import shutil
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

import math


def run_cpp():
    points = []
    cpp_prog = Popen('build/aks_project', shell=True, stdout=PIPE, stdin=PIPE)
    result = cpp_prog.stdout.readlines()

    for i in range(len(result) - 1):
        l = str(result[i]).split()
        print (l)
        points.append(float(l[2][:-1]))
        # ---- 3d -----
        # points.append((float(l[2][:-1]), float(l[3][:-1])))
    return points


def f3(X, Y):
    return np.double(
        -20 * math.exp(-0.2 * math.pow(0.5 * math.pow(X, 2) + math.pow(Y, 2), 0.5)) \
        - math.exp(0.5 * (math.cos(2 * math.pi * X)
                          + math.cos(2 * math.pi * Y))) + math.e + 20)


def f(X):
    return np.double(-(1.4 - 3 * X) * math.sin(18 * X))


def three_dim_plot(points):
    for i in range(len(points)):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        x = np.array(np.linspace(-0.5, 0.5, 100))
        y = np.array(np.linspace(-0.5, 0.5, 100))
        X, Y = np.meshgrid(x, y)

        F = np.vectorize(f3)

        # Z =  -(Y + 47) * math.sin(math.sqrt(abs(X / 2 + (Y + 47)))) - X * math.sin(math.sqrt(abs(X - (Y + 47))))

        ax.plot_surface(X, Y, F(X, Y), rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)

        ax.scatter3D(points[i][0], points[i][0], F(points[i][0], points[i][1]), s=80, c="limegreen");

        plt.savefig('img/img3d{}.png'.format(str(i)))

        plt.show()


if __name__ == "__main__":
    points = run_cpp()
    
#     three_dim_plot(points)

    for i in range(len(points)):
        fig = plt.figure()
        ax = plt.subplot(1, 1, 1)
        X = np.array(np.linspace(0, 1.2))
        F = np.vectorize(f)
        ax.plot(X, F(X))
        ax.scatter(points[i], F(points[i]))
        plt.savefig('img/img{}.png'.format(str(i)))
        plt.show()

