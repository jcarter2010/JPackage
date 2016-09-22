import numpy as np
import sys
import os
import math
from JMath import Vector3
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib
import random
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
matplotlib.use('TkAgg')
from matplotlib.figure import Figure
import JPAstro
if os.name == 'posix':
    if sys.version_info[0] == 2:
        import Tkinter as Tk
    else:
        import tkinter as Tk
if os.name == 'nt':
    import tkinter as Tk

class Main:

    global fig
    global canvas
    global ax
    global astro_eqs
    global particles
    global time_step
    global distance_scale

    def __init__(self):
        root = Tk.Tk()
        root.wm_title('Dark Matter Simulation')
        self.fig = Figure(figsize=(8, 8), dpi=100)
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.canvas = FigureCanvasTkAgg(self.fig, master=root)
        self.canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
        self.canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
        self.astro_eqs = JPAstro.Equations()
        self.time_step = math.pow(10, 13)
        self.distance_scale = math.pow(10, 0)
        positions = []
        masses = []
        forces = []
        velocities = []
        for x in range(-120, 121, 40):
            for y in range(-120, 121, 40):
                for z in range(-120, 121, 40):
                    offset_x = (random.random() - .5) * 20
                    offset_y = (random.random() - .5) * 20
                    offset_z = (random.random() - .5) * 20
                    pos = Vector3(x + offset_x, y + offset_y, z + offset_z)
                    positions.append(pos)
                    masses.append(math.pow(10, -1))
                    forces.append(Vector3(0, 0, 0))
                    velocities.append(Vector3(0, 0, 0))

        data = {'Position' : positions, 'Mass' : masses, 'Force' : forces, 'Velocity' : velocities}

        self.particles = pd.DataFrame(data)

        self.Update()

    def Update(self):
        while True:
            self.Draw()
            print('Updating Gravity')
            self.Update_Gravity()
            print('Updating Velocities')
            self.Update_Velocities()
            print('Updating Positions')
            self.Update_Positions()
            print('Drawing')


    def Update_Gravity(self):
        for i in range(0, len(self.particles['Position']) - 1):
            for j in range(i + 1, len(self.particles['Position'])):
                if not i == j:
                    m_i = self.particles['Mass'][i]
                    m_j = self.particles['Mass'][j]
                    p_i = self.particles['Position'][i]
                    p_j = self.particles['Position'][j]
                    R = p_i.Distance(p_j)
                    direction = p_j.Subtract(p_i).Normalize()
                    F = self.astro_eqs.Solve_Equation('F = (G M m) / (R^2)', 'F', m_i, m_j, R)
                    grav_i = direction.Scalar_Mult(F)
                    grav_j = grav_i.Scalar_Mult(-1)
                    self.particles.set_value(i, 'Force', self.particles['Force'][i].Add(grav_i))
                    self.particles.set_value(j, 'Force', self.particles['Force'][j].Add(grav_j))

    def Update_Velocities(self):
        for i in range(0, len(self.particles['Position'])):
            self.particles.set_value(i, 'Velocity', self.particles['Velocity'][i].Add(self.particles['Force'][i].Scalar_Mult(self.time_step)))

    def Update_Positions(self):
        for i in range(0, len(self.particles['Position'])):
            self.particles.set_value(i, 'Position', self.particles['Position'][i].Add(self.particles['Velocity'][i]))

    def Draw(self):
        xs = []
        ys = []
        zs = []
        us = []
        vs = []
        ws = []
        counter = 0
        for pos in self.particles['Position']:
            velocity = self.particles['Velocity'][counter]
            xs.append(pos.x)
            ys.append(pos.y)
            zs.append(pos.z)
            us.append(velocity.x)
            vs.append(velocity.y)
            ws.append(velocity.z)
            counter = counter + 1
        self.ax.clear()
        if sys.version_info[0] == 3:
            self.ax.quiver(xs, ys, zs, us, vs, ws, length = 10, pivot = 'tail', arrow_length_ratio = .75)
        self.ax.scatter(xs, ys, zs, c='r', s=10, marker = 'o')
        self.ax.set_xlim3d(-120, 120)
        self.ax.set_ylim3d(-120, 120)
        self.ax.set_zlim3d(-120, 120)
        self.ax.set_xlabel('X Label')
        self.ax.set_ylabel('Y Label')
        self.ax.set_zlabel('Z Label')
        self.canvas.draw()

    def rgb_to_hex(self, rgb):
        return '#%02x%02x%02x' % rgb

    def print_full(self, x):
        pd.set_option('display.max_rows', len(x))
        print(x)
        pd.reset_option('display.max_rows')

if __name__ == '__main__':
    Main()
