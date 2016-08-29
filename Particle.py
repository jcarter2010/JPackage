from JMath import Vector3
from tkinter import *
import JPAstro
import random
import math


class Particle:

    global position
    global velocity
    global index
    global mass
    global scale_factor

    def __init__(self, i, m, d):
        self.position = Vector3((random.random() - .5) * d, (random.random() - .5) * d, (random.random() - .5) * d)
        self.velocity = Vector3(0, 0, 0)
        #self.velocity = Vector3(random.random() - .5, random.random() - .5, random.random() - .5)
        self.index = i
        self.mass = m
        self.scale_factor = .001

    def Update_Velocity(self, particles, distance_scale, time_step):
        for i in range(0, len(particles)):
            r = self.position.Subtract(particles[i].position).Length() * distance_scale
            if particles[i].position.Equals(self.position) == False:
                self.velocity = self.velocity.Add(particles[i].position.Subtract(self.position).Normalize().Scalar_Mult(JPAstro.Constants.G * particles[i].mass * self.mass / math.pow(r, 2)).Scalar_Mult(time_step))

    def Update_Position(self):
        self.position = self.position.Add(self.velocity)

    def Rotate(self, projection_matrix, view_projection_transform):
        self.position.Rotate(projection_matrix, view_projection_transform)

    def Set_Display_Values(self, width, height):
        self.position.Set_Display_Values(width, height)

    def Draw_Particle(self, canvas, camera_position):
        size = 1 / (camera_position.Subtract(self.position).Length()) * self.mass * self.scale_factor
        canvas.create_oval(self.position.display_x - size, self.position.display_y - size, self.position.display_x + size, self.position.display_y + size)
