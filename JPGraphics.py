from tkinter import *
import math
from JMath import Matrix
from JMath import Vector3

class Camera:

    global position
    global theta
    global phi
    global r
    global view_matrix

    def __init__(self, val1, val2, val3):
        self.position = Vector3(val1, val2, val3)
        self.theta = 0
        self.phi = 0
        self.r = 20
        self.view_matrix = Matrix([[0,0,0,0,], [0,0,0,0], [0,0,0,0], [0,0,0,0]])

    def Set_View_Matrix(self):
        self.view_matrix = Matrix.Create_From_Yaw_Pitch_Roll(self.phi, self.theta, 0)

    def Rotate_Camera_Around_Point(self, origin, yaw, pitch, roll):
    	self.position.Rotate_Around_Point(origin, yaw, pitch, roll)
    	self.phi -= yaw
    	self.theta -= math.cos(self.phi) * pitch
    	self.theta -= math.sin(self.phi) * roll

class Graphics_Controller:

    global camera
    global projection_matrix
    global width
    global height
    global translation_matrix
    global view_projection_transform

    def __init__(self, w, h):
        self.width = w
        self.height = h
        fov = 1
        H = math.pi / 4
        V = math.pi / 4
        aspectRatio = math.tan(H / 2) / math.tan(V / 2)
        zFar = 30
        zNear = 1
        q = zFar / (zFar - zNear)
        self.projection_matrix = Matrix([[aspectRatio, 0, 0, 0], [0, fov, 0, 0], [0, 0, q, 1], [0, 0, -q*zNear, 0]])
        self.camera = Camera(0, 0, -450)
        self.translation_matrix = Matrix([[1, 0, 0, -self.camera.position.x], [0, 1, 0, -self.camera.position.y], [0, 0, 1, -self.camera.position.z], [0, 0, 0, 1]]);
        self.view_projection_transform = Matrix.Multiply(self.camera.view_matrix, self.translation_matrix)
        self.root = Tk()
        self.root.bind('<KeyPress>', self.On_Key_Press)
        self.canvas = Canvas(self.root, width = self.width, height = self.height)
        self.canvas.pack()
        print(self.width)

    def On_Key_Press(self, event):
        if event.char == 'a':
            self.camera.Rotate_Camera_Around_Point(Vector3(0, 0, 0), -.1, 0, 0)
        if event.char == 'd':
            self.camera.Rotate_Camera_Around_Point(Vector3(0, 0, 0), .1, 0, 0)
        if event.char == 'w':
            self.camera.Rotate_Camera_Around_Point(Vector3(0, 0, 0), 0, .1, 0)
        if event.char == 's':
            self.camera.Rotate_Camera_Around_Point(Vector3(0, 0, 0), 0, -.1, 0)

    def Prep_Draw(self):
        self.canvas.delete('all')
        self.camera.Set_View_Matrix()
        self.translation_matrix = Matrix([[1, 0, 0, -self.camera.position.x], [0, 1, 0, -self.camera.position.y], [0, 0, 1, -self.camera.position.z], [0, 0, 0, 1]]);
        self.view_projection_transform = Matrix.Multiply(self.camera.view_matrix, self.translation_matrix)

    def Update_Canvas(self):
        self.canvas.update()

class Axis:

    global start_position
    global end_position
    global color

    def __init__(self, s_pos, e_pos, c):
        self.start_position = s_pos
        self.end_position = e_pos
        self.color = c

    def Rotate(self, projection_matrix, view_projection_transform):
        self.start_position.Rotate(projection_matrix, view_projection_transform)
        self.end_position.Rotate(projection_matrix, view_projection_transform)

    def Set_Display_Values(self, width, height):
        self.start_position.Set_Display_Values(width, height)
        self.end_position.Set_Display_Values(width, height)

    def Draw(self, canvas):
        canvas.create_line(self.start_position.display_x, self.start_position.display_y, self.end_position.display_x, self.end_position.display_y, fill=self.color)
