import math

class Vector2:

    global x
    global y
    global rotate_x
    global rotate_y
    global rotate_z
    global rotate_y

    def GetX(self):
        return x

    def GetY(self):
        return y

    def SetX(self, val):
        self.x = val

    def SetY(self, val):
        self.y = val

    def __init__(self, val1, val2):
        self.x = val1
        self.y = val2

    def From_Spherical(self, r, phi, theta):
        self.x = r * math.sin(theta) * math.cos(phi)
        self.y = r * math.sin(theta) * math.sin(phi)
        self.z = r * math.cos(theta)

    def From_Cylindrical(self, s, phi, z):
        self.x = s * math.cos(phi)
        self.y = s * math.sin(phi)
        self.z = z

    def Get_Spherical(self):
        r = math.sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)
        phi = math.atan(self.y / self.x)
        theta = math.acos(self.x / r)
        return (r, phi, theta)

    def Get_Cylindrical(self):
        r = math.sqrt(self.x ** 2 + self.y ** 2)
        phi = math.atan(self.y / self.x)
        z = self.z
        return (r, phi, z)

    def Rotate_Around_Point(self, origin, theta):
        point = [[self.x - origin.x], [self.y - origin.y], [0], [1]]
        rotation_matrix = Matrix.Create_From_Yaw_Pitch_Roll(0, 0, theta)
        resultant = Matrix.Point_Multiply(point, rotation_matrix)
        self.rotate_x = resultant.matrix[0][0]
        self.rotate_y = resultant.matrix[1][0]
        self.x = self.rotate_x + origin.x
        self.y = self.rotate_y + origin.y

    def Equals(self, v2):
        if self.Subtract(v2).Length() <= .001:
            return True
        else:
            return False

    def Length(self):
        return math.sqrt(self.x * self.x + self.y * self.y)

    def Add(self, v2):
        return Vector2(self.x + v2.x, self.y + v2.y)

    def Normalize(self):
        if self.Length() != 0:
            return self.Scalar_Mult(1 / self.Length())
        else:
            return self

    def Subtract(self, v2):
        return Vector2(self.x - v2.x, self.y - v2.y)

    def Scalar_Mult(self, s):
        return Vector2(self.x * s, self.y * s)

    def Dot(self, v2):
        return self.x * v2.x + self.y * v2.y

    def Cross(self, v2):
        return self.x * v2.y - self.y * v2.x

    def Angle_Between_Vectors(self, v2):
        return math.asin(self.Cross(v2) / (self.Length() * v2.Length()))

    def Distance(self, v2):
        temp_vector = self.Subtract(v2)
        return temp_vector.Length()

    def Get_Slope(self):
        if self.x != 0:
            return(self.y / self.x)
        return(float('inf'))

    def __str__(self):
        return '(' + str(self.x) + ',' + str(self.y) + ')'

    def __repr__(self):
        return self.__str__()

class Vector3:
    global x
    global y
    global z
    global display_x
    global display_y
    global rotate_x
    global rotate_y
    global rotate_z
    global rotate_w
    global allowed

    def GetX(self):
        return x

    def GetY(self):
        return y

    def GetZ(self):
        return z

    def SetX(self, val):
        self.x = val

    def SetY(self, val):
        self.y = val

    def SetZ(self, val):
        self.z = val

    def __init__(self, val1, val2, val3):
        self.x = val1
        self.y = val2
        self.z = val3

    def Rotate(self, projection_matrix, view_projection_transform):
        resultant = Matrix.Point_Multiply([[self.x], [self.y], [self.z], [1]], view_projection_transform)
        resultant2 = Matrix.Multiply(projection_matrix, resultant)
        self.rotate_x = resultant2.matrix[0][0]
        self.rotate_y = resultant2.matrix[1][0]
        self.rotate_z = resultant2.matrix[2][0]
        self.rotate_w = resultant2.matrix[3][0]

    def Rotate_Around_Point(self, origin, yaw, pitch, roll):
        point = [[self.x - origin.x], [self.y - origin.y], [self.z - origin.z], [1]]
        rotation_matrix = Matrix.Create_From_Yaw_Pitch_Roll(yaw, pitch, roll)
        resultant = Matrix.Point_Multiply(point, rotation_matrix)
        self.rotate_x = resultant.matrix[0][0]
        self.rotate_y = resultant.matrix[1][0]
        self.rotate_z = resultant.matrix[2][0]
        self.x = self.rotate_x + origin.x
        self.y = self.rotate_y + origin.y
        self.z = self.rotate_z + origin.z

    def Equals(self, v2):
        if self.Subtract(v2).Length() <= .001:
            return True
        else:
            return False

    def Set_Display_Values(self, width, height):
        self.display_x = int((self.rotate_x * width) / (2.0 * self.rotate_w) + (0.5 * width))
        self.display_y = int((self.rotate_y * height) / (2.0 * self.rotate_w) + (0.5 * height))

    def Length(self):
        return math.sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

    def Add(self, v2):
        return Vector3(self.x + v2.x, self.y + v2.y, self.z + v2.z)

    def Normalize(self):
        if self.Length() != 0:
            return self.Scalar_Mult(1 / self.Length())
        else:
            return self

    def Subtract(self, v2):
        return Vector3(self.x - v2.x, self.y - v2.y, self.z - v2.z)

    def Scalar_Mult(self, s):
        return Vector3(self.x * s, self.y * s, self.z * s)

    def Dot(self, v2):
        return self.x * v2.x + self.y * v2.y + self.z * v2.z

    def Cross(self, v2):
        return Vector3(self.y * v2.z - self.z * v2.y, self.z * v2.x - self.x * v2.z, self.x * v2.y - self.y * v2.x)

    def Distance(self, v2):
        temp_vector = self.Subtract(v2)
        return temp_vector.Length()

    def __str__(self):
        return '(' + str(self.x) + ',' + str(self.y) + ',' + str(self.z) + ')'

    def __repr__(self):
        return self.__str__()

class Matrix:

    global matrix
    global identity_matrix

    def __init__(self):
        self.matrix = [[0,0,0,0,], [0,0,0,0], [0,0,0,0], [0,0,0,0]]
        self.identity_matrix = [[1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1]]

    def __init__(self, val):
        self.matrix = val
        self.identity_matrix = [[1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1]]

    def Multiply(matrix_a, matrix_b):
        A = matrix_a.matrix
        B = matrix_b.matrix
        a_rows = len(A)
        a_cols = len(A[0])
        b_rows = len(B)
        b_cols = len(B[0])
        sum_val = 0
        temp_array = []
        resultant_array = []
        for i in range(0, a_rows):
            resultant_array.append([])
            for j in range(0, b_cols):
                for k in range(0, b_rows):
                    sum_val = sum_val + A[i][k] * B[k][j]
                resultant_array[i].append(sum_val)
                sum_val = 0
        resultant_matrix = Matrix(resultant_array)
        return resultant_matrix

    def Create_From_Look_At(position, look_at_position):
        look_vector = look_at_position.subtract(position)
        r = look_vector.Length()
        theta = math.acos(look_vector.y / r)
        phi = math.atan(look_vector.z / look_vector.x)
        resultant = Create_From_Yaw_Pitch_Roll(phi, theta, 0)
        return resultant

    def Create_From_Yaw_Pitch_Roll(yaw, pitch, roll):
        rotation_z_axis = Matrix([[1,0,0,0], [0, math.cos(pitch), -math.sin(pitch), 0], [0, math.sin(pitch), math.cos(pitch), 0], [0, 0, 0, 1]])
        rotation_y_axis = Matrix([[math.cos(yaw), 0, math.sin(yaw), 0], [0, 1, 0, 0], [-math.sin(yaw), 0, math.cos(yaw), 0], [0, 0, 0, 1]])
        rotation_x_axis = Matrix([[math.cos(roll), -math.sin(roll), 0, 0], [math.sin(roll), math.cos(roll), 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
        resultant = Matrix.Multiply(Matrix.Multiply(rotation_x_axis, rotation_y_axis), rotation_z_axis)
        return resultant

    def Point_Multiply(point, mult_matrix):
        resultant = Matrix.Multiply(mult_matrix, Matrix(point))
        return resultant
