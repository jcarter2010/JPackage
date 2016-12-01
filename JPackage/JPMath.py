import math
import random
import os

'''
This is not yet finished
It will eventually differentiate a function and return the equation
that is obtained from differentiating it
'''

class Differentiation:

    global eq
    global parts
    global priorities
    global operators

    def __init__(self, eq):
        self.eq = eq
        self.differentiate_parens(eq)

    def differentiate_parens(self, eq):
        parts = eq.replace('(', '---').replace(')', '---').split('---')
        parts = list(filter(None, parts))
        print(parts)
        for i in range(0, len(parts)):
            part = parts[i]
            parts[i] = self.differentiate_exps(part)
            parts[i] = self.check_for_constant(parts[i])
        print(parts)
        return ' '.join(parts)

    def check_for_constant(self, eq):
        parts = eq.split(' ')
        parts = list(filter(None, parts))
        i = 0
        while i < len(parts):
            if parts[i].replace('.', '', 1).isdigit():
                if i == 0:
                    if parts[i + 1] in ['+', '-']:
                        parts.pop(i)
                elif i < len(parts) - 1:
                    if parts[i - 1] in ['+', '-'] and parts[i + 1] in ['+', '-']:
                        parts.pop(i - 1)
                        parts.pop(i - 1)
                else:
                    if parts[i - 1] in ['+', '-']:
                        parts.pop(i - 1)
                        parts.pop(i - 1)
            i = i + 1
        return ' '.join(parts)

    def differentiate_mult(self, eq):
        parts = eq.split('*')
        parts = list(filter(None, parts))
        for i in range(0, len(parts)):
            parts[i] = parts[i].strip()

    def differentiate_exps(self, eq):
        parts = eq.split(' ')
        i = 0
        while i < len(parts):
            if parts[i] == '**':
                parts[i] = self.exponent(' '.join(parts[i - 1:i + 2]))
                parts.pop(i - 1)
                parts.pop(i)
                i = i - 2
            i = i + 1
        return(' '.join(parts))

    def exponent(self, eq):
        print(eq)
        parts = eq.split('**')
        for i in range(0, len(parts)):
            parts[i] = parts[i].strip()
        print(parts)
        if parts[1].strip().replace('.', '', 1).isdigit():
            parts[0] = parts[1] + ' * ' + parts[0]
            parts[1] = str(float(parts[1]) - 1)
            return '{} ** {}'.format(parts[0], parts[1])
        else:
            return parts[0]

    '''
    def differentiate(self, eq):
        priority = 0
        parts = eq.split(' ')
        var = False
        adjust = 0
        for i in range(0, len(parts)):
            try:
                part = parts[i + adjust]
                if var:
                    if part == '**':
                        c, e = self.exponent(float(parts[i + adjust + 1]))
                        parts.insert(i + adjust - 1, '*')
                        parts.insert(i + adjust - 1, c)
                        parts[i + adjust + 3] = e
                        adjust = adjust + 2
                else:
                    if i + adjust < len(parts) - 1:
                        if part.replace('.','',1).isdigit() and parts[i + adjust + 1] in ['+', '-'] and parts[i + adjust - 1] in ['+', '-']:
                            parts.pop(i + adjust + 1)
                            parts.pop(i + adjust)
                            parts[i + adjust - 1] = '+'
                            adjust = adjust - 2
                    else:
                        if part.replace('.','',1).isdigit() and parts[i + adjust - 1] in ['+', '-']:
                            parts.pop(i + adjust - 1)
                            parts.pop(i + adjust - 1)
                            parts[i + adjust - 1] = '+'
                            adjust = adjust - 2
                part = parts[i + adjust]

                if not part.replace('.','',1).isdigit() and not part in self.operators:
                    var = True
                else:
                    var = False
            except:
                break
        print(parts)
    '''

    #def exponent(self, exp):
    #    return (str(exp), str(exp - 1))

class MonteCarloIntegration:
    # input equation in the form 'x * y + 6'
    # with spaces between each character
    # so if you have the variables x, y, and xy,
    # you would input 'x * y + xy'
    # if a number can be negative on the right hand side
    # of the equation, then you have to encase it in parenthesis
    # such as '1 - ( x ) ** 2' if you are going from
    # -1 to 1 for x
    # it's just a quirk that I found with the eval function
    # not sure why it needs it to be this way
    # f is the function to integrate
    # n is the number of random points to take
    # init vals is a dictionary with the keys being the variable names
    #    and the vales being the inital values of the variables
    # end_vals is the ending bounds of the funciton in the integration
    #    in the same dictionary format
    # deltas is the dx for x in y = m * x + b

    def Integrate(self, f, n, init_vals, end_vals, deltas):
        vals = init_vals.copy()
        previous_vals = init_vals.copy()
        parts = f.split()
        return_vals = []
        return_vals.append(self.Evaluate(f, init_vals))
        for i, key in enumerate(vals):
            val = vals[key]
            delta = deltas[key]
            val = val + delta
            vals[key] = val
        return_vals.append(self.Evaluate(f, vals))
        can_run = True
        counter = 0
        while can_run:
            for i, key in enumerate(vals):
                if key in end_vals.keys():
                    if vals[key] >= end_vals[key]:
                        can_run = False
                    else:
                        val = vals[key]
                        delta = deltas[key]
                        val = val + delta
                        vals[key] = val
                        return_vals.append(self.Evaluate(f, vals))
                        previous_vals = vals
                        counter = counter + 1
        max_val = max(return_vals)
        min_val = min(return_vals)

        num_above_above_0 = 0
        num_below_above_0 = 0
        num_above_below_0 = 0
        num_below_below_0 = 0

        for i in range(0, n):
            point_return = random.random() * (max_val - min_val) + min_val
            params = vals
            for i, key in enumerate(vals):
                if key in end_vals.keys():
                    params[key] = random.random() * (end_vals[key] - init_vals[key]) + init_vals[key]
            val = self.Evaluate(f, params)
            if val < point_return:
                if point_return > 0:
                    num_above_above_0 = num_above_above_0 + 1
                else:
                    num_above_below_0 = num_above_below_0 + 1
            else:
                if point_return > 0:
                    num_below_above_0 = num_below_above_0 + 1
                else:
                    num_below_below_0 = num_below_below_0 + 1
        if max_val > 0:
            area_or_vol_etc_above_0 = max_val
        else:
            area_or_vol_etc_above_0 = 0
        if min_val < 0:
            area_or_vol_etc_below_0 = min_val
        else:
            area_or_vol_etc_below_0 = 0
        for i in range(0, len(end_vals.keys())):
            area_or_vol_etc_above_0 = area_or_vol_etc_above_0 * abs(end_vals[list(end_vals.keys())[i]] - init_vals[list(end_vals.keys())[i]])
            area_or_vol_etc_below_0 = area_or_vol_etc_below_0 * abs(end_vals[list(end_vals.keys())[i]] - init_vals[list(end_vals.keys())[i]])
        num_tot_above_0 = num_above_above_0 + num_below_above_0
        num_tot_below_0 = num_above_below_0 + num_below_below_0
        if num_above_above_0 != 0 and num_below_below_0 != 0:
            ratio_above_0 = num_below_above_0 / num_tot_above_0
            result_above_0 = area_or_vol_etc_above_0 * ratio_above_0
            ratio_below_0 = num_above_below_0 / num_tot_below_0
            result_below_0 = area_or_vol_etc_below_0 * ratio_below_0
            result = result_above_0 + result_below_0
        else:
            result = 'inf'

        return result

        # f is the function in the same format as in the previous method
        # vals are the current values of the variables in dict format

    def Evaluate(self, f, vals):
        temp = f.split('=')[1].strip()
        parts = temp.split(' ')
        for i, key in enumerate(vals):
            val = vals[key]
            while key in parts:
                index = parts.index(key)
                parts[index] = str(val)
        g = ' '.join(parts)
        return(eval(g))

class VectorN:

    global pts

    def Help(self):
        lines = []
        f_path = os.path.realpath(__file__)
        path = f_path
        if os.name == 'nt':
            path = f_path[:f_path.rfind('\\')]
        else:
            path = f_path[:f_path.rfind('/')]
        with open('{}/docs/VectorN_Docs'.format(path), 'r') as f_in:
            text = f_in.read()
            lines = text.split('\n')
        for line in lines:
            print(line)

    def __init__(self, pts):
        self.pts = pts

    def Add(self, vect):
        if len(self.pts) != len(vect):
            return None
        else:
            result = []
            for i in range(0, len(self.pts)):
                result.append(self.pts[i] + vect[i])
            return(result)

    def Subtract(self, vect):
        if len(self.pts) != len(vect):
            return None
        else:
            result = []
            for i in range(0, len(self.pts)):
                result.append(self.pts[i] - vect[i])
            return(result)

    def Scalar_Mult(self, scalar):
        for i in range(0, len(self.pts)):
            self.pts[i] = self.pts[i] * scalar

    def Dot(self, vect):
        if len(self.pts) != len(vect):
            return None
        else:
            result = 0
            for i in range(0, len(self.pts)):
                result = result + self.pts[i] * vect[i]
            return(result)

    def Length(self):
        return math.sqrt(self.Dot(self))

    def Equals(self, v2):
        if self.Subtract(v2).Length() <= .001:
            return True
        else:
            return False

    def Normalize(self):
        if self.Length() != 0:
            return self.Scalar_Mult(1 / self.Length())
        else:
            return self

    def Distance(self, v2):
        temp_vector = self.Subtract(v2)
        return temp_vector.Length()

    def __str__(self):
        vals = ''
        for val in self.pts:
            vals = vals + str(val) + ','

        return '(' + vals[:-1] + ')'

    def __repr__(self):
        return self.__str__()

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

    def From_Spherical(self, r, phi, theta):
        self.x = r * math.sin(theta) * math.cos(phi)
        self.y = r * math.sin(theta) * math.sin(phi)
        self.z = r * math.cos(theta)

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
