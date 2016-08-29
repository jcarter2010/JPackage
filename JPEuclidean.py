import math
from JPMath import Vector2

class Constants:

    '''
    There is nothing here
    '''

class Functions:

    def gcd(self, a, b):
        while b:
            a, b = b, a%b
        return a

    def roots(self, a, b, c):
        d = b**2-4*a*c
        if d < 0:
            return None
        x1 = (-b+math.sqrt((b**2)-(4*(a*c))))/(2*a)
        x2 = (-b-math.sqrt((b**2)-(4*(a*c))))/(2*a)
        return (x1, x2)

class Circle:

    global C
    global r

    def __init__(self, C, r):
        self.C = C
        self.r = r

    def Get_Intersections_With_Circle(self, circ2):
        seg = Line_Segment(self.C, circ2.C)
        perp = seg.Create_Perpendicular().To_Line()
        print(perp)
        intersections = perp.Get_Intersections_With_Circle(circ2)
        if str(intersections) == '((nan,nan), (nan,nan))':
            x = seg.Get_Midpoint().x
            ys = self.Get_Y_Values(x)
            return (Vector2(x, ys[0]), Vector2(x, ys[1]))
        return intersections

    def Get_Y_Values(self, x):
        return (math.sqrt(self.r ** 2 - x ** 2), -math.sqrt(self.r ** 2 - x ** 2))

    def __str__(self):
        return 'C : {}, R : {}'.format(str(self.C), str(self.r))

    def __repr__(self):
        return self.__str__()

class Angle:

    global A
    global B
    global C
    global theta
    global dirBA
    global dirBC

    def __init__(self, A, B, C):
        self.A = A
        self.B = B
        self.C = C
        self.dirBA = A.Subtract(B)
        self.dirBC = C.Subtract(B)
        self.theta = self.dirBA.Angle_Between_Vectors(self.dirBC)

    def Get_Bisected(self):
        theta = self.theta / 2
        bisector = self.dirBC
        bisector.Rotate_Around_Point(self.B, theta)
        ang1 = Angle(self.A, self.B, bisector)
        ang2 = Angle(bisector, self.B, self.C)
        return (ang1, ang2)

    def Get_Bisector(self):
        theta = self.theta / 2
        bisector = self.dirBC
        bisector.Rotate_Around_Point(self.B, theta)
        return bisector

    def __str__(self):
        return 'A : {}, B : {}, C : {}, THETA : {}'.format(str(self.A), str(self.B), str(self.C), self.theta)

    def __repr__(self):
        return self.__str__()

class Line:

    global m
    global b

    def __init__(self, eq):
        #eq in form of 'mx+b'
            #if b < 0, then write(mx+-b)
            #if b = 0, then write (mx+0)
        eq = eq.replace('e+', 'ep')
        self.m = float(eq[:eq.find('x')].replace('ep', 'e+'))
        self.b = float(eq[eq.find('+')+1:].replace('ep', 'e+'))

    def Check_Congruent(self, line2):
        if abs(self.m - line2.m) <= 0.001:
            if abs(self.b - line2.b) <= 0.001:
                return True
        return False

    def Check_Parallel(self, line2):
        if abs(self.m - line2.m) <= 0.001 and self.Check_Congruent(line2) == False:
            return True
        return False

    def Get_Intersection_With_Line(self, line2):
        if Check_Parallel(line2) == False:
            x = (line2.b - self.b) / (self.m - line2.m)
            y = self.m * x + self.b
            return Vector2(x, y)
        return None

    def Get_Y_Value(self, x):
        return self.m * x + self.b

    def Get_Intersections_With_Circle(self, circ):
        r = circ.r
        a = self.m ** 2 + 1
        b = 2 * self.m * self.b
        c = self.b ** 2 - r ** 2
        funcs = Functions()
        roots = funcs.roots(a, b, c)
        r_list = []
        if roots[0] == None:
            r_list.append(None)
        else:
            r_list.append(Vector2(roots[0], self.Get_Y_Value(r)))
        if roots[0] == None:
            r_list.append(None)
        else:
            #if roots[1] != roots[0]:
            r_list.append(Vector2(roots[1], self.Get_Y_Value(r)))
            #else:
            #    r_list.append(Vector2(roots[1], -self.Get_Y_Value(r)))
        return tuple(r_list)

    def Create_Perpendicular(self, point):
        if self.m != 0:
            m_prime = math.pow(self.m, -1) * -1
            b_prime = point.y - m_prime * point.x
            return Line('{}x+{}'.format(m_prime, b_prime))
        else:
            m_prime = float('inf')
            b_prime = float('NaN')
            return Line('{}x+{}'.format(m_prime, b_prime))

    def __str__(self):
        return 'y = {} x + {}'.format(self.m, self.b)

    def __repr__(self):
        return self.__str__()

class Triangle:

    global A
    global B
    global C
    global ABC
    global BCA
    global CAB
    global AB
    global BC
    global CA

    def __init__(self, A, B, C):
        self.A = A
        self.B = B
        self.C = C
        self.ABC = Angle(A, B, C)
        self.BCA = Angle(B, C, A)
        self.CAB = Angle(C, A, B)
        self.AB = Line_Segment(A, B)
        self.BC = Line_Segment(B, C)
        self.CA = Line_Segment(C, A)

    def Check_Congruent(self, tri2):
        if self.A.Equals(tri2.A) and self.B.Equals(tri2.B) and self.C.Equals(tri2.C):
            return True
        return False

    def Find_Center(self):
        return Vector2((self.A.x + self.B.x + self.C.x) / 3, (self.A.y + self.B.y + self.C.y) / 3)

    def __str__(self):
        return 'A : {}, B : {}, C : {}'.format(str(self.A), str(self.B), str(self.C))

    def __repr__(self):
        return self.__str__()

class Line_Segment:

    global A
    global B

    def __init__(self, A, B):
        self.A = A
        self.B = B

    def Check_Congruent(self, seg2):
        if self.A.Equals(seg2.A) and self.B.Equals(seg2.B):
            return True
        return False

    def Get_Midpoint(self):
        return Vector2((self.A.x + self.B.x) / 2, (self.A.y + self.B.y) / 2)

    def Get_Intersection_With_Segment(self, seg2):
        l1 = self.To_Line()
        l2 = seg2.To_Line()
        intersect = l1.Get_Intersection_With_Line(l2)
        if intersect != None:
            if intersect.x >= self.A.x and intersect.x <= self.B.x:
                return intersect
        return None

    def Get_Intersection_With_Line(self, line2):
        l1 = self.To_Line()
        intersect = l1.Get_Intersection_With_Line(line2)
        if intersect != None:
            if intersect.x >= self.A.x and intersect.x <= self.B.x:
                return intersect
        return None

    def Get_Intersections_With_Circle(self, circ):
        l1 = self.To_Line()
        intersections = l1.Get_Intersections_With_Circle(circ)
        if intersections[0].x <= self.B.x and intersections[0].x >= self.A.x:
            if intersections[1].x <= self.B.x and intersections[1].x >= self.A.x:
                return (intersections[0], intersections[1])
            return (intersections[0],)
        if intersections[1].x <= self.B.x and intersections[1].x >= self.A.x:
            return (intersections[1],)
        return None

    def Create_Perpendicular(self):
        mid = self.Get_Midpoint()
        dist = self.A.Distance(self.B) / 2
        dir1 = self.A.Subtract(self.B).Normalize()
        dir2 = self.B.Subtract(self.A).Normalize()
        dir1 = Vector2(dir1.y, -dir1.x)
        dir2 = Vector2(dir2.y, -dir2.x)
        A_prime = mid.Add(dir1.Scalar_Mult(dist))
        B_prime = mid.Add(dir2.Scalar_Mult(dist))
        return Line_Segment(A_prime, B_prime)


    def To_Line(self):
        v = Vector2(self.A.x - self.B.x, self.A.y - self.B.y)
        m = v.Get_Slope()
        b = self.B.y - m * self.B.x
        return Line('{}x+{}'.format(m, b))

    def __str__(self):
        return 'A : {}, B : {}'.format(str(self.A), str(self.B))

    def __repr__(self):
        return self.__str__()

class Ray:

    global start
    global direction

    def __init__(self, start, direction): #start and point are both JPMath.Vector2's
        self.start = start
        self.direction = direction.Normalize()

    def Check_Congruent(self, ray2):
        if self.start.Equals(ray2.start) and self.dir.Equals(ray2.dir):
            return True
        return False

    def Get_Intersections_With_Circle(self, circ):
        l = self.To_Line()
        intersections = self.Get_Intersections_With_Circle(circ)
        i_list = []
        if self.direction.x > 0:
            for i in intersections:
                if i > self.start.x:
                    i_list.append(Vector2(i, l.Get_Y_Value(i)))
        if self.direction.x == 0:
            for i in intersections:
                if i == self.start.x:
                    i_list.append(Vector2(i, l.Get_Y_Value(i)))
        if self.direction.x < 0:
            for i in intersections:
                if i < self.start.x:
                    i_list.append(Vector2(i, l.Get_Y_Value(i)))
        return tuple(i_list)

    def Get_Intersection_With_Line(self, line):
        l = self.To_Line()
        intersection = self.Get_Intersection_With_Line(line)
        if self.direction.x > 0:
            if intersection.x > self.start.x:
                return(intersection)
        if self.direction.x == 0:
            if intersection.x == self.start.x:
                return(intersection)
        if self.direction.x < 0:
            if intersection.x < self.start.x:
                return(intersection)
        return None

    def To_Line(self):
        seg = Line_Segment(self.start, self.start.Add(self.direction))
        return seg.To_Line()

    def __str__(self):
        return 'P : {} , DIR : '.format(str(self.start), str(self.direction))

    def __repr__(self):
        return self.__str__()

class Equations:

    '''
    There is nothing here
    '''
