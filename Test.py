from JPMath import Vector2
import math
from JPEuclidean import Circle
from JPEuclidean import Line_Segment
from JPEuclidean import Line
from JPEuclidean import Angle
from JPEuclidean import Ray

r = Ray(Vector2(0, 0), Vector2(1, 1))
#print(r.To_Line())

circ1 = Circle(Vector2(0, 0), 1)
circ2 = Circle(Vector2(1, 0), 1)

#print(circ1)
#print(circ2)

print(circ1.Get_Intersections_With_Circle(circ2))
#segm = Line_Segment(Vector2(-1, -1), Vector2(1, 1))
#line = Line('1x+1')

#print(segm.Get_Intersections_With_Circle(circ))
#print(segm.Create_Perpendicular())
#print(line.Create_Perpendicular(Vector2(1, 0)))
#ang = Angle(Vector2(0, 1), Vector2(0, 0), Vector2(1, 0))
#print(ang.Get_Bisected())
