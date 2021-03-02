"""utils

    This module should contain small usefull functions.

"""
import math

def _cartesian_distance(x1:float,x2:float,y1:float,y2:float,z1:float,z2:float)->float:

    return math.sqrt((x1 - x2)**2+(y1 - y2)**2+(z1 - z2)**2)