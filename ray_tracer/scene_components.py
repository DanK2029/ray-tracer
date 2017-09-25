class Ray():
    def __init__(self, x0, y0, z0, dx, dy, dz):
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.dx = dx
        self.dy = dy
        self.dz = dz
    
    def get_origin(self):
        return [self.x0, self.y0, self.z0]
    
    def get_dir(self):
        return [self.dx, self.dy, self.dz]
    
    def getPoint(self, t):
        return [self.x0 + self.dx*t, self.y0 + self.dy*t, self.z0 + self.dz*t]
    
    def getX(self, t):
        return self.x0 + self.dx*t
    
    def getY(self, t):
        return self.y0 + self.dy*t
    
    def getZ(self, t):
        return self.z0 + self.dz*t

class My_Sphere():
    def __init__(self, radius, cx, cy, cz, surface):
        self.radius = radius
        self.cx = cx
        self.cy = cy
        self.cz = cz
        self.surface = surface
        
        
    def get_N(self, px, py, pz):
        vx = px - self.cx
        vy = py - self.cy
        vz = pz - self.cz
        return [vx, vy, vz]
    
class My_Triangle():
    def __init__(self, v0, v1, v2, surface):
        self.v0 = v0
        self.v1 = v1
        self.v2 = v2
        self.surface = surface
        self.N = get_triangle_N(v0, v1, v2)

def get_triangle_N(v0, v1, v2):
    r1 = [v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2]]
    r2 = [v2[0]-v0[0], v2[1]-v0[1], v2[2]-v0[2]]
    cp = cross_product(r1, r2)
    l = vec_len(cp)
    return [cp[0]/l, cp[1]/l, cp[2]/l]

def vec_len(v):
    return sqrt(sq(v[0])+sq(v[1])+sq(v[2]))
    
def cross_product(A, B):
    d0 = A[1]*B[2]-A[2]*B[1]
    d1 = A[2]*B[0]-A[0]*B[2]
    d2 = A[0]*B[1]-A[1]*B[0]
    return [d0, d1, d2]
        
class hit():
    def __init__(self, t, object, L):
        self.t = t
        self.object = object
        self.L = L
        
    def getL(self):
        return self.L

class Surface():
    def __init__(self, Cdr, Cdg, Cdb, Car, Cag, Cab, Csr, Csg, Csb, P, Krefl):
         self.Cdr   = Cdr
         self.Cdg   = Cdg
         self.Cdb   = Cdb
         self.Car   = Car
         self.Cag   = Cag
         self.Cab   = Cab
         self.Csr   = Csr
         self.Csg   = Csg
         self.Csb   = Csb
         self.P     = P
         self.Krefl = Krefl
         
    def get_P(self):
        return self.P
         
    def get_diffuse_color_vec(self):
        return [self.Cdr, self.Cdg, self.Cdb]
    
    def get_ambient_color_vec(self):
        return [self.Car, self.Cag, self.Cab]
    
    def get_specular_color_vec(self):
        return [self.Csr, self.Csg, self.Csb]
         
class Light():
    def __init__(self, x, y, z, r, g, b):
        self.x = x
        self.y = y
        self.z = z
        self.r = r
        self.g = g
        self.b = b
        
    def get_Pos_Vec(self):
        return [self.x, self.y, self.z]
    
    def get_Color_Vec(self):
        return [self.r, self.g, self.b]
        