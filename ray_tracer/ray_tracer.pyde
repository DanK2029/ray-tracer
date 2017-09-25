from scene_components import *

sphere_list = []
vertex_list = []
triangle_list = []
light_list = []
FOV = 0
b1, b2, b3 = (0,)*3
cx, cy, cz = (0,)*3

def setup():
    size(500, 500) 
    noStroke()
    colorMode(RGB, 1.0)  # Processing color values will be in [0, 1]  (not 255)
    background(0, 0, 0)

# read and interpret the appropriate scene description .cli file based on key press
def keyPressed():
    if key == '1':
        interpreter("i1.cli")
    elif key == '2':
        interpreter("i2.cli")
    elif key == '3':
        interpreter("i3.cli")
    elif key == '4':
        interpreter("i4.cli")
    elif key == '5':
        interpreter("i5.cli")
    elif key == '6':
        interpreter("i6.cli")
    elif key == '7':
        interpreter("i7.cli")
    elif key == '8':
        interpreter("i8.cli")
    elif key == '9':
        interpreter("i9.cli")
    elif key == '0':
        interpreter("i10.cli")
    elif key == 't':
        interpreter("iT.cli")

def interpreter(fname):
    global b1, b2, b3, FOV
    b1, b2, b3 = (0,)*3
    surface = None
    fname = "data/" + fname
    # read in the lines of a file
    with open(fname) as f:
        lines = f.readlines()

    # parse each line in the file in turn
    for line in lines:
        words = line.split()  # split the line into individual tokens
        if len(words) == 0:   # skip empty lines
            continue
        if words[0] == 'sphere':
            radius = float(words[1])
            x = float(words[2])
            y = float(words[3])
            z = float(words[4])
            create_sphere(radius,x,y,z, surface)
        elif words[0] == 'fov':
            FOV = radians(float(words[1]))
        elif words[0] == 'background':
            b1 = float(words[1])
            b2 = float(words[2])
            b3 = float(words[3])
        elif words[0] == 'light':
            l1 = float(words[1])
            l2 = float(words[2])
            l3 = float(words[3])
            l4 = float(words[4])
            l5 = float(words[5])
            l6 = float(words[6])
            add_light(l1, l2, l3, l4, l5, l6)
        elif words[0] == 'surface':
            Cdr = float(words[1])
            Cdg = float(words[2])
            Cdb = float(words[3])
            Car = float(words[4])
            Cag = float(words[5])
            Cab = float(words[6])
            Csr = float(words[7])
            Csg = float(words[8])
            Csb = float(words[9])
            P   = float(words[10])
            Krefl = float(words[11])
            surface = Surface(Cdr, Cdg, Cdb, Car, Cag, Cab, Csr, Csg, Csb, P, Krefl)
        elif words[0] == 'begin':
            beginShape()
        elif words[0] == 'vertex':
            vx = float(words[1])
            vy = float(words[2])
            vz = float(words[3])
            add_vertex(vx, vy, vz, surface)
        elif words[0] == 'end':
            endShape()
        elif words[0] == 'write':
            render_scene()    # render the scene
            save(words[1])  # write the image to a file
            pass

# render the ray tracing scene
def render_scene():
    println("Rendering...")
    #println(light_list)
    #println(sphere_list)
    global FOV, b1, b2, b3, vertex_list
    set_triangle_list()
    
    for j in range(height):
        for i in range(width):
            # create an eye ray for pixel (i,j) and cast it into the scene
            z = -1.0
            k = tan(FOV/2)
            xp = i / abs(z)
            yp = (height - j) / abs(z)
            xpp = (xp - (width/2) )  * ( (2*k) / width)
            ypp = (yp - (height/2) ) * ( (2*k) / height)
            zpp = -1.0
            dx = xpp - cx 
            dy = ypp - cy
            dz = zpp - cz
            d = [dx, dy, dz]
            dlen = vec_len(d)
            d = [dx/dlen, dy/dlen, dz/dlen]
            ray = Ray(cx,cy,cz, d[0], d[1], d[2])
            
            L = cast_ray(ray, i, j, 100)
            Cr = L[0]
            Cg = L[1]
            Cb = L[2]
            pix_color = color(Cr, Cg, Cb)  # you should calculate the correct pixel color here
            set (i, j, pix_color)         # fill the pixel with the calculated color
    reset()
    println("Done")
    pass
    
# resets the all global lists for other scenes
def reset():
    global light_list, sphere_list, vertex_list, b1, b2, b3, triangle_list
    sphere_list = []
    light_list = []
    vertex_list = []
    triangle_list = []
    b1, b2, b3 = (0,)*3
    
def set_triangle_list():
    global triangle_list, vertex_list
    list_length = len(vertex_list)
    if (list_length%3 != 0):
        println("vertex list does not have 3 vertices per triangle")
    for i in range(0, list_length, 3):
        tri = My_Triangle(vertex_list[i], vertex_list[i+1], vertex_list[i+2], vertex_list[i][3])
        triangle_list.append(tri)
    
def add_vertex(x, y, z, surface):
    v = [x, y, z, surface]
    vertex_list.append(v)
    
def create_sphere(radius, x, y, z, surface):
    my_sphere = My_Sphere(radius, x, y, z, surface)
    sphere_list.append(my_sphere)
    
def add_light(x, y, z, r, g, b):
    light = Light(x, y, z, r, g, b)
    light_list.append(light)
    
def cast_ray(ray, i, j, maxDepth):
    global  b1, b2, b3
    L = [b1, b2, b3]
    if(maxDepth == 0):
        return L
    ray_intersect_info = ray_intersect_scene(ray, i, j, 0)
    
    T = ray_intersect_info[0]
    type = ray_intersect_info[1]
    objID = ray_intersect_info[2]
    
    if (T != float("inf") or type != -1):
        if(type == 0):
            hx = ray.getX(T)
            hy = ray.getY(T)
            hz = ray.getZ(T)
            Ka = sphere_list[objID].surface.get_ambient_color_vec()
            Kd = sphere_list[objID].surface.get_diffuse_color_vec()
            Ks = sphere_list[objID].surface.get_specular_color_vec()
            P  = sphere_list[objID].surface.get_P()
            N  = sphere_list[objID].get_N(hx, hy, hz)
            Nlen = vec_len(N)
            N = [N[0]/Nlen, N[1]/Nlen, N[2]/Nlen]
            multiple_light_sum = [0,0,0]
            for k in range(len(light_list)):
                lPos = light_list[k].get_Pos_Vec()
                l = [lPos[0]-hx, lPos[1]-hy, lPos[2]-hz]
                lLen = vec_len(l)
                l = [l[0]/lLen, l[1]/lLen, l[2]/lLen]
                NdotL = dot_product(N, l)
                I = light_list[k].get_Color_Vec()
                V = vec_sub(ray.get_origin(), [hx, hy, hz])
                Vlen = vec_len(V)
                V = [V[0]/Vlen, V[1]/Vlen, V[2]/Vlen]
                H = vec_add(V,l)
                Hlen = vec_len(H)
                H = [H[0]/Hlen, H[1]/Hlen, H[2]/Hlen]
                specular_value = vec_mult(Ks, vec_const_mult(I, pow(max(0.0, dot_product(N, H)), P)))
                diffuse_value = vec_mult(Kd, vec_const_mult(I, max(0, NdotL)))
                diff_and_spec_value = vec_add(diffuse_value, specular_value)
                #if(i==410 and j==330):
                    #println("Normal Vec: "+str(N))
                #e = 0.00001
                #hit_vector = [hx + N[0]*e, hy + N[1]*e, hz + N[2]*e]
                shadow_val = cast_shadow_ray([hx,hy,hz], lPos, i, j)
                
                diff_and_spec_value = vec_const_mult(diff_and_spec_value, shadow_val)
                multiple_light_sum = vec_add(multiple_light_sum, diff_and_spec_value)
            
            Krefl = sphere_list[objID].surface.Krefl
            if(Krefl != 0):
                d = ray.get_dir()
                refld = vec_sub(d, vec_const_mult(N , 2*dot_product(d, N)))
                refld_len = vec_len(refld)
                refld = [refld[0]/refld_len, refld[1]/refld_len, refld[2]/refld_len]
                e = vec_const_mult(refld, 0.0001)
                refl_ray = Ray(hx+e[0], hy+e[1], hz+e[2], refld[0], refld[1], refld[2])
                Rc = cast_ray(refl_ray, i, j, maxDepth-1)
                #if(i==410 and j==330):
                    #println([ray.get_dir(), refl_ray.get_dir()])
                Rc = vec_const_mult(Rc, Krefl)
                multiple_light_sum = vec_add(multiple_light_sum, Rc)
                
            L = vec_add(Ka, multiple_light_sum)
            return L
                
        if(type == 1):#triangle hit
            hx = ray.getX(T)
            hy = ray.getY(T)
            hz = ray.getZ(T)
            Ka = triangle_list[objID].surface.get_ambient_color_vec()
            Kd = triangle_list[objID].surface.get_diffuse_color_vec()
            Ks = triangle_list[objID].surface.get_specular_color_vec()
            P  = triangle_list[objID].surface.get_P()
            N = triangle_list[objID].N
            Nlen = vec_len(N)
            N = [-N[0]/Nlen, -N[1]/Nlen, -N[2]/Nlen]
            multiple_light_sum = [0,0,0]
            for k in range(len(light_list)):
                lPos = light_list[k].get_Pos_Vec()
                l = [lPos[0]-hx, lPos[1]-hy, lPos[2]-hz]
                lLen = vec_len(l)
                l = [l[0]/lLen, l[1]/lLen, l[2]/lLen]
                NdotL = dot_product(N, l)
                I = light_list[k].get_Color_Vec()
                V = vec_sub(ray.get_origin(), [hx, hy, hz])
                Vlen = vec_len(V)
                V = [V[0]/Vlen, V[1]/Vlen, V[2]/Vlen]
                H = vec_add(V,l)
                Hlen = vec_len(H)
                H = [H[0]/Hlen, H[1]/Hlen, H[2]/Hlen]
                NdotH = dot_product(N, H)
                specular_value = vec_mult(Ks, vec_const_mult(I, pow(max(0.0, NdotH), P)))
                diffuse_value = vec_mult(Kd, vec_const_mult(I, max(0, NdotL)))
                diff_and_spec_value = vec_add(diffuse_value, specular_value)
                shadow_val = cast_shadow_ray([hx, hy, hz], lPos, i, j)
                diff_and_spec_value = vec_const_mult(diff_and_spec_value, shadow_val)
                multiple_light_sum = vec_add(multiple_light_sum, diff_and_spec_value)
                
            Krefl = triangle_list[objID].surface.Krefl
            if(Krefl != 0):
                d = ray.get_dir()
                refld = vec_sub(d, vec_const_mult(N , 2*dot_product(d, N)))
                refld_len = vec_len(refld)
                refld = [refld[0]/refld_len, refld[1]/refld_len, refld[2]/refld_len]
                e = vec_const_mult(refld, 0.00001)
                refl_ray = Ray(hx+e[0], hy+e[1], hz+e[2], refld[0], refld[1], refld[2])
                Rc = cast_ray(refl_ray, i, j, maxDepth-1)
                Rc = vec_const_mult(Rc, Krefl)
                multiple_light_sum = vec_add(multiple_light_sum, Rc) 

            L = vec_add(Ka, multiple_light_sum)
            return L
    return L
    

def ray_intersect_scene(ray, i, j, sh):
    T = float("inf")
    dx = ray.dx
    dy = ray.dy
    dz = ray.dz
    
    x0 = ray.x0
    y0 = ray.y0
    z0 = ray.z0
    
    t = float("inf")
    
    type = -1 #sphere == 0  triangle == 1
    objID = -1
    for k in range(len(sphere_list)):
        r  = sphere_list[k].radius
        cx = sphere_list[k].cx
        cy = sphere_list[k].cy
        cz = sphere_list[k].cz
        
        a = sq(dx) + sq(dy) + sq(dz)
        b = 2 *((dx*(x0-cx))+(dy*(y0-cy))+(dz*(z0-cz)))
        c = sq(x0-cx) + sq(y0-cy) + sq(z0-cz) - sq(r)
        
        disc = sq(b) - 4*a*c
        
        if(disc > 0):
            t1 = (-b + sqrt(sq(b)-(4*a*c)))/(2*a)
            t2 = (-b - sqrt(sq(b)-(4*a*c)))/(2*a)
            
            if(t1 >= 0 and t2 >= 0):
                t = min(t1, t2)
                
            elif(t1 >= 0 and t2 < 0):
                t = t1
                
            elif(t1 < 0 and t2 >= 0):
                t = t2
                
        elif(disc == 0):
            t = -b / (2*a)
            type = 0
            
        elif(disc < 0):
            t = float("inf")
        
        if(t < T and t > 0.0001):
            T = t
            objID = k
            type = 0
    
            
    for k in range(len(triangle_list)):
        tri_t = line_intersect_triangle(ray, triangle_list[k], i, j, sh, k)
        if(tri_t != float("inf") and tri_t <= T and tri_t > 0):
            T = tri_t
            type = 1
            objID = k

    return([T, type, objID])
    

def cast_shadow_ray(ray_origin, light_origin, i, j):
    #find t for ray intesect light
    shadow_vec_dir = [light_origin[0]-ray_origin[0], light_origin[1]-ray_origin[1], light_origin[2]-ray_origin[2]]
    shadow_vec_length = vec_len(shadow_vec_dir)
    shadow_vec_dir = [shadow_vec_dir[0]/shadow_vec_length, shadow_vec_dir[1]/shadow_vec_length, shadow_vec_dir[2]/shadow_vec_length]
    e = vec_const_mult(shadow_vec_dir, 0.0001)
    ray_origin = vec_add(ray_origin, e)
    shadow_ray = Ray(ray_origin[0], ray_origin[1], ray_origin[2],  shadow_vec_dir[0], shadow_vec_dir[1], shadow_vec_dir[2])
    hit_object = ray_intersect_scene(shadow_ray, i, j, 1)
    
    
    if(hit_object[0] == float("inf") or hit_object[0] >= shadow_vec_length):
        return 1

    else:
        return 0
    

def line_intersect_triangle(ray, tri, i, j, sh, k):
    #does the ray hit the plane
    R = [ray.dx, ray.dy, ray.dz]
    Rlen = vec_len(R)
    R = [R[0]/Rlen, R[1]/Rlen, R[2]/Rlen]
    A = tri.v0
    B = tri.v1
    C = tri.v2
    P = [ray.x0, ray.y0, ray.z0]
    N = tri.N
    NdotR = dot_product(N, R)
    if(NdotR == 0):
        return float("inf")
    #is the hit point inside the triangle
    t_longer = (dot_product(tri.N, A) - dot_product(P, tri.N))/(NdotR)
    t = (dot_product(vec_sub(A, P), N))/(NdotR)
    H = ray.getPoint(t)
    
    if(dot_product(cross_product(vec_sub(B, A), vec_sub(H, A)), N) < -0.000001):
        return float("inf")
    if(dot_product(cross_product(vec_sub(C, B), vec_sub(H, B)), N) < -0.000001):
        return float("inf")
    if(dot_product(cross_product(vec_sub(A, C), vec_sub(H, C)), N) < -0.000001):
        return float("inf")

    return t
        
def dot_product(v1, v2):
    return (v1[0] * v2[0]) + (v1[1] * v2[1]) + (v1[2] * v2[2])

def vec_mult(v1, v2):
    return [(v1[0] * v2[0]), (v1[1] * v2[1]), (v1[2] * v2[2])]

def vec_const_mult(v, c):
    return [v[0]*c, v[1]*c, v[2]*c]

def vec_add(v1, v2):
    return [(v1[0] + v2[0]), (v1[1] + v2[1]), (v1[2] + v2[2])]

def vec_sub(v1, v2):
    return [(v1[0] - v2[0]), (v1[1] - v2[1]), (v1[2] - v2[2])]

def vec_div(v1, v2):
    if(v2[0] == 0 or v2[1] == 0 or v2[2] == 0):
        println("vec_div: Divide by Zero")
    return [(v1[0]/v2[0]), (v1[1]/v2[1]), (v1[2]/v2[2])]

    
# should remain empty
def draw():
    pass