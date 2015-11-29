add_library('peasycam')
add_library('physics')

from itertools import combinations
import hyperbolic_tiling as ht

def setup():
    size(640, 640, P3D)
    
    global halfWidth, halfHeight
    halfWidth = width / 2
    halfHeight = height / 2
    PeasyCam(this, halfWidth, halfHeight, 0.0, 500)
    
    #hyperbolic tiling
    p, q, layers = 6, 4, 4 
    verts_out, edges_out, faces_out = ht.poincare_tiling(p, q, layers) 
    global v2, e2, f2
    v2, e2, f2 = ht.cone_tris(p, verts_out, edges_out, faces_out)
    
    print(len(f2))
    
    #physics spring system  
    # system parameters
    mass = 1.1 # mass of particles
    restlength = 1.0  # Rest length of springs between nodes
    strength = 0.99  # Strength of springs between nodes
    damping = 0.0   # Damping of springs between nodes
    drag = 0.5  # Physics System drag (friction) 
    repulsion = -0.007  # Repulsion force between nodes
    repulsion_mindist = 0.1  # Distance from node where repulsion force begins to decrease
    
    global physics
    physics = ParticleSystem( 0, drag )
    #particles at verticies
    for v in v2:
        physics.makeParticle(mass, *v)
    #springs along edges
    for e in e2:
        p1 = physics.getParticle(e[0])
        p2 = physics.getParticle(e[1])
        stg = strength
        physics.makeSpring(p1, p2, stg, damping, restlength) 
    # repulsion between all nodes 
    N = physics.numberOfParticles()
    ats = set(combinations(range(N),2))
    e22 = set(map(tuple, e2)) 
    ats2 = ats - e22
    e23 = {e[::-1] for e in e22}
    ats2 = ats2 - e23
    print('ats len', len(ats2), len(e2), len(ats))
    ninner = 0
    for i, j in ats:
         p1 = physics.getParticle(i)
         p2 = physics.getParticle(j)   
         rep = repulsion
         physics.makeAttraction(p1, p2, rep, repulsion_mindist) 
          
    print('no of edges', len(e2))
    print('no of attractions', physics.numberOfAttractions() ) 
    print('ninner', ninner)   
                  
def draw():    
    background(255)
    lights()
    smooth()
    stroke(0)
    fill(128)
    translate(halfWidth, halfHeight, -100)
    rotateX(PI/4)
    strokeWeight(0.01) # set thin because it gets scaled too
    scale(50) 
    
    physics.tick(1)
    
    beginShape(TRIANGLES)
    for face in f2:
        for v in face:
            p1 = physics.getParticle(v)            
            vertex( p1.position().x(), p1.position().y(), p1.position().z() ) 
    endShape()
    
    # color edges by lengths
    global lengths
    global n_long
    global n_short
    d_long = 1.2
    d_short = 0.8
    n_long = 0
    n_short = 0
    lengths = []
    for e in e2:
        p1 = physics.getParticle(e[0])
        p2 = physics.getParticle(e[1])
        d = ht.distance((p1.position().x(), p1.position().y(), p1.position().z() ), 
                        (p2.position().x(), p2.position().y(), p2.position().z() ))
        lengths.append(d)
        if d > d_long:
            stroke(0, 255, 0)
            n_long = n_long +1
        elif d < d_short:
            stroke(255, 0, 0)
            n_short = n_short + 1
        else: 
            stroke(0, 0, 0)
                
        line(p1.position().x(), p1.position().y(), p1.position().z(), p2.position().x(), p2.position().y(), p2.position().z()) 

def keyPressed():
    if key == 'w':
        #write obj file
        verts = []
        N = physics.numberOfParticles()
        for i in range(N):
            p1 = physics.getParticle(i)
            verts.append( (p1.position().x(), p1.position().y(), p1.position().z() ) )
        filename = 'ht_01.obj'
        write_obj(filename, verts, f2)
        print('obj file written: ' + filename)
    if key == 'e':
        #calculate mesh length data and print to console
        edge_length_data(e2)
        
def write_obj(filename, verts, faces):
        """
        minimal obj format
        verts is list of coordinate tuples
        faces lists of vert indicies in each face 
        """
        output = createWriter(filename)
        for v in verts:
                output.print( 'v ' + nf(v[0],0,5) + " " + nf(v[1], 0, 5) + " " + nf(v[2], 0, 5) + "\n")
        
        for face in faces:
                output.print('f ')             
                for vi in face:
                        # obj indices start with 1 not zero
                        output.print(nf(vi + 1) + " ")
                        
                output.print("\n")
        output.flush()
        output.close()    
        return

   
            
def edge_length_data(edges):
        
    mean = sum(lengths)/float(len(lengths))
    ss = sum((x-mean)**2 for x in lengths)
    stdev = (ss/float(len(lengths) - 1))**0.5
    
    print(mean, stdev, max(lengths), min(lengths), n_long, n_short)
    
        