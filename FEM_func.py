##### HOLDS ALL FEM FUNCTIONS USED IN PROJECT.IPYNB #####

import numpy as np


#### MESH TOOL ####
def mesh(bot,top,left,right):

    # returns [XYZ, CON, DOF] 
    # XYZ - array of nodal coordinates [number of elements x 2]
    # CON - array of node numbers for elements (linear QUADS [number of elements x 4])
    # DOF - array of element DOFs (4-node (linear) quadrilateral element => [number of elements x 8]); 
    
    # number of nodes and elements in the domain
    nnodesx = len(bot)                     # number of horizontal nodes
    nnodesy = len(left)                    # number of vertical nodes 
    nelx = nnodesx-1                       # number of horizontal elements
    nely = nnodesy-1                       # number of vertical elements
    nnodes = nnodesx*nnodesy               # total number of nodes    

    # dimensions of the domain
    lx = bot[nnodesx-1] - bot[0]           # length of the domain in x-direction (horizontal)
    ly = left[nnodesy-1] - left[0]         # length of the domain in y-direction (vertical)

    # GENERATE COORDINATES OF NODES 'XYZ'
    XYZ = np.zeros((nnodes,2))            # two-column array [nnodes x 2] containing all nodal coordinates  
    for i in range(nnodesy):              # loop over all nodes on the vertical sides 
        yl = left[i] - left[0]            # distance between node 'i' and left-bottom node '1'
        dy = right[i] - left[i]           # distance between the corresponing nodes j on top and bottom 
        for j in range(nnodesx):          # loop over all nodes on the horizontal sides
            xb = bot[j] - bot[0]          # distance between node 'j' and bottom-left node '1' 
            dx = top[j] - bot[j]          # distance between nodes 'j' on opposite sides (top and bottom)

            x = (dx*yl+xb*ly)/(ly-dx*dy/lx) # x-coordinate (horizontal) of a node in the interior of the domain
            y = dy/lx*x+yl                  # y-coordinate (vertical) of a node in the interior of the domain

            XYZ[j+i*nnodesx, 0] = x + bot[0]  # coordinate 'x' in the global coordinate system 
            XYZ[j+i*nnodesx, 1] = y + left[0] # coordinate 'y' in the global coordinate system

    # NODE NUMBERS FOR ELEMENTS 
    nel = nelx*nely                              # total number of elements in the domain
    CON = np.zeros((nel,4), dtype=int)           # [nel*4] array of node number for each element
    for i in range(nely):                        # loop over elements in the vertical direction 
        for j in range(nelx):                    # loop over elements in the horizontal direction 
            # element 'el' and corresponding node numbers
            CON[j+i*nelx, :] = [j+i*nnodesx, j+i*nnodesx+1,j+(i+1)*nnodesx+1, j+(i+1)*nnodesx] 

    # Global DOF for each element (4-node (linear) quadrilateral element)
    DOF = np.zeros((nel,2*4), dtype=int)
    for i in range(nel):
        # defines single row of DOF for each element 'i'
        DOF[i,:] = [CON[i,0]*2, CON[i,1]*2-1, CON[i,1]*2, CON[i,1]*2+1,CON[i,2]*2, CON[i,2]*2+1, CON[i,3]*2, CON[i,3]*2+1]
        
    return XYZ, CON, DOF, nelx, nely


### CREATE NODAL POSISTIONS/COORDINATES
def create_grid(delta_x, delta_y):
    ### GRID DIMESNSIONS
    length = 8.0
    height = 4.0
    ### NODAL STEP X
    step = length/delta_x
    bot = np.arange(0.0,length+step,step); top = bot
    ### NODAL STEP Y
    step = height/delta_y
    left = np.arange(0.0,height+step,step); right = left
    
    return bot, top, right, left


### ELASTICITY MATRIX

def plane_strain(E, nu):
    
    const = E/((1+nu)*(1-2*nu))
    C = const*np.array([[1.0-nu, nu, 0.0],
                        [nu, 1.0-nu, 0.0], 
                        [0.0, 0.0, 0.5*(1-2*nu)]])

    return C

### JACOBIAN DETERMINANT AND INVERSE TO OUTPUT SHAPE FUNCTION DERIVATIVE
def Jacobian_dN_dx(pos,gp):
    # Jacobian
    xi = gp[0]
    eta = gp[1]
    ### SHAPE FUNCTION DERIVATIVE NATURAL COORDS
    dN_dx_nat = np.array([[eta-1,1-eta,1+eta,-eta-1],[xi-1,-xi-1,1+xi,1-xi]])*0.25
    J = np.dot(dN_dx_nat,pos)
    
    ### J DETERMINANT
    detJ = np.linalg.det(J)
    
    ### SHAPE FUNCTION DERIVATIVE REAL COORDS
    invJ = np.linalg.inv(J)
    dN_dx_real = np.dot(invJ,dN_dx_nat)
    
    return detJ, dN_dx_real


### STARIN-DISPLACEMENT MATRIX - B
def strain_displacement(dN_dx):
    B1 = np.array([[dN_dx[0,0],0],[0,dN_dx[1,0]],[dN_dx[1,0],dN_dx[0,0]]])
    B2 = np.array([[dN_dx[0,1],0],[0,dN_dx[1,1]],[dN_dx[1,1],dN_dx[0,1]]])
    B3 = np.array([[dN_dx[0,2],0],[0,dN_dx[1,2]],[dN_dx[1,2],dN_dx[0,2]]])
    B4 = np.array([[dN_dx[0,3],0],[0,dN_dx[1,3]],[dN_dx[1,3],dN_dx[0,3]]])
    B = np.concatenate((B1,B2,B3,B4),axis=1)
    Bt = np.transpose(B)
    
    return B, Bt
