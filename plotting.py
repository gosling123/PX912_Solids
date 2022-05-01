#### PLOTTING FUNCTIONS ####


#hello
import numpy as np
import matplotlib.pyplot as plt


### INITIAL GRID

def initial_grid(XYZ, CON):

    nnodes = len(XYZ)
    nels = len(CON)

    for i in range(nnodes):
        ### EXCLUDE NON-SOLID NODES
        if 2.0<XYZ[i,0]<6.0 and 0.0<=XYZ[i,1]<2.0:
            continue
        plt.plot(XYZ[i, 0], XYZ[i, 1], 'sk')


    for i in range(nels):
    
        ### EXCLUDE NON-SOLID ELEMENTS
        j = 3
        if 2.0<=XYZ[CON[i, j],0] <6.0 and 0<XYZ[CON[i, j],1]<=2:
            continue 
        plt.fill(XYZ[CON[i, :], 0], XYZ[CON[i, :], 1], edgecolor='grey', color = 'grey', fill=True)

    for i in range(4):                             #loop over all nodes within an element
        for j in range(nels):                  #loop over all elements
            sh=0.01
            plt.text(XYZ[CON[j,i],0]+sh,XYZ[CON[j,i],1]+sh, CON[j,i])


    plt.axhline(y = 2.0, color = 'red', linestyle = '--')

    plt.axvline(x = 2.0, color = 'blue', linestyle = '--')
    plt.axvline(x = 6.0, color = 'blue', linestyle = '--')


    plt.xlabel("$x_1$ (m)", fontsize=25)

    plt.ylabel("$x_2$ (m)", fontsize=25)

    plt.tick_params(axis='both', which='major', labelsize=20)


    plt.show()



### NEW GRID PLOT 

def new_grid_plot(XYZ, CON, d_reshape, scalar, mesh_density_lvl):
    
    scalar = 5e7 # Displacements tiny so scaled for visualisation


    newpos = XYZ + d_reshape*scalar

    nnodes = len(XYZ)
    nels = len(CON)



    # Original object
    for i in range(nnodes):
        if 2.0<XYZ[i,0]<6.0 and 0.0<=XYZ[i,1]<2.0:
            continue
        plt.plot(XYZ[i, 0], XYZ[i, 1], 'sk')
    
    for i in range(nels):
    
        j = 3
    
        if 2.0<=XYZ[CON[i, j],0] <6.0 and 0<XYZ[CON[i, j],1]<=2:
            continue 
        
        plt.fill(XYZ[CON[i, :], 0], XYZ[CON[i, :], 1], color = 'grey', edgecolor='black', fill=False)

    # New object
    for i in range(nnodes):
        if 2.0<XYZ[i,0]<6.0 and 0.0<=XYZ[i,1]<2.0:
            continue
        plt.plot(newpos[i, 0], newpos[i, 1], 'or')
        
    for i in range(nels):
    
        j = 3
    
        if 2.0<=XYZ[CON[i, j],0] <6.0 and 0<XYZ[CON[i, j],1]<=2:
            continue 
        
        plt.fill(newpos[CON[i, :], 0], newpos[CON[i, :], 1], edgecolor='r', fill=False)


    fname = "Deform_Pics/deform_" + str(mesh_density_lvl) + ".png"
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.xlabel(r'$x_{1}$', fontsize = 25)
    plt.ylabel(r'$x_{2}$', fontsize = 25)
    # plt.savefig(fname)
    plt.show()

### MIDPOINTS FOR STRESS/STRAIN CONTOUR

def midpoint(XYZ, CON):
    
    ### RETURNS MIDPOINTS OF ELEMENTS

    nels = len(CON)
    
    mid_points = []
    for i in range(nels):
        
        points = XYZ[CON[i, :]]
    
        delta_x = np.abs(points[0, 0] - points[1,0])/2
        delta_y = np.abs(points[0, 1] - points[len(CON[i,:])-1 ,1])/2
        
        x = points[0, 0] + delta_x
        y = points[0, 1] + delta_y
        
        mid = [x, y]
        
        mid_points.append(mid)
        
    mid_points = np.array(mid_points)

    return mid_points


def stress_strain_contour(x1, x2, zeps_data, zsig_data, tol, mesh_density_lvl, cmap = 'jet'):


    fig = plt.figure()


    alabel_size = 20
    title_size = 25
    ticks_size = 20
    cbar_size = 20

    N  = 200


    ax1 = plt.subplot(3, 2, 1)
    plt.tricontourf(x1, x2, zeps_data[0]*1e11, levels = np.linspace(tol*np.max(zeps_data[0])*1e11, np.max(zeps_data[0])*1e11, N) , cmap = cmap)
    ax1.set_title(r'$\epsilon_{11} (\times 10 ^{-11} m)$', fontsize = title_size)
    ax1.set_xlabel(r'$x_1$', fontsize = alabel_size)
    ax1.set_ylabel(r'$x_2$', fontsize = alabel_size)
    cbar = plt.colorbar(ax = plt.gca())
    plt.xticks(fontsize= ticks_size)
    plt.yticks(fontsize= ticks_size)
    cbar.ax.tick_params(labelsize=cbar_size)


    ax2 = plt.subplot(3, 2, 3)
    plt.tricontourf(x1, x2, zeps_data[1]*1e11, levels = np.linspace(tol*np.max(zeps_data[1])*1e11, np.max(zeps_data[1])*1e11, N) , cmap = cmap)
    ax2.set_title(r'$\epsilon_{22} (\times 10 ^{-11} m)$', fontsize = title_size)
    ax2.set_xlabel(r'$x_1$', fontsize = alabel_size)
    ax2.set_ylabel(r'$x_2$', fontsize = alabel_size)
    cbar = plt.colorbar(ax = plt.gca())
    plt.xticks(fontsize= ticks_size)
    plt.yticks(fontsize= ticks_size)
    cbar.ax.tick_params(labelsize=cbar_size)

    ax3 = plt.subplot(3, 2, 5)
    plt.tricontourf(x1, x2, zeps_data[2]*1e11, levels = np.linspace(tol*np.max(zeps_data[2])*1e11, np.max(zeps_data[2])*1e11, N) , cmap = cmap)
    ax3.set_title(r'$\gamma_{12} (\times 10 ^{-11} m)$', fontsize = title_size)
    ax3.set_xlabel(r'$x_1$', fontsize = alabel_size)
    ax3.set_ylabel(r'$x_2$', fontsize = alabel_size)
    cbar = plt.colorbar(ax = plt.gca())
    plt.xticks(fontsize= ticks_size)
    plt.yticks(fontsize= ticks_size)
    cbar.ax.tick_params(labelsize=cbar_size)


    ax4 = plt.subplot(3, 2, 2)
    plt.tricontourf(x1, x2, zsig_data[0], levels = np.linspace(tol*np.max(zsig_data[0]), np.max(zsig_data[0]), N) , cmap = cmap)
    ax4.set_title(r'$\sigma_{11} (Pa)$', fontsize = title_size)
    ax4.set_xlabel(r'$x_1$', fontsize = alabel_size)
    ax4.set_ylabel(r'$x_2$', fontsize = alabel_size)
    cbar = plt.colorbar(ax = plt.gca())
    plt.xticks(fontsize= ticks_size)
    plt.yticks(fontsize= ticks_size)
    cbar.ax.tick_params(labelsize=cbar_size)

    ax5 = plt.subplot(3, 2, 4)
    plt.tricontourf(x1, x2, zsig_data[1], levels = np.linspace(tol*np.max(zsig_data[1]), np.max(zsig_data[1]), N) , cmap = cmap)
    ax5.set_title(r'$\sigma_{22} (Pa) $', fontsize = title_size)
    ax5.set_xlabel(r'$x_1$', fontsize = alabel_size)
    ax5.set_ylabel(r'$x_2$', fontsize = alabel_size)
    cbar = plt.colorbar(ax = plt.gca())
    plt.xticks(fontsize= ticks_size)
    plt.yticks(fontsize= ticks_size)
    cbar.ax.tick_params(labelsize=cbar_size)

    ax6 = plt.subplot(3, 2, 6)
    plt.tricontourf(x1, x2, zsig_data[2], levels = np.linspace(tol*np.max(zsig_data[2]), np.max(zsig_data[2]), N) , cmap = cmap)
    ax6.set_title(r'$\sigma_{12} (Pa)$', fontsize = title_size)
    ax6.set_xlabel(r'$x_1$', fontsize = alabel_size)
    ax6.set_ylabel(r'$x_2$', fontsize = alabel_size)
    cbar = plt.colorbar(ax = plt.gca())
    plt.xticks(fontsize= ticks_size)
    plt.yticks(fontsize= ticks_size)
    cbar.ax.tick_params(labelsize=cbar_size)

    plt.subplots_adjust(left=0.1,
                        bottom=0.1, 
                        right=0.9, 
                        top=0.9, 
                        wspace=0.2, 
                        hspace=0.4)

    plt.gcf().set_size_inches(20, 20)

    fname = 'Contour_Plots/contour_panel_'+str(mesh_density_lvl)+'.png'

    # plt.savefig(fname)


def displacement_contour(x1, x2, zdis_data , tol, mesh_density_lvl, cmap='jet'):

    fig = plt.figure()


    alabel_size = 20
    title_size = 25
    ticks_size = 20
    cbar_size = 20

    N  = 200




    ax1 = plt.subplot(1, 2, 1)
    plt.tricontourf(x1, x2, zdis_data[0]*1e9, levels = np.linspace(tol*np.max(zdis_data[0])*1e9, np.max(zdis_data[0])*1e9, N) , cmap = cmap)
    ax1.set_title(r'$d_{x1} (nm)$', fontsize = title_size)
    ax1.set_xlabel(r'$x_1$', fontsize = alabel_size)
    ax1.set_ylabel(r'$x_2$', fontsize = alabel_size)
    cbar = plt.colorbar(ax = plt.gca())
    plt.xticks(fontsize= ticks_size)    
    plt.yticks(fontsize= ticks_size)
    cbar.ax.tick_params(labelsize=cbar_size)

    ax2 = plt.subplot(1, 2, 2)
    plt.tricontourf(x1, x2, zdis_data[1]*1e9, levels = np.linspace(tol*np.max(zdis_data[1])*1e9, np.max(zdis_data[1])*1e9, N) , cmap = cmap)
    ax2.set_title(r'$d_{x2} (nm) $', fontsize = title_size)
    ax2.set_xlabel(r'$x_1$', fontsize = alabel_size)
    ax2.set_ylabel(r'$x_2$', fontsize = alabel_size)
    cbar = plt.colorbar(ax = plt.gca())
    plt.xticks(fontsize= ticks_size)
    plt.yticks(fontsize= ticks_size)
    cbar.ax.tick_params(labelsize=cbar_size)


    plt.subplots_adjust(left=0.1,
                        bottom=0.1, 
                        right=0.9, 
                        top=0.9, 
                        wspace=0.2, 
                        hspace=0.4)

    plt.gcf().set_size_inches(20, 6)

    fname = 'Contour_Plots/contour_panel_d_'+str(mesh_density_lvl)+'.png'

    # plt.savefig(fname)
