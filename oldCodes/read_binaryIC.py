"""

http://memosisland.blogspot.com/2013/04/reading-binary-data-files-written-in-c.html

"""

import numpy as np
import matplotlib.pylab as plt



def load_file(fileHead):
    dt = np.dtype("f4, f4, f4, f4, f4, f4, f4, i8")
    #dt = np.dtype([('x', f4), ('vx', f4), ('y', f4), ('vy', f4), ('z', f4), ('vz', f4), ('phi',f4), ('idx', i8)])
    num_part = 128**3
    tot_entry = 0
    xyz = np.zeros(shape=(3))
    vxvyvz = np.zeros(shape=(3))

    x = []
    y = []
    z = []
    vx = []
    vy = []
    vz = []


    phi = []
    idx = []

    for fidx in range(8):
        np_file =  np.fromfile(fileHead + str(fidx),  dtype=dt)
        tot_entry = tot_entry + np.shape(np_file)[0]
        x = np.append( x, np.array(np_file['f0'])  )
        y = np.append( y, np.array(np_file['f2'])  )
        z = np.append( z, np.array(np_file['f4'])  )
        vx = np.append( vx, np.array(np_file['f1'])  )
        vy = np.append( vy, np.array(np_file['f3'])  )
        vz = np.append( vz, np.array(np_file['f5'])  )
        idx = np.append( idx, np.array(np_file['f7']) )
        phi = np.append( phi, np.array(np_file['f6']) )


    print tot_entry
    print(10*'-')
    return idx, x, y, z, vx, vy, vz






def load_file_za(fileHead):
    dt = np.dtype("f4, f4, f4, f4, f4, f4, f4, i8")
    #dt = np.dtype([('x', f4), ('vx', f4), ('y', f4), ('vy', f4), ('z', f4), ('vz', f4), ('phi', f4), ('idx', i8)])
    num_part = 128**3
    tot_entry = 0

    xyz = np.zeros(shape=(3))
    vxvyvz = np.zeros(shape=(3))

    x = []
    y = []
    z = []
    vx = []
    vy = []
    vz = []
    qx = []
    qy = []
    qz = []
    sx = []
    sy = []
    sz = []



    phi = []
    idx = []

    for fidx in range(8):
        np_file =  np.fromfile(fileHead + str(fidx),  dtype=dt)
        tot_entry = tot_entry + np.shape(np_file)[0]
        x = np.append( x, np.array(np_file['f0'])  )
        y = np.append( y, np.array(np_file['f2'])  )
        z = np.append( z, np.array(np_file['f4'])  )
        vx = np.append( vx, np.array(np_file['f1'])  )
        vy = np.append( vy, np.array(np_file['f3'])  )
        vz = np.append( vz, np.array(np_file['f5'])  )
        
        qx = np.append( vx, np.array(np_file['f1'])  )
        qy = np.append( vy, np.array(np_file['f3'])  )
        qz = np.append( vz, np.array(np_file['f5'])  )

        sx = np.append( vx, np.array(np_file['f1'])  )
        sy = np.append( vy, np.array(np_file['f3'])  )
        sz = np.append( vz, np.array(np_file['f5'])  )
        idx = np.append( idx, np.array(np_file['f7']) )


    print tot_entry

    #print np_file[0]

    print(10*'-')
    #print(np.shape(np_file))
    #print((128**3)*8)

    return idx, x, y, z, vx, vy, vz

def dtfe3d(x, y, z, ngr):
    import sys
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.ndimage.filters import gaussian_filter
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.animation as animation

    from pydtfe import map_dtfe3d 

    rho = map_dtfe3d(x,y,z, ngr) 

    #rho = np.zeros(shape=(ngr, ngr, ngr))
    return rho



#==============================================================
# 4 point Cloud-in-cell method to determine density
def cic3d(x, y, z, L, size_fact):   #  0 <  x, y, z < nGr*refFactor in 1D 


    #x = (x-x_l)*Nx/(x_h - x_l)
    #y = (y-y_l)*Ny/(y_h - y_l) 

    x = x*size_fact/L
    y = y*size_fact/L
    z = z*size_fact/L
    
    Np = np.size(x)
    macro = np.zeros([size_fact, size_fact, size_fact])

    for particle in range(Np):
        i = int(x[particle]) 
        j = int(y[particle]) 
        k = int(z[particle])

	#print 'max i', np.max(i)

        dx = 1
        dy = 1
        dz = 1
        
        a1 = np.around(x[particle], decimals = 4) - i*dx
        b1 = np.around(y[particle], decimals = 4) - j*dy
        c1 = np.around(z[particle], decimals = 4) - k*dz
        
        a2 = dx - a1
        b2 = dy - b1
        c2 = dz - c1
        
        wx1 = a1/dx
        wx2 = a2/dx
        wy1 = b1/dy
        wy2 = b2/dy
        wz1 = c1/dz
        wz2 = c2/dz
        
        #print a1, a2, wx1, wz2
        
        
        macro[i, j, k] += (wx1 * wy1 * wz1)
        macro[np.mod(i+1,size_fact), j, k] += (wx2 * wy1 * wz1)
        macro[i, np.mod(j+1,size_fact), k] += (wx1 * wy2 * wz1)
        macro[np.mod(i+1,size_fact), np.mod(j+1,size_fact), k] += (wx2 * wy2 * wz1)
        
        macro[i, j, np.mod(k+1,size_fact)] += (wx1 * wy1 * wz2)
        macro[np.mod(i+1,size_fact), j, np.mod(k+1,size_fact)] += (wx2 * wy1 * wz2)
        macro[i, np.mod(j+1,size_fact), np.mod(k+1,size_fact)] += (wx1 * wy2 * wz2)
        macro[np.mod(i+1,size_fact), np.mod(j+1,size_fact), np.mod(k+1,size_fact)] += (wx2 * wy2 * wz2)

    return macro

#==============================================================


if __name__== "__main__":
    dirIn = "../HACC/sims/sim932/output/"
    fileIn = "m000.full.hcosmo.IC."

    _, x, y, z, _, _, _ =  load_file(dirIn + fileIn)

    #rho3d = dtfe3d(x, y, z, ngr=128)
    L = 20.
    ngr = 128
    rho3d = cic3d(x, y, z, L, size_fact= ngr)
    np.save('../HACC/Data/cic_'+str(int(L))+'_'+str(ngr)+'.npy', rho3d)

    print x
