"""

http://memosisland.blogspot.com/2013/04/reading-binary-data-files-written-in-c.html

"""

import numpy as np
import matplotlib.pylab as plt


def load_file(fileHead):
    dt = np.dtype("f4, f4, f4, f4, f4, f4, f4, f4, f4, f4, f4, f4, i8")
    # dt = np.dtype([('x', f4), ('vx', f4), ('y', f4), ('vy', f4), ('z', f4), ('vz', f4), ('phi', f4), ('idx', i8)])
    num_part = 128 ** 3
    tot_entry = 0

    # xyz = np.zeros(shape=(3))
    # vxvyvz = np.zeros(shape=(3))
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

    # phi = []
    idx = []

    for fidx in range(8):
        np_file = np.fromfile(fileHead + str(fidx), dtype=dt)
        tot_entry = tot_entry + np.shape(np_file)[0]

        x = np.append(x, np.array(np_file['f0']))
        y = np.append(y, np.array(np_file['f4']))
        z = np.append(z, np.array(np_file['f8']))

        vx = np.append(vx, np.array(np_file['f1']))
        vy = np.append(vy, np.array(np_file['f5']))
        vz = np.append(vz, np.array(np_file['f9']))

        qx = np.append(qx, np.array(np_file['f2']))
        qy = np.append(qy, np.array(np_file['f6']))
        qz = np.append(qz, np.array(np_file['f10']))

        sx = np.append(sx, np.array(np_file['f3']))
        sy = np.append(sy, np.array(np_file['f7']))
        sz = np.append(sz, np.array(np_file['f11']))

        idx = np.append(idx, np.array(np_file['f12']))
        # phi = np.append( phi, np.array(np_file['f6']) )

    print tot_entry

    # print np_file[0]

    print(10 * '-')
    # print(np.shape(np_file))
    # print((128**3)*8)

    return idx, x, y, z, vx, vy, vz, sx, sy, sz, qx, qy, qz


# ==============================================================


if __name__ == "__main__":
    z_arr = [1, 5, 159]
    growth_arr = [0.624913, 0.219980, 0.008534]

    #z_arr = [0.01, 1, 5, 10]
    #growth_arr = [0.995227, 0.624913, 0.219980, 0.120347]




    z0_ind = 0
    z1_ind = 1


    z0 = z_arr[z0_ind]
    dirIn = "../sims/sim_z"+str(z0)+"/output/"
    DataDir = "../HACC/"  # "../"
    # dirIn = DataDir + "sims/simi1Gpc_z" + str(z0) + "/output/"

    # fileIn = "m000.full.hcosmo.IC."

    # _, x, y, z, _, _, _, sx, sy, sz, qx, qy, qz = load_file(dirIn + fileIn)

    # rho3d = dtfe3d(x, y, z, ngr=128)
    L = 20.
    ngr = 128
    # rho3d = cic3d(x, y, z, L, size_fact= ngr)
    rho_z0 = np.load(DataDir + 'Data/cic_'+str(int(L))+'_'+str(ngr)+'_z_' +str(z0) + '.npy')

    dsxdqx = np.load(DataDir + 'Data/dsxdqx_' + str(int(L)) + '_' + str(ngr) + '_z_' + str(z0) + '.npy')
    dsxdqy = np.load(DataDir + 'Data/dsxdqy_' + str(int(L)) + '_' + str(ngr) + '_z_' + str(z0) + '.npy')
    dsxdqz = np.load(DataDir + 'Data/dsxdqz_' + str(int(L)) + '_' + str(ngr) + '_z_' + str(z0) + '.npy')

    dsydqx = np.load(DataDir + 'Data/dsydqx_' + str(int(L)) + '_' + str(ngr) + '_z_' + str(z0) + '.npy')
    dsydqy = np.load(DataDir + 'Data/dsydqy_' + str(int(L)) + '_' + str(ngr) + '_z_' + str(z0) + '.npy')
    dsydqz = np.load(DataDir + 'Data/dsydqz_' + str(int(L)) + '_' + str(ngr) + '_z_' + str(z0) + '.npy')

    dszdqx = np.load(DataDir + 'Data/dszdqx_' + str(int(L)) + '_' + str(ngr) + '_z_' + str(z0) + '.npy')
    dszdqy = np.load(DataDir + 'Data/dszdqy_' + str(int(L)) + '_' + str(ngr) + '_z_' + str(z0) + '.npy')
    dszdqz = np.load(DataDir + 'Data/dszdqz_' + str(int(L)) + '_' + str(ngr) + '_z_' + str(z0) + '.npy')




    rho_z1 = np.zeros( (ngr, ngr, ngr) )
    delta_ij = np.zeros((3, 3), int)
    np.fill_diagonal(delta_ij, 1.0)


    for xi in range(ngr):
        for yi in range(ngr):
            for zi in range(ngr):

                DsDq = np.array([ [dsxdqx[xi, yi, zi], dsydqx[xi, yi, zi], dszdqx[xi, yi, zi]],
                         [dsydqx[xi, yi, zi], dsydqy[xi, yi, zi], dsydqz[xi, yi, zi]],
                         [dszdqx[xi, yi, zi], dszdqy[xi, yi, zi], dszdqz[xi, yi, zi]] ] )



                zel_scaling1 = np.linalg.det( delta_ij +  ((growth_arr[z0_ind]/growth_arr[
                    z1_ind])*(DsDq.T) ))

                numer = np.linalg.det( delta_ij +  (growth_arr[z0_ind]*(DsDq) ) )
                denom = np.linalg.det( delta_ij +  (growth_arr[z1_ind]*(DsDq) ) )
                zel_scaling = numer/denom


                rho_z1[xi, yi, zi] = zel_scaling1*rho_z0[xi, yi, zi]







    ######## PLOTTING ############
    # import matplotlib as mpl
    # mpl.pyplot.viridis()
    # mpl.pyplot.magma()
    # mpl.pyplot.inferno()
    # mpl.pyplot.plasma()
    # mpl.pyplot.summer()



    z1 = z_arr[z1_ind]
    rho_z1_real = np.load(DataDir + 'Data/cic_'+str(int(L))+'_'+str(ngr)+'_z_' +str(z1) + '.npy')

    cmap = plt.get_cmap('afmhot')

    fig = plt.figure(132, figsize=(8,8))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax4 = fig.add_subplot(224)

    ax1.imshow( rho_z0[:,:,53], vmin= 0, vmax= 5, cmap = cmap, alpha = 0.8)
    ax2.imshow( rho_z1_real[:,:,53], vmin= 0, vmax= 5, cmap = cmap, alpha = 0.8)
    # ax4.imshow( rho_z1[:,:,53], vmin= 0, vmax= 5, cmap = cmap , alpha = 0.8)
    ax4.imshow( rho_z1[:,:,53], vmin= 0, vmax= 5, cmap = cmap , alpha = 0.8)

    ax1.set_title(r'$\rho$(z =' + str(z0) + ')')
    ax2.set_title(r'$\rho$(z =' + str(z1) + ')')
    ax4.set_title(r'$\rho$ (ZA mapped)')

    # plt.colorbar(ax1)
    plt.show()


    # from matplotlib.colors import ListedColormap
    #
    # viridis = ListedColormap(_viridis_data, name='viridis')
    #
    # plt.register_cmap(name='viridis', cmap=viridis)
    # plt.set_cmap(viridis)


