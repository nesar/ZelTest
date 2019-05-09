import numpy as np
import matplotlib.pylab as plt



if __name__== "__main__":
    #dirIn = "../sims/sim932_z1/output/"
    #fileIn = "m000.full.hcosmo.IC."

    zin = 5 #1

    #rho3d = dtfe3d(x, y, z, ngr=128)
    L = 20.
    ngr = 128
    cic_file = '../HACC/Data/cic_'+str(int(L))+'_'+str(ngr)+'_z_' +str(zin) + '.npy'

    rho3d = np.load(cic_file)

    
    plt.figure(132)
    #plt.imshow( np.sum(rho3d[:,:,1], axis = 2) , cmap=plt.get_cmap('CMRmap_r', 3) )
    plt.imshow( rho3d[:,:,53] , cmap=plt.get_cmap('CMRmap_r', 30) )
    plt.colorbar()
    plt.show()


    plt.figure(245)
    plt.hist(np.ravel(rho3d), bins = 100)
    plt.show()

