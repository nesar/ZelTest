import numpy as np
import matplotlib.pylab as plt



if __name__== "__main__":
    dirIn = "../HACC/sims/sim932/output/"
    fileIn = "m000.full.hcosmo.IC."



    #rho3d = dtfe3d(x, y, z, ngr=128)
    L = 20.
    ngr = 128
    cic_file = '../HACC/Data/cic_'+str(int(L))+'_'+str(ngr)+'.npy'

    rho3d = np.load(cic_file)

    
    plt.figure(132)
    #plt.imshow( np.sum(rho3d[:,:,1], axis = 2) , cmap=plt.get_cmap('CMRmap_r', 3) )
    plt.imshow( rho3d[:,:,1] , cmap=plt.get_cmap('CMRmap_r', 3) )
    plt.colorbar()
    plt.show()


    plt.figure(245)
    plt.hist(np.ravel(rho3d), bins = 4)
    plt.show()
