#################################################
###### Authors: Jai Kumar, David Zarzoso ########
#################################################

"""
Class used to defined PiNN method for Vlasov equation
Verify the solution of simulated data.
"""

# pylint: disable=invalid-name
# pylint: disable=no-member
# pylint: disable=import-error

import os
import random

import h5py as h5
import matplotlib.pylab as plt
import numpy as np
import tensorflow as tf
from sklearn.utils import shuffle
from tensorflow.keras.models import load_model

import PINNread as PINNrd
import PINNutils as PINNut

class PINN():
    '''
    Class used for PINN method
    '''
    def __init__(self, inputFile, pathToData):

        self.input_file = inputFile

        # Read the PINN parameters defined in the YAML input file
        input_data = PINNrd.Read_PINN_parameters(inputFile)
        self.__dict__.update(input_data)

        # Create the names of the output files
        self.path_to_data = pathToData
        if pathToData[-1] == '/':
            pathToData = pathToData[:-1]

        # pylint: disable=access-member-before-definition
        path_forPINN_results = os.path.join(self.path_to_results, os.path.basename(pathToData))
        os.makedirs(self.path_to_results, exist_ok=True)
        os.makedirs(path_forPINN_results, exist_ok=True)

        self.path_to_results = path_forPINN_results
        self.path_to_summary = self.path_to_results + "/PINNoutput.dat"
        self.path_to_model = self.path_to_results + "/" + self.model_file
        self.path_to_training_summary = self.path_to_results + "/training_summary.h5"

        # Prepare the variables for PINN analysis
        self.opt = tf.keras.optimizers.Adam(learning_rate=self.learning_rate)
        self.train_loss = []
        self.L1 = []
        self.L2 = []
        self.L3 = []
        self.L4 = []
        self.L5 = []
        self.L6 = []
        self.test_loss = []
        self.no_batches = int(100/self.batch_perc)
        self.batch = 0
        self.epoch = 0
        self.Coeff_1 = []
        self.Coeff_2 = []
        self.Coeff_3 = []

        # Initialise to None all the array that will be used for Indexes
        self.N_fun = None
        self.N_fun_train = None
        self.N_fun_test = None
        self.batch_size_fun = None
        self.indexes_fun = None
        self.ind_t_fun_train = None
        self.ind_x_fun_train = None
        self.ind_v_fun_train = None
        self.ind_t_fun_test = None
        self.ind_x_fun_test = None
        self.ind_v_fun_test = None
        self.ind_t_fun_batch = None
        self.ind_x_fun_batch = None
        self.ind_v_fun_batch = None
        self.tf_felec_fun_test = None
        self.tf_Efield_fun_test = None
        self.tf_t_fun_test = None
        self.tf_x_fun_test = None
        self.tf_v_fun_test = None
        self.tf_felec_fun_train = None
        self.tf_Efield_fun_train = None
        self.tf_t_fun_train = None
        self.tf_x_fun_train = None
        self.tf_v_fun_train = None
        self.model_fun = None

        self.N_PDE = None
        self.N_PDE_train = None
        self.N_PDE_test = None
        self.batch_size_PDE = None
        self.indexes_PDE = None
        self.ind_t_PDE_train = None
        self.ind_x_PDE_train = None
        self.ind_v_PDE_train = None
        self.ind_t_PDE_test = None
        self.ind_x_PDE_test = None
        self.ind_v_PDE_test = None
        self.ind_t_PDE_batch = None
        self.ind_x_PDE_batch = None
        self.ind_v_PDE_batch = None
        self.tf_felec_PDE_train = None
        self.tf_Efield_PDE_train = None
        self.tf_t_PDE_train = None
        self.tf_x_PDE_train = None
        self.tf_v_PDE_train = None
        self.model_PDE = None

        self.N_Efield = None
        self.N_Efield_train = None
        self.N_Efield_test = None
        self.indexes_Efield = None
        self.ind_t_Efield_train = None
        self.ind_x_Efield_train = None
        self.ind_v_Efield_train = None
        self.ind_t_Efield_test = None
        self.ind_x_Efield_test = None
        self.ind_v_Efield_test = None

        self.N_init = None
        self.N_init_train = None
        self.N_init_test = None
        self.batch_size_init = None
        self.indexes_init = None
        self.ind_t_init_train = None
        self.ind_x_init_train = None
        self.ind_v_init_train = None
        self.ind_t_init_test = None
        self.ind_x_init_test = None
        self.ind_v_init_test = None
        self.ind_t_init_batch = None
        self.ind_x_init_batch = None
        self.ind_v_init_batch = None
        self.tf_felec_init_train = None
        self.tf_Efield_init_train = None
        self.tf_t_init_train = None
        self.tf_x_init_train = None
        self.tf_v_init_train = None
        self.model_init = None

        self.N_pbc_0 = None
        self.N_pbc_0_train = None
        self.N_pbc_0_test = None
        self.batch_size_pbc_0 = None
        self.indexes_pbc_0 = None
        self.ind_t_pbc_0_train = None
        self.ind_x_pbc_0_train = None
        self.ind_v_pbc_0_train = None
        self.ind_t_pbc_0_test = None
        self.ind_x_pbc_0_test = None
        self.ind_v_pbc_0_test = None
        self.ind_t_pbc_0_batch = None
        self.ind_x_pbc_0_batch = None
        self.ind_v_pbc_0_batch = None
        self.tf_felec_pbc_0_train = None
        self.tf_Efield_pbc_0_train = None
        self.tf_t_pbc_0_train = None
        self.tf_x_pbc_0_train = None
        self.tf_v_pbc_0_train = None
        self.model_pbc_0 = None

        self.N_pbc_1 = None
        self.N_pbc_1_train = None
        self.N_pbc_1_test = None
        self.batch_size_pbc_1 = None
        self.indexes_pbc_1 = None
        self.ind_t_pbc_1_train = None
        self.ind_x_pbc_1_train = None
        self.ind_v_pbc_1_train = None
        self.ind_t_pbc_1_test = None
        self.ind_x_pbc_1_test = None
        self.ind_v_pbc_1_test = None
        self.ind_t_pbc_1_batch = None
        self.ind_x_pbc_1_batch = None
        self.ind_v_pbc_1_batch = None
        self.tf_felec_pbc_1_train = None
        self.tf_Efield_pbc_1_train = None
        self.tf_t_pbc_1_train = None
        self.tf_x_pbc_1_train = None
        self.tf_v_pbc_1_train = None
        self.model_pbc_1 = None

        self.N_dbc_0 = None
        self.N_dbc_0_train = None
        self.N_dbc_0_test = None
        self.batch_size_dbc_0 = None
        self.indexes_dbc_0 = None
        self.ind_t_dbc_0_train = None
        self.ind_x_dbc_0_train = None
        self.ind_v_dbc_0_train = None
        self.ind_t_dbc_0_test = None
        self.ind_x_dbc_0_test = None
        self.ind_v_dbc_0_test = None
        self.ind_t_dbc_0_batch = None
        self.ind_x_dbc_0_batch = None
        self.ind_v_dbc_0_batch = None
        self.tf_felec_dbc_0_train = None
        self.tf_Efield_dbc_0_train = None
        self.tf_t_dbc_0_train = None
        self.tf_x_dbc_0_train = None
        self.tf_v_dbc_0_train = None
        self.model_dbc_0 = None

        self.N_dbc_1 = None
        self.N_dbc_1_train = None
        self.N_dbc_1_test = None
        self.batch_size_dbc_1 = None
        self.indexes_dbc_1 = None
        self.ind_t_dbc_1_train = None
        self.ind_x_dbc_1_train = None
        self.ind_v_dbc_1_train = None
        self.ind_t_dbc_1_test = None
        self.ind_x_dbc_1_test = None
        self.ind_v_dbc_1_test = None
        self.ind_t_dbc_1_batch = None
        self.ind_x_dbc_1_batch = None
        self.ind_v_dbc_1_batch = None
        self.tf_felec_dbc_1_train = None
        self.tf_Efield_dbc_1_train = None
        self.tf_t_dbc_1_train = None
        self.tf_x_dbc_1_train = None
        self.tf_v_dbc_1_train = None
        self.model_dbc_1 = None

        self.L1_batch = None
        self.L2_batch = None
        self.L3_batch = None
        self.L4_batch = None
        self.L5_batch = None
        self.L6_batch = None
        self.L6_batch = None
        self.loss_batch = None

        self.training_summary_file = None

        # Defining the coefficients for verification
        if self.verify_simu == 1:
            self.Coeff = [tf.Variable(tf.abs(tf.random.normal([1]))) for i in range(2)]
            self.append_Coeff()
        if self.verify_simu == 0:
            self.Coeff = [tf.Variable(1.0), tf.Variable(1.0)]

        if not self.restart:
            self.make_model()
        else:
            self.model = load_model(self.path_to_model)

    #end def __init__


    #-----------------------------------------------
    # Load VOICEXX results
    #-----------------------------------------------
    def load_VOICEXX_resu(self):
        '''
        Load VOICEXX results
        '''
        resudir = self.path_to_data
        pinn = PINNrd.Read_PINN_fromVOICEXX(resudir)
        self.__dict__.update(pinn.__dict__)


    #-----------------------------------------------
    # Build the indexes to define which data
    #  will be used in PINN
    #-----------------------------------------------
    def buildIndexes(self):
        '''
        Build the indexes to define which data will be used in PINN
        '''

        def indexes(it1, it2, ix1, ix2, iv1, iv2):
            ind_t = np.arange(it1, it2)
            ind_x = np.arange(ix1, ix2)
            ind_v = np.arange(iv1, iv2)
            ind_tt, ind_xx, ind_vv = np.meshgrid(ind_t, ind_x, ind_v, indexing='ij')
            ind_tt_ = ind_tt.flatten()
            ind_xx_ = ind_xx.flatten()
            ind_vv_ = ind_vv.flatten()
            myIndex = (ind_tt_, ind_xx_, ind_vv_)
            lenIndex = len(ind_tt_)
            return lenIndex, myIndex

        def indexes_E(it1, it2, ix1, ix2):
            ind_t = np.arange(it1, it2)
            ind_x = np.arange(ix1, ix2)
            ind_tt, ind_xx = np.meshgrid(ind_t, ind_x, indexing='ij')
            ind_tt_ = ind_tt.flatten()
            ind_xx_ = ind_xx.flatten()
            myIndex = (ind_tt_, ind_xx_)
            lenIndex = len(ind_tt_)
            return lenIndex, myIndex

        # Function indexes in the inner region
        self.N_fun, self.indexes_fun = indexes(1, self.nt,
                                               1, self.nx - 1,
                                               1, self.nv - 1)

        # PDE or Efield indexes depending on use_fun_data=0 or 1
        if self.use_fun_data == 1:
            #--> PDE indexes in the inner region
            self.N_PDE, self.indexes_PDE = indexes(1, self.nt,
                                                   1, self.nx - 1,
                                                   1, self.nv - 1)
        if self.use_fun_data == 0:
            #--> Efield indexes in the inner region
            self.N_Efield, self.indexes_Efield = indexes_E(1, self.nt,
                                                           1, self.nx - 1)

        # Initial condition (IC) indexes
        self.N_init, self.indexes_init = indexes(0, 1,
                                                 0, self.nx,
                                                 0, self.nv)

        #VG# [TODO] : ATTENTION verifier condition de periodicité dans VOICEXX
        #VG#  --> Est-ce que le point periodique est sauve ?
        # Periodic boundary conditions (PBC) indexes
        self.N_pbc_0, self.indexes_pbc_0 = indexes(1, self.nt,
                                                   0, 1,
                                                   1, self.nv - 1)
        self.N_pbc_1, self.indexes_pbc_1 = indexes(1, self.nt,
                                                   self.nx - 1, self.nx,
                                                   1, self.nv - 1)

        # Dirichtlet boundary conditions (DBC) indexes, where f(v_i) = 0, f(v_f) = 0
        self.N_dbc_0, self.indexes_dbc_0 = indexes(1, self.nt,
                                                   1, self.nx,
                                                   0, 1)
        self.N_dbc_1, self.indexes_dbc_1 = indexes(1, self.nt,
                                                   1, self.nx,
                                                   self.nv-1, self.nv)
    #end def buildIndexes


    #-----------------------------------------------
    # Select randomly the indexes for data used
    #  for train and test
    #-----------------------------------------------
    def selectRandomIndexesTrainTest(self):
        '''
        Select randomly the indexes for data used for train and test
        '''

        def splitTrainTest(indexes, N, perc, train_perc):
            Ntotal = (N*perc)//100
            Ntrain = (Ntotal*train_perc)//100
            Ntest = (Ntotal*(100-train_perc))//100
            Nsel = Ntrain + Ntest
            (shuffled_indexes_0,
             shuffled_indexes_1,
             shuffled_indexes_2) = shuffle(indexes[0],
                                           indexes[1],
                                           indexes[2],
                                           random_state=0)
            (ind_t_train, ind_x_train,
             ind_v_train) = \
                 (shuffled_indexes_0[:Ntrain],
                  shuffled_indexes_1[:Ntrain],
                  shuffled_indexes_2[:Ntrain])
            (ind_t_test, ind_x_test,
             ind_v_test) = \
                 (shuffled_indexes_0[Ntrain:Nsel],
                  shuffled_indexes_1[Ntrain:Nsel],
                  shuffled_indexes_2[Ntrain:Nsel])
            return Ntrain, Ntest, \
                ind_t_train, ind_x_train, ind_v_train, \
                ind_t_test, ind_x_test, ind_v_test


        def splitTrainTest_E(indexes, N, x_t_perc, v_perc, train_perc):
            Ntotal = (N*x_t_perc)//100
            Ntrain = (Ntotal*train_perc)//100
            Ntest = (Ntotal*(100-train_perc))//100
            Nsel = Ntrain + Ntest
            Nvel = (self.nv*v_perc)//100
            ind_t_train = []
            ind_x_train = []
            ind_v_train = []
            ind_t_test = []
            ind_x_test = []
            ind_v_test = []
            (shuffled_indexes_0,
             shuffled_indexes_1) = shuffle(indexes[0],
                                           indexes[1],
                                           random_state=0)
            (ind_t_train_, ind_x_train_) = \
                 (shuffled_indexes_0[:Ntrain],
                  shuffled_indexes_1[:Ntrain])

            for (ti, xi) in zip(ind_t_train_, ind_x_train_):
                ind_t_train.append(np.repeat(ti, Nvel))
                ind_x_train.append(np.repeat(xi, Nvel))
                ind_v_train.append(random.sample(sorted(np.arange(1, self.nv)), Nvel))
            ind_t_train = np.array(ind_t_train).flatten()
            ind_x_train = np.array(ind_x_train).flatten()
            ind_v_train = np.array(ind_v_train).flatten()

            (ind_t_test_, ind_x_test_) = \
                 (shuffled_indexes_0[Ntrain:Nsel],
                  shuffled_indexes_1[Ntrain:Nsel])

            for (ti, xi) in zip(ind_t_test_, ind_x_test_):
                ind_t_test.append(np.repeat(ti, Nvel))
                ind_x_test.append(np.repeat(xi, Nvel))
                ind_v_test.append(random.sample(sorted(np.arange(1, self.nv)), Nvel))
            ind_t_test = np.array(ind_t_test).flatten()
            ind_x_test = np.array(ind_x_test).flatten()
            ind_v_test = np.array(ind_v_test).flatten()

            return [Ntrain, Ntest,
                    ind_t_train, ind_x_train, ind_v_train,
                    ind_t_test, ind_x_test, ind_v_test]

        [self.N_fun_train, self.N_fun_test,
         self.ind_t_fun_train, self.ind_x_fun_train, self.ind_v_fun_train,
         self.ind_t_fun_test, self.ind_x_fun_test, self.ind_v_fun_test] = \
             splitTrainTest(self.indexes_fun, self.N_fun,
                            self.fun_perc, self.train_perc)

        if self.use_fun_data == 1:
            [self.N_PDE_train, self.N_PDE_test,
             self.ind_t_PDE_train, self.ind_x_PDE_train, self.ind_v_PDE_train,
             self.ind_t_PDE_test, self.ind_x_PDE_test, self.ind_v_PDE_test] = \
                 splitTrainTest(self.indexes_PDE, self.N_PDE,
                                self.PDE_perc, self.train_perc)

        if self.use_fun_data == 0:
            [self.N_Efield_train, self.N_Efield_test,
             self.ind_t_Efield_train, self.ind_x_Efield_train,
             self.ind_v_Efield_train, self.ind_t_Efield_test,
             self.ind_x_Efield_test, self.ind_v_Efield_test] = \
                 splitTrainTest_E(self.indexes_Efield, self.N_Efield,
                                  self.x_t_data_perc, self.v_data_perc, self.train_perc)

            self.N_PDE_train = len(self.ind_t_Efield_train)
            self.N_PDE_test = len(self.ind_t_Efield_test)

            (self.ind_t_PDE_train,
             self.ind_x_PDE_train,
             self.ind_v_PDE_train,
             self.ind_t_PDE_test,
             self.ind_x_PDE_test,
             self.ind_v_PDE_test) = (self.ind_t_Efield_train,
                                     self.ind_x_Efield_train,
                                     self.ind_v_Efield_train,
                                     self.ind_t_Efield_test,
                                     self.ind_x_Efield_test,
                                     self.ind_v_Efield_test)

        [self.N_init_train, self.N_init_test,
         self.ind_t_init_train, self.ind_x_init_train, self.ind_v_init_train,
         self.ind_t_init_test, self.ind_x_init_test, self.ind_v_init_test] = \
             splitTrainTest(self.indexes_init, self.N_init,
                            self.init_perc, self.train_perc)

        [self.N_pbc_0_train, self.N_pbc_0_test,
         self.ind_t_pbc_0_train, self.ind_x_pbc_0_train, self.ind_v_pbc_0_train,
         self.ind_t_pbc_0_test, self.ind_x_pbc_0_test, self.ind_v_pbc_0_test] = \
             splitTrainTest(self.indexes_pbc_0, self.N_pbc_0,
                            self.pbc_0_perc, self.train_perc)
        [self.N_pbc_1_train, self.N_pbc_1_test,
         self.ind_t_pbc_1_train, self.ind_x_pbc_1_train, self.ind_v_pbc_1_train,
         self.ind_t_pbc_1_test, self.ind_x_pbc_1_test, self.ind_v_pbc_1_test] = \
             splitTrainTest(self.indexes_pbc_1, self.N_pbc_1,
                            self.pbc_1_perc, self.train_perc)
        [self.N_dbc_0_train, self.N_dbc_0_test,
         self.ind_t_dbc_0_train, self.ind_x_dbc_0_train, self.ind_v_dbc_0_train,
         self.ind_t_dbc_0_test, self.ind_x_dbc_0_test, self.ind_v_dbc_0_test] = \
             splitTrainTest(self.indexes_dbc_0, self.N_dbc_0,
                            self.dbc_0_perc, self.train_perc)
        [self.N_dbc_1_train, self.N_dbc_1_test,
         self.ind_t_dbc_1_train, self.ind_x_dbc_1_train, self.ind_v_dbc_1_train,
         self.ind_t_dbc_1_test, self.ind_x_dbc_1_test, self.ind_v_dbc_1_test] = \
             splitTrainTest(self.indexes_dbc_1, self.N_dbc_1,
                            self.dbc_1_perc, self.train_perc)

        print("No. of points for training:")
        print("      To fit d.f: ", self.N_fun_train)
        print("      To fit PDE: ", self.N_PDE_train)
        print("      To fit I.C: ", self.N_init_train)
        print("      To fit B.C-1: ", self.N_pbc_0_train)
        print("      To fit B.C-1: ", self.N_pbc_1_train)
        print("      To fit B.C-2: ", self.N_dbc_0_train)
        print("      To fit B.C-2: ", self.N_dbc_1_train)
        print("No. of points for testing:")
        print("      To fit d.f: ", self.N_fun_test)

        self.batch_size_fun = self.N_fun_train//self.no_batches
        self.batch_size_PDE = self.N_PDE_train//self.no_batches
        self.batch_size_init = self.N_init_train//self.no_batches
        self.batch_size_pbc_0 = self.N_pbc_0_train//self.no_batches
        self.batch_size_pbc_1 = self.N_pbc_1_train//self.no_batches
        self.batch_size_dbc_0 = self.N_dbc_0_train//self.no_batches
        self.batch_size_dbc_1 = self.N_dbc_1_train//self.no_batches

    #end def selectRandomIndexesTrainTest


    #-----------------------------------------------
    # Select randomly the indexes of data for train
    #-----------------------------------------------
    def randomizeTrainIndexes(self):
        '''
        Select randomly the indexes of data for train
        '''
        (self.ind_t_fun_train, self.ind_x_fun_train,
         self.ind_v_fun_train) = shuffle(self.ind_t_fun_train,
                                         self.ind_x_fun_train,
                                         self.ind_v_fun_train,
                                         random_state=0)
        (self.ind_t_PDE_train, self.ind_x_PDE_train,
         self.ind_v_PDE_train) = shuffle(self.ind_t_PDE_train,
                                         self.ind_x_PDE_train,
                                         self.ind_v_PDE_train,
                                         random_state=0)
        (self.ind_t_init_train, self.ind_x_init_train,
         self.ind_v_init_train) = shuffle(self.ind_t_init_train,
                                          self.ind_x_init_train,
                                          self.ind_v_init_train,
                                          random_state=0)
        (self.ind_t_pbc_0_train, self.ind_x_pbc_0_train,
         self.ind_v_pbc_0_train) = shuffle(self.ind_t_pbc_0_train,
                                           self.ind_x_pbc_0_train,
                                           self.ind_v_pbc_0_train,
                                           random_state=0)
        (self.ind_t_pbc_1_train, self.ind_x_pbc_1_train,
         self.ind_v_pbc_1_train) = shuffle(self.ind_t_pbc_1_train,
                                           self.ind_x_pbc_1_train,
                                           self.ind_v_pbc_1_train,
                                           random_state=0)
        (self.ind_t_dbc_0_train, self.ind_x_dbc_0_train,
         self.ind_v_dbc_0_train) = shuffle(self.ind_t_dbc_0_train,
                                           self.ind_x_dbc_0_train,
                                           self.ind_v_dbc_0_train,
                                           random_state=0)
        (self.ind_t_dbc_1_train, self.ind_x_dbc_1_train,
         self.ind_v_dbc_1_train) = shuffle(self.ind_t_dbc_1_train,
                                           self.ind_x_dbc_1_train,
                                           self.ind_v_dbc_1_train,
                                           random_state=0)
    #end def randomizeTrainIndexes


    #-----------------------------------------------
    # Get all the indexes depending of the batch
    #-----------------------------------------------
    def getBatch(self):
        '''
        Get all the indexes depending of the batch
        '''
        myInd = np.arange(self.batch_size_fun*self.batch,
                          self.batch_size_fun*(self.batch+1))
        self.ind_t_fun_batch = self.ind_t_fun_train[myInd]
        self.ind_x_fun_batch = self.ind_x_fun_train[myInd]
        self.ind_v_fun_batch = self.ind_v_fun_train[myInd]

        myInd = np.arange(self.batch_size_PDE*self.batch,
                          self.batch_size_PDE*(self.batch+1))
        self.ind_t_PDE_batch = self.ind_t_PDE_train[myInd]
        self.ind_x_PDE_batch = self.ind_x_PDE_train[myInd]
        self.ind_v_PDE_batch = self.ind_v_PDE_train[myInd]

        myInd = np.arange(self.batch_size_init*self.batch,
                          self.batch_size_init*(self.batch+1))
        self.ind_t_init_batch = self.ind_t_init_train[myInd]
        self.ind_x_init_batch = self.ind_x_init_train[myInd]
        self.ind_v_init_batch = self.ind_v_init_train[myInd]

        myInd = np.arange(self.batch_size_pbc_0*self.batch,
                          self.batch_size_pbc_1*(self.batch+1))
        self.ind_t_pbc_0_batch = self.ind_t_pbc_0_train[myInd]
        self.ind_x_pbc_0_batch = self.ind_x_pbc_0_train[myInd]
        self.ind_v_pbc_0_batch = self.ind_v_pbc_0_train[myInd]

        myInd = np.arange(self.batch_size_pbc_1*self.batch,
                          self.batch_size_pbc_1*(self.batch+1))
        self.ind_t_pbc_1_batch = self.ind_t_pbc_1_train[myInd]
        self.ind_x_pbc_1_batch = self.ind_x_pbc_1_train[myInd]
        self.ind_v_pbc_1_batch = self.ind_v_pbc_1_train[myInd]

        myInd = np.arange(self.batch_size_dbc_0*self.batch,
                          self.batch_size_dbc_0*(self.batch+1))
        self.ind_t_dbc_0_batch = self.ind_t_dbc_0_train[myInd]
        self.ind_x_dbc_0_batch = self.ind_x_dbc_0_train[myInd]
        self.ind_v_dbc_0_batch = self.ind_v_dbc_0_train[myInd]

        myInd = np.arange(self.batch_size_dbc_1*self.batch,
                          self.batch_size_dbc_1*(self.batch+1))
        self.ind_t_dbc_1_batch = self.ind_t_dbc_1_train[myInd]
        self.ind_x_dbc_1_batch = self.ind_x_dbc_1_train[myInd]
        self.ind_v_dbc_1_batch = self.ind_v_dbc_1_train[myInd]

    #end def getBatch


    #-----------------------------------------------
    # Create TensorFlow variable from the test data
    #-----------------------------------------------
    def createTFvblesTest(self):
        '''
        Create TensorFlow variable from the test data
        '''
        ind_t = self.ind_t_fun_test
        ind_x = self.ind_x_fun_test
        ind_v = self.ind_v_fun_test
        self.tf_felec_fun_test = tf.Variable(
            self.felec[(ind_t, ind_x, ind_v)].flatten()[:, None], dtype=tf.float32)
        self.tf_Efield_fun_test = tf.Variable(
            self.Efield[(ind_t, ind_x)].flatten()[:, None], dtype=tf.float32)

        self.tf_t_fun_test = tf.Variable(
            self.timegrid[(ind_t)].flatten()[:, None], dtype=tf.float32)
        self.tf_x_fun_test = tf.Variable(
            self.xgrid[(ind_x)].flatten()[:, None], dtype=tf.float32)
        self.tf_v_fun_test = tf.Variable(
            self.vgrid[(ind_v)].flatten()[:, None], dtype=tf.float32)

    #end def createTFvblesTest


    #-----------------------------------------------
    # Create TensorFlow variables from the train data evaluated
    # at the randomly selected indexes
    #-----------------------------------------------
    def createTFvblesTrain(self):
        '''
        Create TensorFlow variables from the train data evaluated
        at the randomly selected indexes
        '''
        def createTFvbles(felec, Efield, timegrid, xgrid, vgrid,
                          ind_t, ind_x, ind_v):
            tf_felec = tf.Variable(
                felec[(ind_t, ind_x, ind_v)].flatten()[:, None], dtype=tf.float32)
            tf_Efield = tf.Variable(
                Efield[(ind_t, ind_x)].flatten()[:, None], dtype=tf.float32)
            tf_t = tf.Variable(
                timegrid[(ind_t)].flatten()[:, None], dtype=tf.float32)
            tf_x = tf.Variable(
                xgrid[(ind_x)].flatten()[:, None], dtype=tf.float32)
            tf_v = tf.Variable(
                vgrid[(ind_v)].flatten()[:, None], dtype=tf.float32)
            return [tf_felec, tf_Efield, tf_t, tf_x, tf_v]

        [self.tf_felec_fun_train, self.tf_Efield_fun_train,
         self.tf_t_fun_train, self.tf_x_fun_train, self.tf_v_fun_train] = \
             createTFvbles(self.felec, self.Efield, self.timegrid,
                           self.xgrid, self.vgrid, self.ind_t_fun_batch,
                           self.ind_x_fun_batch, self.ind_v_fun_batch)

        [self.tf_felec_PDE_train, self.tf_Efield_PDE_train,
         self.tf_t_PDE_train, self.tf_x_PDE_train, self.tf_v_PDE_train] = \
             createTFvbles(self.felec, self.Efield, self.timegrid,
                           self.xgrid, self.vgrid, self.ind_t_PDE_batch,
                           self.ind_x_PDE_batch, self.ind_v_PDE_batch)
        [self.tf_felec_init_train, self.tf_Efield_init_train,
         self.tf_t_init_train, self.tf_x_init_train, self.tf_v_init_train] = \
             createTFvbles(self.felec, self.Efield, self.timegrid,
                           self.xgrid, self.vgrid, self.ind_t_init_batch,
                           self.ind_x_init_batch, self.ind_v_init_batch)
        [self.tf_felec_pbc_0_train, self.tf_Efield_pbc_0_train,
         self.tf_t_pbc_0_train, self.tf_x_pbc_0_train, self.tf_v_pbc_0_train] = \
             createTFvbles(self.felec, self.Efield, self.timegrid,
                           self.xgrid, self.vgrid, self.ind_t_pbc_0_batch,
                           self.ind_x_pbc_0_batch, self.ind_v_pbc_0_batch)
        [self.tf_felec_pbc_1_train, self.tf_Efield_pbc_1_train,
         self.tf_t_pbc_1_train, self.tf_x_pbc_1_train, self.tf_v_pbc_1_train] = \
             createTFvbles(self.felec, self.Efield, self.timegrid,
                           self.xgrid, self.vgrid, self.ind_t_pbc_1_batch,
                           self.ind_x_pbc_1_batch, self.ind_v_pbc_1_batch)
        [self.tf_felec_dbc_0_train, self.tf_Efield_dbc_0_train,
         self.tf_t_dbc_0_train, self.tf_x_dbc_0_train, self.tf_v_dbc_0_train] = \
             createTFvbles(self.felec, self.Efield, self.timegrid,
                           self.xgrid, self.vgrid, self.ind_t_dbc_0_batch,
                           self.ind_x_dbc_0_batch, self.ind_v_dbc_0_batch)
        [self.tf_felec_dbc_1_train, self.tf_Efield_dbc_1_train,
         self.tf_t_dbc_1_train, self.tf_x_dbc_1_train, self.tf_v_dbc_1_train] = \
             createTFvbles(self.felec, self.Efield, self.timegrid,
                           self.xgrid, self.vgrid, self.ind_t_dbc_1_batch,
                           self.ind_x_dbc_1_batch, self.ind_v_dbc_1_batch)

    #end def createTFvblesTrain


    #-----------------------------------------------
    # Define the PINN model
    #-----------------------------------------------
    def make_model(self):
        '''
        Define the PINN model
        '''
        self.model = PINNut.myModel(hidden_dim=self.layers,
                                    input_dim=self.input_dim,
                                    output_dim=self.output_dim,
                                    hidden_act_func=self.hidden_act_func,
                                    output_act_func=self.output_act_func,
                                    file_for_summary=self.path_to_summary)


    #-----------------------------------------------
    # Evaluate the PINN model
    #-----------------------------------------------
    def eval_model(self, t, x, v):
        '''
        Evaluate the PINN model
        '''
        u = tf.concat([t, x, v], axis=1)
        return self.model(u)


    #-----------------------------------------------
    # Save the PINN model
    #-----------------------------------------------
    def save_model(self):
        '''
        Save the PINN model
        '''
        self.model.save(self.path_to_model)


    #-----------------------------------------------
    # Construct the physical model depending on
    #  the equation to be treated
    #-----------------------------------------------
    def PhysFun(self):
        '''
        Construct the physical model depending on the equation to be treated
        '''

        with tf.GradientTape(persistent=True) as tape:
            model_txv = self.eval_model(self.tf_t_PDE_train,
                                        self.tf_x_PDE_train,
                                        self.tf_v_PDE_train)

        dmodel_dt = tape.gradient(model_txv, self.tf_t_PDE_train)
        dmodel_dx = tape.gradient(model_txv, self.tf_x_PDE_train)
        dmodel_dv = tape.gradient(model_txv, self.tf_v_PDE_train)

        if self.use_fun_data == 1:
            self.model_fun = self.eval_model(self.tf_t_fun_train,
                                             self.tf_x_fun_train,
                                             self.tf_v_fun_train)
        self.model_init = self.eval_model(self.tf_t_init_train,
                                          self.tf_x_init_train,
                                          self.tf_v_init_train)
        self.model_pbc_0 = self.eval_model(self.tf_t_pbc_0_train,
                                           self.tf_x_pbc_0_train,
                                           self.tf_v_pbc_0_train)
        self.model_pbc_1 = self.eval_model(self.tf_t_pbc_1_train,
                                           self.tf_x_pbc_1_train,
                                           self.tf_v_pbc_1_train)
        self.model_dbc_0 = self.eval_model(self.tf_t_dbc_0_train,
                                           self.tf_x_dbc_0_train,
                                           self.tf_v_dbc_0_train)
        self.model_dbc_1 = self.eval_model(self.tf_t_dbc_1_train,
                                           self.tf_x_dbc_1_train,
                                           self.tf_v_dbc_1_train)
        del tape

        #VG# [TODO] Is this computation useful if verify_simu=1 ?
        C1, C2 = self.Coeff
        self.model_PDE = dmodel_dt + \
            C1*self.tf_v_PDE_train*dmodel_dx - \
            C2*self.tf_Efield_PDE_train*dmodel_dv

    #end def PhysFun


    #-----------------------------------------------
    # Compute losses associated to the batch
    #-----------------------------------------------
    def lossBatch(self):
        '''
        Compute the loss associated to the batch
        '''
        self.L1_batch = tf.reduce_mean(tf.square(self.model_PDE))
        self.L2_batch = tf.reduce_mean(tf.square(self.tf_felec_init_train - self.model_init))
        self.L3_batch = tf.reduce_mean(tf.square(self.model_pbc_0 - self.model_pbc_1))
        self.L4_batch = tf.reduce_mean(tf.square(self.model_dbc_0))
        self.L5_batch = tf.reduce_mean(tf.square(self.model_dbc_1))
        if self.use_fun_data == 1:
            self.L6_batch = tf.reduce_mean(tf.square(self.tf_felec_fun_train - self.model_fun))
        if self.use_fun_data == 0:
            self.L6_batch = tf.constant(0.0)
        self.loss_batch = self.L1_batch + self.L2_batch + \
            self.L3_batch + self.L4_batch + self.L5_batch + \
            self.L6_batch


    #-----------------------------------------------
    # Save losses associated to each epoch
    # by adding losses of batches
    #-----------------------------------------------
    def lossEpoch(self):
        '''
        Save losses associated to each epoch by adding losses of batches
        '''
        self.L1.append(self.L1_batch)
        self.L2.append(self.L2_batch)
        self.L3.append(self.L3_batch)
        self.L4.append(self.L4_batch)
        self.L5.append(self.L5_batch)
        self.L6.append(self.L6_batch)
        self.train_loss.append(self.loss_batch)


    #-----------------------------------------------
    # Train the PINN model
    #-----------------------------------------------
    def train(self):
        '''
        Train the PINN model
        '''
        with tf.GradientTape(persistent=True) as tape_train:
            self.PhysFun()
            self.lossBatch()
            train_loss = self.loss_batch

        grad_train_loss = tape_train.gradient(train_loss, self.model.trainable_variables)
        self.opt.apply_gradients(zip(grad_train_loss, self.model.trainable_variables))

        if self.verify_simu == 1:
            grad_train_Coeff = tape_train.gradient(train_loss, self.Coeff)
            self.opt.apply_gradients(zip(grad_train_Coeff, self.Coeff))

        del tape_train


    #-----------------------------------------------
    # Test the PINN model
    #-----------------------------------------------
    def test(self):
        '''
        Test the PINN model
        '''
        self.test_loss.append(
            tf.reduce_mean(tf.square(
                self.tf_felec_fun_test-self.eval_model(self.tf_t_fun_test,
                                                       self.tf_x_fun_test,
                                                       self.tf_v_fun_test))))


    #-----------------------------------------------
    # Save the training
    #-----------------------------------------------
    def save_training(self):
        '''
        Save the training
        '''

        print('Save training in :', self.path_to_training_summary)
        self.training_summary_file = h5.File(self.path_to_training_summary, 'w')
        PINNut.save_data(self.training_summary_file, 'L1',
                         [dum.numpy() for dum in self.L1])
        PINNut.save_data(self.training_summary_file, 'L2',
                         [dum.numpy() for dum in self.L2])
        PINNut.save_data(self.training_summary_file, 'L3',
                         [dum.numpy() for dum in self.L3])
        PINNut.save_data(self.training_summary_file, 'L4',
                         [dum.numpy() for dum in self.L4])
        PINNut.save_data(self.training_summary_file, 'L5',
                         [dum.numpy() for dum in self.L5])
        PINNut.save_data(self.training_summary_file, 'L6',
                         [dum.numpy() for dum in self.L6])
        PINNut.save_data(self.training_summary_file, 'mse_test',
                         [dum.numpy() for dum in self.test_loss])
        self.training_summary_file.close()


    #-----------------------------------------------
    # Add the computation of the coefficients
    # which should correspond to the coefficients
    # of each operators in the Vlasov equation
    #-----------------------------------------------
    def append_Coeff(self):
        '''
        Add the computation of the coefficients if verify_simu = 1
        They should correspond to the coefficients of each operators in Vlasov equation
        '''
        self.Coeff_1.append(self.Coeff[0].numpy()[0])
        self.Coeff_2.append(self.Coeff[1].numpy()[0])

    #end def append_Coeff


    #-----------------------------------------------
    # Plot coefficient values which should be equal
    # to see if the different operators of the
    # equations are well solved
    #-----------------------------------------------
    def plot_Coeff(self):
        '''
        Plot coefficient values which should be equal
        to see if the different operators of the
        equations are well solved
        '''
        plt.figure(dpi=200)
        epochs = np.arange(0, len(self.Coeff_1))*self.freq_save_data
        plt.scatter(epochs, self.Coeff_1, label='C:1', s=0.5)
        plt.scatter(epochs, self.Coeff_2, label='C:2', s=0.5)
        plt.legend()
        plt.xlabel('Coefficient Value')
        plt.ylabel('no. epochs')
        figureFile = self.path_to_results+'/Coefficient_train.png'
        print("Figure for Coefficient train saved in :", figureFile)
        plt.savefig(figureFile)
        plt.close()


    #-----------------------------------------------
    # Plot the evolution of losses
    #-----------------------------------------------
    def plot_losses(self):
        '''
        Plot the evolution of losses
        '''
        train_loss = np.asarray([tl.numpy() for tl in self.train_loss])
        test_loss = np.asarray([tl.numpy() for tl in self.test_loss])
        L1 = np.asarray([l.numpy() for l in self.L1])
        L2 = np.asarray([l.numpy() for l in self.L2])
        L3 = np.asarray([l.numpy() for l in self.L3])
        L4 = np.asarray([l.numpy() for l in self.L4])
        L5 = np.asarray([l.numpy() for l in self.L5])
        L6 = np.asarray([l.numpy() for l in self.L6])

        epochs = np.arange(0, len(self.train_loss))*self.freq_save_data

        plt.figure(dpi=200)
        plt.yscale('log')
        plt.plot(epochs, test_loss, label='test error')
        plt.grid()
        plt.plot(epochs, train_loss, label='train error')
        plt.plot(epochs, L1, label='PDE')
        plt.plot(epochs, L2, label='IC')
        plt.plot(epochs, L3, label='PBC')
        plt.plot(epochs, L4, label='DBC@vmin')
        plt.plot(epochs, L5, label='DBC@vmax')
        plt.plot(epochs, L6, label='Dist. function')
        plt.legend()
        plt.xlabel('no. epoch')
        plt.ylabel('Loss')

        figureFile = self.path_to_results+'/Losses.png'
        print("Figure for losses saved in :", figureFile)
        plt.savefig(figureFile)
        plt.close()

    #end def plot_losses


    #-----------------------------------------------
    # Plot the predicted results of PINN
    #-----------------------------------------------
    def plot_prediction(self):
        '''
        Plot the predicted results of PINN
        '''
        ix_array = np.linspace(self.nx//8, self.nx, 4, dtype=int, endpoint=False)
        it_array = np.linspace(0, self.nt-1, 4, dtype=int, endpoint=True)
        iv_array = np.arange(0, self.nv)
        style_pred = ['b-', 'c-', 'r-', 'k-']
        style_true = ['b--', 'c--', 'r--', 'k--']
        _, axs = plt.subplots(1, len(ix_array), figsize=(30, 8))
        for iix, ix in enumerate(ix_array):
            for iit, it in enumerate(it_array):
                ind_tt, ind_xx, ind_vv = np.meshgrid(it, ix, iv_array, indexing='ij')
                ind_tt_ = ind_tt.flatten()
                ind_xx_ = ind_xx.flatten()
                ind_vv_ = ind_vv.flatten()
                tf_t = tf.Variable(
                    self.timegrid[(ind_tt_)].flatten()[:, None], dtype=tf.float32)
                tf_x = tf.Variable(
                    self.xgrid[(ind_xx_)].flatten()[:, None], dtype=tf.float32)
                tf_v = tf.Variable(
                    self.vgrid[(ind_vv_)].flatten()[:, None], dtype=tf.float32)
                tf_prediction = self.eval_model(tf_t, tf_x, tf_v)
                axs[iix].plot(self.vgrid,
                              self.felec[it, ix, :],
                              style_true[iit])
                label_t = 't = ' + str(self.timegrid[it])
                axs[iix].plot(self.vgrid[0:self.nv],
                              tf_prediction.numpy(),
                              style_pred[iit],
                              label=label_t)
                subtitle = 'x = '+str(round(self.xgrid[ix], 2))
                axs[iix].set_title(subtitle, fontsize=16)
                axs[iix].set_xlabel('vel', fontsize=16)
                if iix == 0:
                    axs[iix].set_ylabel('Dist. function', fontsize=16)

                axs[iix].legend()

        fig_title = f"/Prediction_distribfunc_epoch-{self.epoch}.png"
        figureFile = self.path_to_results+fig_title
        print("Figure for prediction saved in :", figureFile)
        plt.savefig(figureFile)
        plt.close()
