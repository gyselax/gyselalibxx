'''
Main program to apply PINN method to Vlasov-Poisson equation
'''

# pylint: disable=import-error
# pylint: disable=invalid-name
# pylint: disable=no-member

import getopt
import sys
import time

import PINNclass as PINNcl

#===============================================
# Read the arguments given in the command line
#-----------------------------------------------
def read_args(argv):
    '''
    Read the arguments given in the command line
    '''

    try:
        opts, _ = getopt.getopt(argv, "hi:p:",
                                ["ifile=", "pathtodata="])
    except getopt.GetoptError:
        print('python3 PINNmain.py -i <inputfile> -p <pathtodata>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('python3 PINNmain.py -i <inputfile> -p <pathtodata>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            arg_inputFile = arg
        elif opt in ("-p", "--pathtodata"):
            arg_pathToData = arg

    return [arg_inputFile, arg_pathToData]


#-----------------------------------------------------------
#  MAIN PROGRAM
#-----------------------------------------------------------
if __name__ == "__main__":

    # Read command line
    [inputFile, pathToData] = read_args(sys.argv[1:])

    # Instantiate the PiNN object
    pinn = PINNcl.PINN(inputFile, pathToData)

    # Load the data from VOICEXX
    pinn.load_VOICEXX_resu()

    # Build the index vectors for collocation, initial and boundary conditions
    pinn.buildIndexes()

    # Randomly select a percentage of the indexes for train and test
    pinn.selectRandomIndexesTrainTest()

    # Create TensorFlow variables from the test data
    # evaluated at the randomly selected indexes
    pinn.createTFvblesTest()

    # Training Loop
    tstart_training = time.time()
    for epoch in range(pinn.no_epochs):
        # Shuffle train indexes so that the batches are different for each epoch
        pinn.randomizeTrainIndexes()
        for batch in range(pinn.no_batches):
            pinn.batch = batch
            pinn.epoch = epoch
            pinn.getBatch()
            # Create TensorFlow variables from the train data evaluated
            # at the randomly selected indexes
            pinn.createTFvblesTrain()
            pinn.train()
            print("...training loss over batches = {}".format(pinn.loss_batch.numpy()), end='\r')

        if epoch%pinn.freq_save_data == 0:
            pinn.lossEpoch()
            pinn.test()
            t3 = time.time()
            print("\n", end='\r')
            print("Time = {} -> Epoch {}/{}:".format(t3-tstart_training,
                                                     epoch+1, pinn.no_epochs))
            print("             Train loss = {}, Test loss = {}".\
                  format(pinn.train_loss[-1].numpy(), pinn.test_loss[-1].numpy()))
            print("             L1 = {}, L2 = {}, L3 = {}, L4 = {}, L5 = {}, L6 = {}".\
                  format(pinn.L1[-1].numpy(), pinn.L2[-1].numpy(),
                         pinn.L3[-1].numpy(), pinn.L4[-1].numpy(),
                         pinn.L5[-1].numpy(), pinn.L6[-1].numpy()))
            if pinn.verify_simu == 1:
                pinn.append_Coeff()
                print("             Coefficients:[{},{}]".\
                      format(pinn.Coeff[0].numpy(), pinn.Coeff[1].numpy()))

        if epoch%pinn.freq_save_pred == 0:
            # Plot the losses and the prediction
            pinn.plot_losses()
            pinn.plot_prediction()
            if pinn.verify_simu == 1:
                pinn.plot_Coeff()
            # Save the model every 100 epochs
            pinn.save_model()

    tend_training = time.time()

    # Save the training and test losses
    pinn.save_training()
    tend_save_training = time.time()

    print("Time for training: {}".format(tend_training-tstart_training))
    print("Time for saving the training losses: {}".\
          format(tend_save_training-tend_training))
    print("\n")

    # Save final model
    pinn.save_model()
    print("Model Saved")

    # Plot the losses and the prediction
    pinn.plot_losses()
    pinn.plot_prediction()
