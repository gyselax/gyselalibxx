'''
General functions needed for PINN method
'''

# pylint: disable=invalid-name
# pylint: disable=import-error

from contextlib import redirect_stdout
from tensorflow.keras.layers import Dense, Input
from tensorflow.keras.models import Model

#-----------------------------------------------
# Define Tensorflow model
#-----------------------------------------------
def myModel(hidden_dim=(20, 20, 20), input_dim=2, output_dim=1,
            hidden_act_func='relu', output_act_func='relu',
            file_for_summary='summary.dat'):
    '''
    Define Tensorflow model
    '''

    # prepare input placeholder
    num_layers = len(hidden_dim)
    input_shape = (input_dim,)
    inputs = Input(shape=input_shape)
    x = Dense(hidden_dim[0], activation=hidden_act_func,
              input_shape=input_shape, dtype='float32')(inputs)

    for layer_count in range(1, num_layers):
        print("Adding layer: ", layer_count)
        x = Dense(hidden_dim[layer_count], activation=hidden_act_func)(x)

    output = Dense(output_dim, activation=output_act_func, dtype='float32')(x)

    model = Model(inputs, output)
    # print out network description
    with open(file_for_summary, 'w', encoding="utf-8") as f:
        with redirect_stdout(f):
            model.summary()

    return model


#-----------------------------------------------
# Save dataset as dataname in filename
#-----------------------------------------------
def save_data(filename, dataname, dataset):
    '''
    Save dataset as dataname in filename
    '''
    filename[dataname] = dataset
