The code is now capable of:
1. Infering the distribution function using some percentage(tested for 10%) of
   distibution function data from simulation.
2. Solving Vlasov equation meaning infering the distribution function using
   some percentage(tested for 40%) of only Electric field data from simulation data. 
3. Verfying the solution of PDE by obtaining the coefficients(same as POPE method)

################ RUNNING THE CODE ###############################################
# To run the code following libraries are needed:   
# 1. Tensorflow
# 2. numpy
# 3. sklearn
# 4. matplotlib
# Command to run the code: 
#  python PINNmain.py -i <YAML_input_PINN> -p <VOICEXX_directory>
#   . <YAML_input_PINN> : File used to give input variables for training
#   . <VOICEXX_directory> : Directory of the VOICEXX results
# To run the code on jean-zay submit the script:
#  sbatch jean-zay-gpu-pinn.slurm <YAML_input_PINN> <VOICEXX_directory>
##################################################################################

#################### 'input_PINN.yaml' #################################################
Following are the parameters used in input file:
----------------------------------------------------------------------------------
restart : False      	#When restart the training from the last point it uses the
	  		saved neural network from previous training
learning_rate : 0.0005	#Learning rate for the optimizer used in the training for
	      		now we are using 'Adam' optimizer
train_perc : 95		#Percentage of points used for training and remaining for
	     		testing
batch_perc : 10		#Percentage of training points used in each batch
use_fun_data : 0	#0: To solve Vlasov equation and using only Efield data
	       		#1: To infer distribution function using fdata
verify_simu : 1		#To find coefficients of PDE, to verify that the PDE is
	      		solved correctly or not.
fun_perc : 10		#Percentage of data-points used in training for distribution
	   		function from simulation (when, use_fun_data = 1,this
			percentage only used for calculating testing error.)
PDE_perc : 99		#Percentage of points used in (x,v,t) for training for
	   		PDE(not used when use_fun_data = 1)
x_t_data_perc : 40	#Percentage of Electric field data used to solve Vlasov
	      		 equation (only valid if use_Efield_data = 0)
v_data_perc : 40	#Percentage of v-domain used for training to solve Vlasov
	      		 equation (only valid if use_Efield_data = 0)
init_perc : 99		#Percentage of points used in Initial Condition in training
pbc_0_perc : 99		#Percentage of points used in PBC Condition(x=0,v,t) in
	     		training
pbc_1_perc : 99		#Percentage of points used in PBC Condition(x=x_final,v,t)
	     		in training
dbc_0_perc : 99		#Percentage of points used in DBC Condition(x,v=0,t) in
	     		training
dbc_1_perc : 99		#Percentage of points used in DBC Condition(x,v=v_final,t)
	     		in training
no_epochs : 1000	#No. of epochs for training
input_dim : 3		#Input dimension in Neural Network (currently 3 for [x,v,t])
output_dim : 1		#Output dimension of Neural Network (currently 1 for f,
	     		distribution function)
layers : [100,100,100,100] #Defining the number of nodes in each hidden layer of
       	 		   neural network, in this case 100 nodes in 4 hidden layers
hidden_act_func : tanh	#Activation function of each hidden layer of Neural Network
output_act_func : None	#Activation function of last layer of Neural Network
path_to_results : PINN_Analysis/	#Path to store the results of PINN analysis
model_file : felec_NN.h5		#Neural Network to be used when restarting
	     				the training (restart = True)
freq_save_data : 10			#Freq of saving the data for loss
freq_save_pred : 100			#Frequency of saving the plots for
	       	 			prediction(must be less than 'freq_save_data')
-----------------------------------------------------------------------------------
###################################################################################
