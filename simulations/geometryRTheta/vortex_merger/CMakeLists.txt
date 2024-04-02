

add_executable(vortex_merger
     vortex_merger.cpp
)
target_compile_features(vortex_merger PUBLIC cxx_std_17)
target_link_libraries(vortex_merger
    PUBLIC       
        DDC::DDC
	    DDC::PDI_Wrapper
	    paraconf::paraconf
	    PDI::pdi
	        
	    Eigen3::Eigen
	        
	    sll::SLL
	    gslx::paraconfpp
	       
	    gslx::poisson_RTheta
	    gslx::geometry_RTheta
	    gslx::interpolation_2D_rp
	    gslx::advection_rp
	    gslx::time_integration_rp
	    gslx::quadrature
	    gslx::utils
	    gslx::vortex_merger_init_eq
		gslx::simulation_utils
)

