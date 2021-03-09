#include <vlasov.hpp>
/*
subroutine solve_vlasov_xvx( geom, dt, E, all_f, krook_op )

    type(geometry)                , intent(in)    :: geom
    real(F64)                     , intent(in)    :: dt
    real(F64)      , dimension(0:), intent(in)    :: E
    type(all_fdistribu_t)         , intent(inout) :: all_f
    type(krook_t)                 , intent(inout) :: krook_op
    
    real(F64) :: half_dt, nu
    
    half_dt = 0.5_F64*dt
    nu = nu_sink + nu_source

    if (nu.ne.0) &
      call solve_Krook_operator(krook_op,geom,all_f,half_dt)
    if (.not.two_species) then
       if (diff_coeff.ne.0) &
            call solve_Diffusion_operator(all_f%f_elec)
       call advec1D_x(geom,half_dt,all_f%f_elec)
       call advec1D_v(geom,dt,E,all_f%f_elec)
       call advec1D_x(geom,half_dt,all_f%f_elec)
       if (diff_coeff.ne.0) &
            call solve_Diffusion_operator(all_f%f_elec)
    else
       if (diff_coeff.ne.0) then
          call solve_Diffusion_operator(all_f%f_elec)
          call solve_Diffusion_operator(all_f%f_ion)
       end if
       call advec1D_x(geom,half_dt,all_f%f_elec)
       call advec1D_v(geom,dt,E,all_f%f_elec)
       call advec1D_x(geom,half_dt,all_f%f_elec)
       
       call advec1D_x(geom,half_dt,all_f%f_ion)
       call advec1D_v(geom,dt,E,all_f%f_ion)
       call advec1D_x(geom,half_dt,all_f%f_ion)
       if (diff_coeff.ne.0) then
          call solve_Diffusion_operator(all_f%f_elec)
          call solve_Diffusion_operator(all_f%f_ion)
       end if
    end if
    if (nu.ne.0) &
         call solve_Krook_operator(krook_op,geom,all_f,half_dt)

  end subroutine solve_vlasov_xvx*/



void Vlasov::operator()(const Field2D& cur, Field2D& next, double dt) const
{
	m_advec_x
	
}
