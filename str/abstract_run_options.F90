module abstract_run_options_mod

   use constants_mod
   use option_vars_mod

   implicit none


   !simulations settings that change from run to tun
   !all one d sims will need this info, not that generic :*)
   type, abstract :: run_info
      real(kind=uk) :: dtime      ,&
                       total_time ,&
                       dz         ,&
                       q_infil    ,&
                       q_drain    ,&
                       initial_wtd
                    
      integer :: ntime            ,&
                 nlevel           ,&
                 nhorz            

        real(kind=uk), dimension(11,5) :: soil_props = transpose(reshape(    &   
                              [0.430_uk,0.109_uk,0.000333_uk,-794.8_uk,4.48_uk,    &   
                               0.437_uk,0.035_uk,0.017000_uk,-205.8_uk,1.81_uk,    &   
                               0.453_uk,0.041_uk,0.007200_uk,-302.2_uk,2.64_uk,    &   
                               0.463_uk,0.027_uk,0.001900_uk,-401.2_uk,3.97_uk,    &   
                               0.501_uk,0.015_uk,0.003670_uk,-508.7_uk,4.27_uk,    &   
                               0.398_uk,0.068_uk,0.001200_uk,-594.1_uk,3.13_uk,    &   
                               0.464_uk,0.075_uk,0.000639_uk,-564.3_uk,4.13_uk,    &   
                               0.471_uk,0.075_uk,0.000417_uk,-703.3_uk,5.65_uk,    &   
                               0.430_uk,0.109_uk,0.000333_uk,-794.8_uk,4.48_uk,    &   
                               0.479_uk,0.056_uk,0.000250_uk,-765.4_uk,6.67_uk,    &   
                               0.475_uk,0.090_uk,0.000167_uk,-856.0_uk,6.06_uk]    &    
                               ,[5,11]))

      integer, dimension(:), pointer :: isoils,isoils_beg

      real(kind=uk), dimension(:), pointer :: layer_thickness

      character(len=:), allocatable :: namelist_file 

      character(len=len('run_info')) :: my_name='run_info'

      logical :: nml_done = .false.

      contains
         procedure(info_char_iface), pass, deferred :: read_namelist

   end type run_info
   abstract interface 
      subroutine info_char_iface(self)
         import run_info
         class(run_info), intent(inout) :: self
      end subroutine info_char_iface
   end interface



end module abstract_run_options_mod
