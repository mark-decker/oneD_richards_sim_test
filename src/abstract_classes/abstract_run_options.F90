module abstract_run_options_mod

   use constants_mod

   implicit none


   !temporary variables that comprise the comps of run_info that need to
   !be read from nml.  ugly, but removes need for custom IO
   integer :: isbg_one=1,isbg_two=-1,isbg_three=-1,isbg_four=-1
   integer :: isoil_one=1,isoil_two=-1,isoil_three=-1,isoil_four=-1
   integer :: nhorz=1, nlevel=10, ntime = 100
   integer :: max_iterations = 1000
   real(kind=uk) :: dtime=1800._uk , dz= 5._uk ,q_infil = 0._uk, q_drain = 0._uk, &
                    initial_wtd = 2500._uk, dtime_max = 1800._uk, dtime_min=10._uk,&
                    iteration_tolerance=0.1_uk
   real(kind=uk) :: snd_one=0.3_uk,snd_two=-1._uk,snd_three=-1._uk,snd_four=-1._uk
   real(kind=uk) :: slt_one=0.4_uk,slt_two=-1._uk,slt_three=-1._uk,slt_four=-1._uk
   real(kind=uk) :: cly_one=0.3_uk,cly_two=-1._uk,cly_three=-1._uk,cly_four=-1._uk
   real(kind=uk) :: rel_sat_max = 1._uk + 1.0e-5_uk, rel_sat_min = 1._uk - 1.0e-5_uk
   real(kind=uk) :: delta_liq_max = 0.05_uk, delta_se_max = 0.02_uk


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

      integer, dimension(:), allocatable :: isoils,isoils_beg,layer_thickness

      character(len=:), allocatable :: namelist_file 

       character(len=len('run_info')) :: my_name='run_info'

      contains
         procedure(info_char_iface), pass, deferred :: read_namelist

   end type run_info
   abstract interface 
      subroutine info_char_iface(self,nml_file)
         import run_info
         class(run_info), intent(inout) :: self
         character(len=*) :: nml_file
      end subroutine info_char_iface
   end interface



end module abstract_run_options_mod
