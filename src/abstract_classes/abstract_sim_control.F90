module abstract_sim_controller_mod

   use abstract_soil_parameter_mod
   use abstract_soil_state_mod
   use matrix_utils
   use abstract_run_options_mod

   implicit none
   
   !defines procedures soil_water_sim class must have
   !defines how the solution is found, typically created
   !for classes that extend ab_soil_state

   type, abstract  ::  ab_sim_controller
      integer :: npts
      contains
         procedure(sim_iface), pass, deferred :: update   !update all diagnostic vars
         procedure(sim_iface), pass, deferred :: initialize  !set IC
         procedure(sim_vec_iface)  , pass, deferred :: get_soln  !find change in state
         !procedure(sim_char_all_iface)  , pass, deferred :: setup  !find change in state
         !procedure(sim_all_iface)  , pass, deferred :: alloc  !find change in state

   end type ab_sim_controller
   abstract interface 
      subroutine sim_all_iface(self,rinfo,param,state,tri)
         import ab_sim_controller
         import ab_soil_param
         import ab_soil_state
         import tridiag_matrix
         import run_info
         class(ab_sim_controller), intent(inout) :: self
         class(run_info), intent(inout) :: rinfo
         class(ab_soil_param), intent(inout) :: param
         class(ab_soil_state), intent(inout) :: state
         class(tridiag_matrix), intent(inout) :: tri
      end subroutine
   end interface


!   abstract interface 
!      subroutine sim_char_all_iface(self,nml_file,rinfo,param,state,tri)
!         import ab_sim_controller
!         import ab_soil_param
!         import ab_soil_state
!         import tridiag_matrix
!         import run_info
!         class(ab_sim_controller), intent(inout) :: self
!         character(len=*), intent(in) :: nml_file
!
!         class(run_info), intent(inout) :: rinfo
!         class(ab_soil_param), intent(inout) :: param
!         class(ab_soil_state), intent(inout) :: state
!         class(tridiag_matrix), intent(inout) :: tri
!      end subroutine
!   end interface


   abstract interface 
      function sim_vec_iface(self)
         import ab_sim_controller
         import uk
         class(ab_sim_controller), intent(inout) :: self
         real(kind=uk), dimension(self%npts) :: sim_vec_iface
      end function
   end interface

   abstract interface 
      subroutine sim_iface(self)
         import ab_sim_controller
         class(ab_sim_controller), intent(inout) :: self
      end subroutine sim_iface
   end interface


   abstract interface 
      subroutine new_sim_iface(self,new)
         import ab_sim_controller
         class(ab_sim_controller), intent(in) :: self
         class(ab_sim_controller), intent(out) :: new
      end subroutine
   end interface


   abstract interface
      subroutine sim_int_iface(self,npts)
         import ab_sim_controller
         class(ab_sim_controller)    , intent(inout) :: self
         integer, intent(in) :: npts
      end subroutine sim_int_iface
   end interface

   abstract interface 
      subroutine sim_sim_iface(self,other)
         import ab_sim_controller
         class(ab_sim_controller), intent(inout) :: self
         class(ab_sim_controller), intent(inout) :: other
      end subroutine sim_sim_iface
   end interface

end module abstract_sim_controller_mod

