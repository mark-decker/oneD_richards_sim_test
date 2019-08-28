module abstract_soil_parameter_mod

   use constants_mod, only : smp_min, eps, zero , one, two, three, ten, half , half_eps , close_val
   use matrix_utils
   use abstract_run_options_mod

   type oned_grid
      integer :: npts
      real(kind=uk), dimension(:), allocatable :: z, dz, zi
   end type

   !defines procedures soil_param class must have
   !holds all parameters including for hyd cond and SWF
   type, extends(oned_grid), abstract :: ab_soil_param
      integer, allocatable :: isoils(:)
      contains
         procedure(param_int_iface)     , pass, deferred :: alloc
         procedure(param_info_iface)     , pass, deferred :: set_params
         procedure(new_param_iface)   , pass, deferred :: clone
   end type ab_soil_param

   !ab_soil_param interfaces
   abstract interface 
      subroutine new_param_iface(self,new)
         import ab_soil_param
         class(ab_soil_param), intent(inout) :: self
         class(ab_soil_param),allocatable, intent(out) :: new
      end subroutine
   end interface

   abstract interface 
      subroutine param_info_iface(self,rinfo)
         import ab_soil_param
         import run_info
         class(ab_soil_param), intent(inout) :: self
         class(run_info), intent(in) :: rinfo
      end subroutine param_info_iface
   end interface



   abstract interface 
      subroutine param_int_iface(self,npts)
         import ab_soil_param
         class(ab_soil_param), intent(inout) :: self
         integer, intent(in) :: npts
      end subroutine param_int_iface
   end interface


   abstract interface 
      subroutine param_iface(self)
         import ab_soil_param
         class(ab_soil_param), intent(inout) :: self
      end subroutine param_iface
   end interface


end module abstract_soil_parameter_mod

