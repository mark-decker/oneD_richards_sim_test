module abstract_soil_state_mod

   use constants_mod
   !defines procedures soil_state class must have
   !defines the prognostic variables, dep variables

   type, abstract :: ab_soil_state
      integer :: npts
      real(kind=uk) :: wtd
      contains
         procedure(state_int_iface)     , pass, deferred :: alloc

   end type ab_soil_state

   abstract interface 
      subroutine state_int_iface(self,npts)
         import ab_soil_state
         class(ab_soil_state), intent(inout) :: self
         integer, intent(in) :: npts
      end subroutine state_int_iface
   end interface


end module abstract_soil_state_mod
