module BrooksCorey_parameters_mod

   use constants_mod, only : smp_min, eps, zero , one, two, three, ten, half , half_eps , close_val, nhorz_max
   use matrix_utils
   use abstract_soil_parameter_mod

   !currently two options: Brooks-Corey and Cosby using the same hydrualic cond params
   !Brooks Corey
   type, extends(ab_soil_param) ::  bc_param
      integer :: kind=uk
      character(len=:), allocatable :: my_name

      real(kind=uk), allocatable, dimension(:) :: liq_sat,smp_sat,lam,bsw,  &
                                              hk_exp,swc_exp,invswc_exp,&
                                              hk_sat, liq_res,liq_avl

      contains
         procedure :: set_params => set_bc_params
         procedure :: get_isoil => get_bc_isoil
         procedure :: alloc       => alloc_bc
         procedure :: clone      => clone_bc
   end type bc_param

   interface bc_param
      module procedure make_bc_param
   end interface

  contains 

     !Brooks-corey subroutines are all first
     subroutine set_bc_params(self,rinfo)
        class(bc_param), intent(inout) :: self
        class(run_info), intent(in) :: rinfo
        integer :: i,j,k,end_lev

        call self%get_isoil(rinfo)

        do i=1,self%npts
           k = self%isoils(i)
           k = max(1, min(11, k))

           self%hk_sat(i)  = rinfo%soil_props(k,3)
           self%smp_sat(i) = rinfo%soil_props(k,4)
           self%bsw(i)     = rinfo%soil_props(k,5)
           self%liq_sat(i) = rinfo%soil_props(k,1)
           self%liq_res(i) = rinfo%soil_props(k,2)
        end do


        self%liq_avl(:) = self%liq_sat(:) - self%liq_res(:)
        self%lam(:) = 1._uk / self%bsw(:)

        self%hk_exp(:) = 2._uk*self%bsw(:) + 3._uk

        self%zi(0) = 0.
        self%dz(:) = rinfo%dz
        do i=1,self%npts
           self%zi(i) = self%zi(i-1) + self%dz(i)
           self%z(i)  = self%zi(i-1) + 0.5*self%dz(i)
        end do


     end subroutine set_bc_params

     !contructors
     function make_bc_param(npts)
        type(bc_param) :: make_bc_param
        integer, intent(in) :: npts
        call alloc_bc(make_bc_param,npts)
     end function
      
     subroutine clone_bc(self,new) 
        class(bc_param), intent(inout) :: self
        class(ab_soil_param), allocatable, intent(out) :: new

        allocate(new, source=self)

     end subroutine

     subroutine alloc_bc(self,npts)
        class(bc_param), intent(inout) :: self
        integer, intent(in)    :: npts
      
        if (.not.allocated(self%my_name)) then
           self%my_name='soil_param'
        end if

        self%npts=npts
        allocate(self%liq_avl(npts))
        allocate(self%liq_sat(npts))
        allocate(self%liq_res(npts))
        allocate(self%hk_sat(npts))
        allocate(self%hk_exp(npts))
        allocate(self%smp_sat(npts))
        allocate(self%lam(npts))
        allocate(self%bsw(npts))
        allocate(self%z(npts))
        allocate(self%dz(npts))
        allocate(self%zi(0:npts))
        allocate(self%isoils(npts))
        self%isoils(:)  = -1
        self%liq_avl(:) = zero
        self%liq_sat(:) = zero
        self%liq_res(:) = zero
        self%hk_sat(:)  = zero
        self%hk_exp(:)  = zero
        self%smp_sat(:) = zero
        self%lam(:)     = zero
        self%bsw(:)     = zero
        self%z(:)       = zero
        self%dz(:)      = zero
        self%zi(:)      = zero

     end subroutine alloc_bc

     subroutine get_bc_isoil(self,rinfo)
        class(bc_param),intent(inout) :: self 
        class(run_info), intent(in) :: rinfo

        integer :: nhrz,npts,i,ii,j,jj,k,kk,nhrz_max

        nhrz_max = size(rinfo%isoils,dim=1)
        npts = self%npts

        nhrz=0
        do i=1,nhrz_max
           if (rinfo%isoils(i) .ge. 1 .and. rinfo%isoils(i) .le. 11 .and. &
               rinfo%isoils_beg(i) .ge. 1 .and. rinfo%isoils_beg(i) .le. npts) then 

              nhrz = nhrz + 1
           end if
        end do

        i=1
        do k=1,nhrz
          if (k .eq. nhrz) then 
             ii = npts 
          else
             ii = i + rinfo%isoils_beg(k+1)-1
          end if
          do j=i,ii
             self%isoils(j) = rinfo%isoils(k)
          end do
        end do

     end subroutine get_bc_isoil

end module BrooksCorey_parameters_mod

