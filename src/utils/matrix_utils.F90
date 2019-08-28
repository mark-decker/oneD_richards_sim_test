module matrix_utils

   use constants_mod

   type ::  tridiag_matrix
      integer :: npts
      real(kind=uk), allocatable, dimension(:) :: amx,bmx,cmx,rmx,gam,soln
      contains
         procedure :: get_soln => get_tridiag_soln
         procedure :: setup_solve_tridiag_vectors
         procedure :: error   => tridiag_error
         procedure :: set     => set_tridiag_vectors
         procedure :: alloc    => tridiag_alloc
         procedure :: destroy => tridiag_destroy
         procedure :: clone   => tridiag_clone
   end type tridiag_matrix
   
   interface tridiag_matrix
       procedure make_tridiag
   end interface

   
  contains 

     function get_tridiag_soln(self)
        class(tridiag_matrix) :: self
        real(kind=uk), dimension(self%npts) :: get_tridiag_soln

        real(uk) :: bet
        integer :: j

        bet = self%bmx(1)
        get_tridiag_soln(1) = self%rmx(1) / bet
        do j=2,self%npts

           self%gam(j) = self%cmx(j-1)/bet
           bet = self%bmx(j)-self%amx(j)*self%gam(j)
           get_tridiag_soln(j) = (self%rmx(j)-self%amx(j)*get_tridiag_soln(j-1))/bet
        end do
        get_tridiag_soln(self%npts) = self%rmx(self%npts)/self%bmx(self%npts)
        do j=(self%npts-1),1,-1
           get_tridiag_soln(j) = get_tridiag_soln(j)-self%gam(j+1)*get_tridiag_soln(j+1)
        end do

     end function

     subroutine tridiag_error(self,message)
        class(tridiag_matrix) :: self
        character(len=*), intent(in) :: message

        write(*,*) 'ERROR '//trim(message)
        stop
    end subroutine tridiag_error

     subroutine set_tridiag_vectors(self,a,b,c,r)
        class(tridiag_matrix) :: self
        real(uk), dimension(:), intent(in) :: a,b,c,r

        self%amx(:) = a(:)
        self%bmx(:) = b(:)
        self%cmx(:) = c(:)
        self%rmx(:) = r(:)

     end subroutine set_tridiag_vectors

     subroutine setup_solve_tridiag_vectors(self,a,b,c,r,soln)
        class(tridiag_matrix) :: self
        real(uk), dimension(:), intent(in) :: a,b,c,r
        real(uk), dimension(:), intent(inout) :: soln

        call self%set(a,b,c,r)
     
        soln(:) = self%get_soln()  

     end subroutine setup_solve_tridiag_vectors

!     subroutine solve(self) 
!        class(tridiag_matrix) :: self
!
!        real(uk) :: bet
!        integer :: j
!
!        bet = self%bmx(1)
!        self%soln(1) = self%rmx(1) / bet
!        do j=2,self%npts
!
!           self%gam(j) = self%cmx(j-1)/bet
!           bet = self%bmx(j)-self%amx(j)*self%gam(j)
!           self%soln(j) = (self%rmx(j)-self%amx(j)*self%soln(j-1))/bet
!        end do
!        self%soln(self%npts) = self%rmx(self%npts)/self%bmx(self%npts)
!        do j=(self%npts-1),1,-1
!           self%soln(j) = self%soln(j)-self%gam(j+1)*self%soln(j+1)
!        end do
!
!     end subroutine solve

     function make_tridiag(npts)
        type(tridiag_matrix) :: make_tridiag
        integer, intent(in) :: npts

        call make_tridiag%alloc(npts)

     end function

     subroutine tridiag_alloc(self,npts)
        class(tridiag_matrix),intent(out) :: self
        integer, intent(in)   :: npts

        self%npts=npts

        if (.not.allocated(self%amx) )  allocate(self%amx(npts))
        if (.not.allocated(self%bmx) )  allocate(self%bmx(npts))
        if (.not.allocated(self%cmx) )  allocate(self%cmx(npts))
        if (.not.allocated(self%rmx) )  allocate(self%rmx(npts))
        if (.not.allocated(self%gam) )  allocate(self%gam(npts))
        if (.not.allocated(self%soln))  allocate(self%soln(npts))

     end subroutine tridiag_alloc

     subroutine tridiag_destroy(self)
        class(tridiag_matrix), intent(out) :: self
     
        if (allocated(self%amx) ) deallocate(self%amx) 
        if (allocated(self%bmx) ) deallocate(self%bmx) 
        if (allocated(self%cmx) ) deallocate(self%cmx) 
        if (allocated(self%rmx) ) deallocate(self%rmx) 
        if (allocated(self%gam) ) deallocate(self%gam) 
        if (allocated(self%soln)) deallocate(self%soln)

     end subroutine tridiag_destroy

     function tridiag_clone(in_tridiag) result(new_tridiag)
        class(tridiag_matrix), intent(in)  :: in_tridiag
        type(tridiag_matrix)  :: new_tridiag

        allocate( new_tridiag%amx, source=in_tridiag%amx )
        allocate( new_tridiag%bmx, source=in_tridiag%bmx )   
        allocate( new_tridiag%cmx, source=in_tridiag%cmx ) 
        allocate( new_tridiag%rmx, source=in_tridiag%rmx )  
        allocate( new_tridiag%gam, source=in_tridiag%gam ) 
        allocate(new_tridiag%soln, source=in_tridiag%soln)
   
     end function

end module matrix_utils
