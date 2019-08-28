module pan_control_mod

   use constants_mod
   use matrix_utils
   use BrooksCorey_parameters_mod
   use ClappHorn_parameters_mod

   use head_soil_state_mod
   use pan_run_options_mod
   use abstract_sim_controller_mod

   implicit none

   !concrete class implementing the transformed pressure heat approach
   !with picard iterations from pan

   type, extends(ab_sim_controller) :: pan_controller
      class(pan_opts), allocatable :: opts 
      class(h_soil_state), allocatable :: st
      class(iter_h_soil_state), allocatable :: it
      class(ab_soil_param), allocatable :: sp
      class(tridiag_matrix), allocatable :: tri 

      integer :: iteration_count, ntime
      real(kind=uk) :: dtime
      real(kind=uk), dimension(:), allocatable :: pt, hk_star, hk_facestar,&                 
                                                  hk_face
      real(kind=uk), dimension(:), allocatable :: residuals
      real(kind=uk) :: residual
      contains
         !procedure :: interfacial_conductivity => pan_internode_hk_pt_cor
         procedure :: setup => pan_setup
         procedure :: step_time                => pan_time_step
         !procedure :: update_dep_vars          => pan_update_deps

         procedure :: alloc                     => pan_alloc  
         procedure :: update                   => pan_update
         procedure :: initialize  => pan_initialize
         procedure :: get_soln                 => pan_get_soln !none

         procedure :: pan_hk_pt_correction
         procedure :: continue_iterating
         procedure :: iterate => pan_iterate
         procedure :: pan_calc_fluxes, set_tridiag,pan_internode_conductivity

   end type pan_controller

  contains 

      subroutine pan_setup(self,npts,nml_file,opts,param,state,tri)
         class(pan_controller), intent(inout) :: self
         character(len=*), intent(in) :: nml_file
         integer, intent(in) :: npts

         class(pan_opts) , intent(inout) :: opts
         class(h_soil_state)       , intent(inout) :: state
         class(ab_soil_param)    , intent(inout) :: param
         class(tridiag_matrix)   , intent(inout) :: tri 

         self%npts = npts

         call self%alloc(opts,param,state,tri)

         call self%opts%read_namelist(nml_file)

         call self%sp%set_params(self%opts)

      end subroutine pan_setup
      


      subroutine pan_alloc(self,opts,param,state,tri)
         class(pan_controller), intent(inout) :: self
         class(pan_opts) , intent(inout) :: opts
         class(h_soil_state)       , intent(inout) :: state
         class(ab_soil_param)    , intent(inout) :: param
         class(tridiag_matrix)   , intent(inout) :: tri 
         
         self%npts = opts%nlevel
      
         call param%alloc(self%npts)
         call state%alloc(self%npts)
         call tri%alloc(self%npts)

         !call move_alloc(opts,self%opts)
         !call move_alloc(param,self%sp)
         !call move_alloc(state,self%st)
         !call move_alloc(tri,self%tri)


         allocate(self%opts , source= opts)
         allocate(self%sp , source= param)
         allocate(self%st , source= state)
         allocate(self%it , source=clone(self%st))
         allocate(self%tri , source= tri)
      
         allocate(self%residuals(self%npts) )
         allocate(self%pt(self%npts) )
         allocate(self%hk_star(self%npts) )
         allocate(self%hk_facestar(self%npts) )
         allocate(self%hk_face(self%npts) )
      
      end subroutine pan_alloc
      
      !return true if continue iterations
      function continue_iterating(self)
         class(pan_controller), intent(inout) :: self
         logical :: continue_iterating
      
         integer :: i,j,k
         real(kind=uk) :: max_residual
      
         self%iteration_count = self%iteration_count + 1
      
         self%residual = maxval(abs(self%residuals(:)),dim=1)
      
         if ( (self%iteration_count .gt. self%opts%max_iterations ) .or. &
              (self%residual .le. self%opts%iteration_tolerance) )  then
      
              continue_iterating = .false.
              !leaving iteration loop, zero out counter for next time step
              self%iteration_count = 0
         else
              continue_iterating = .true.
         end if
      
      end function
      
      subroutine pan_iterate(self)
         class(pan_controller), intent(inout) :: self
     
         real(kind=uk), dimension(self%npts) :: residuals
 
         associate(param => self%sp)

         select type (param)
         class is (bc_param)
 
         self%tri%rmx(:) = one

         self%residuals(:) = 1.0e20_uk
         self%iteration_count = 0
         do while (self%continue_iterating())
      
            self%it%h = self%it%head_from_pt(param)  !iteration state
            
            !calc_fluxes does what?  needed for hk and facial hk
            call self%pan_calc_fluxes()
      
            call self%it%specific_wat_cap(param)
      
            call self%set_tridiag()
      
            self%it%pt(:) = self%it%pt(:) + self%tri%get_soln()

            self%residuals(:) = self%tri%rmx(:)
      
         end do
      
         class is (ch_param)
 
         self%tri%rmx(:) = one

         self%iteration_count = 0
         do while (self%continue_iterating())
      
            self%it%h = self%it%head_from_pt(param)  !iteration state
            
            !calc_fluxes does what?  needed for hk and facial hk
            call self%pan_calc_fluxes()
      
            call self%it%specific_wat_cap(param)
      
            call self%set_tridiag()
      
            self%it%pt(:) = self%it%pt(:) + self%tri%get_soln()
      
            self%residuals(:) = self%tri%rmx(:)

         end do

         end select

         end associate
      
      end subroutine pan_iterate
      
      
      subroutine pan_time_step(self)
         class(pan_controller), intent(inout) :: self
     
         associate(param => self%sp)

         select type (param)
         class is (bc_param)
 
         self%st%h = self%st%head_from_pt(param)
      
         call self%st%update_diags(param)
      
         call self%iterate()
      
         self%st%pt = self%it%pt  !save new result in state
      
         call self%it%update_diags(param)
         call self%st%update_diags(param)
      
         class is (ch_param)
 
         self%st%h = self%st%head_from_pt(param)
      
         call self%st%update_diags(param)
      
         call self%iterate()
      
         self%st%pt = self%it%pt  !save new result in state
      
         call self%it%update_diags(param)
         call self%st%update_diags(param)
      
         end select

         end associate

      end subroutine pan_time_step
      
      
      subroutine pan_calc_fluxes(self)
         class(pan_controller), intent(inout) :: self
      
      
         real(kind=uk), dimension(self%npts) :: h_it
         integer :: i,j,k
      
         associate(param => self%sp, &
                   state => self%st, &
                   opts => self%opts, &
                   tri => self%tri, &
                   it_state => self%it)
      
         call it_state%update_diags(param)   !n+1 time level
      
         call it_state%h_update_iteration_diags(state,param)  !it state at half n + 1/2
      
         call it_state%hydraulic_conductivity(param)  !use half time lev _a
      
         call self%pan_internode_conductivity()  !use half time lev _a
      
         end associate
      
      end subroutine pan_calc_fluxes

      subroutine pan_update(self)
         class(pan_controller), intent(inout) :: self

         integer :: i,j,k
      
         associate(param => self%sp, &
                   state => self%st, &
                   opts => self%opts, &
                   tri => self%tri, &
                   it_state => self%it)

         call state%update_diags(param)

         it_state%pt = state%pt
         it_state%h = state%h
         it_state%smp = state%smp
         it_state%se = state%se
         it_state%liq_vol = state%liq_vol

         it_state%pt_a = it_state%pt
         it_state%h_a = it_state%h
         it_state%smp_a = it_state%smp
         it_state%se_a = it_state%se

         end associate

      end subroutine pan_update
      
      subroutine set_tridiag(self)
         class(pan_controller), intent(inout) :: self
      
         real(kind=uk) :: Klo_star,Khi_star
         integer :: j
      
         associate(tri => self%tri, &
                   state => self%st, &
                   if_state => self%it, &
                   param => self%sp, &
                   opts => self%opts)
      
         select type (param)
         class is (bc_param)
 
         do j=1,self%npts
            if (j .gt. 1 .and. j .lt. self%npts) then
      
               Klo_star = self%hk_facestar(j-1)/param%dz(j)
               Khi_star = self%hk_facestar(j)/param%dz(j)
      
               tri%amx(j) = -Klo_star
               tri%cmx(j) = -Khi_star
               tri%bmx(j) = param%dz(j)*state%dsmpds(j)/self%dtime +Klo_star+Khi_star
               tri%rmx(j) = (Khi_star*(state%pt(j+1)-state%pt(j)) & 
                         - Klo_star*(state%pt(j)-state%pt(j-1))) & 
                             - (self%hk_star(j)-self%hk_star(j-1)) -&
                               param%dz(j)*(self%it%liq_vol(j)-state%liq_vol(j))/self%dtime
           !       - (Khi-Klo) -param%dz(j)*(self%it%liq_vol(j)-self%liq_vol(j))/self%dtime
      
            elseif (j == 1) then
               !Klo = self%q_infil
               !Khi_ = face_hk(j)
               Klo_star = self%hk_star(j) / param%dz(j)
               Khi_star = self%hk_facestar(j) / param%dz(j)
      
               !a(j) = 0.0
               !c(j) = -Khi_star
               tri%amx(j) = 0.0
               tri%cmx(j) = -Khi_star
               tri%bmx(j) = param%dz(j)*self%it%dsmpds(j)/self%dtime +Klo_star+Khi_star
               tri%rmx(j) = (Khi_star*(state%pt(j+1)-state%pt(j))) & 
                           - (self%hk_star(j)-self%opts%q_infil) -&
                             param%dz(j)*(self%it%liq_vol(j)-state%liq_vol(j))/self%dtime
                 !- (Khi-Klo) -param%dz(j)*(self%it%liq_vol(j)-self%liq_vol(j))/self%dtime
      
      
            elseif (j == self%npts) then
               !Khi = self%q_drain
               !Klo = face_hk(j-1)
               Klo_star = self%hk_facestar(j-1) / param%dz(j)
               Khi_star = self%hk_star(j) / param%dz(j)
      
               !a(j) =  -Klo_star
               !c(j) = 0.0
               !a(j) = -self%hk_facestar(j-1)/self%dz(j)
               !c(j) = 0.0
      
               tri%amx(j) = -Klo_star
               tri%cmx(j) = 0.0
               tri%bmx(j) = param%dz(j)*self%it%dsmpds(j)/self%dtime +Klo_star+Khi_star
               tri%rmx(j) = (-Klo_star*(state%pt(j)-state%pt(j-1))) & 
                           - (self%opts%q_drain-self%hk_star(j-1)) -&
                             param%dz(j)*(self%it%liq_vol(j)-state%liq_vol(j))/self%dtime
            end if
      
         end do
      
         class is (ch_param)
         do j=1,self%npts
            if (j .gt. 1 .and. j .lt. self%npts) then
      
               Klo_star = self%hk_facestar(j-1)/param%dz(j)
               Khi_star = self%hk_facestar(j)/param%dz(j)
      
               tri%amx(j) = -Klo_star
               tri%cmx(j) = -Khi_star
               tri%bmx(j) = param%dz(j)*self%it%dsmpds(j)/self%dtime +Klo_star+Khi_star
               tri%rmx(j) = (Khi_star*(state%pt(j+1)-state%pt(j)) & 
                         - Klo_star*(state%pt(j)-state%pt(j-1))) & 
                             - (self%hk_star(j)-self%hk_star(j-1)) -&
                               param%dz(j)*(self%it%liq_vol(j)-state%liq_vol(j))/self%dtime
           !       - (Khi-Klo) -param%dz(j)*(self%it%liq_vol(j)-self%liq_vol(j))/self%dtime
      
            elseif (j == 1) then
               !Klo = self%q_infil
               !Khi_ = face_hk(j)
               Klo_star = self%hk_star(j) /param%dz(j)
               Khi_star = self%hk_facestar(j) /param%dz(j)
      
               !a(j) = 0.0
               !c(j) = -Khi_star
               tri%amx(j) = 0.0
               tri%cmx(j) = -Khi_star
               tri%bmx(j) = param%dz(j)*self%it%dsmpds(j)/self%dtime +Klo_star+Khi_star
               tri%rmx(j) = (Khi_star*(state%pt(j+1)-state%pt(j))) & 
                           - (self%hk_star(j)-self%opts%q_infil) -&
                             param%dz(j)*(self%it%liq_vol(j)-state%liq_vol(j))/self%dtime
                 !- (Khi-Klo) -param%dz(j)*(self%it%liq_vol(j)-self%liq_vol(j))/self%dtime
      
      
            elseif (j == self%npts) then
               !Khi = self%q_drain
               !Klo = face_hk(j-1)
               Klo_star = self%hk_facestar(j-1) / param%dz(j)
               Khi_star = self%hk_star(j) / param%dz(j)
      
               !a(j) =  -Klo_star
               !c(j) = 0.0
               !a(j) = -self%hk_facestar(j-1)/self%dz(j)
               !c(j) = 0.0
      
               tri%amx(j) = -Klo_star
               tri%cmx(j) = 0.0
               tri%bmx(j) = param%dz(j)*self%it%dsmpds(j)/self%dtime +Klo_star+Khi_star
               tri%rmx(j) = (-Klo_star*(state%pt(j)-state%pt(j-1))) & 
                           - (self%opts%q_drain-self%hk_star(j-1)) -&
                             param%dz(j)*(self%it%liq_vol(j)-state%liq_vol(j))/self%dtime
            end if
      
         end do
     
         end select 

         end associate
      
      end subroutine set_tridiag

      subroutine pan_internode_conductivity(self)
         class(pan_controller), intent(inout) :: self

         integer :: j

         associate(param => self%sp)

         select type (param)
         class is (bc_param)
 
         !first set the values at nodes
         call self%it%hydraulic_conductivity(param)
         call self%st%hydraulic_conductivity(param)


         do j=1,self%npts-1
            self%hk_facestar(j) = 0.5*(self%it%hk_star(j) + self%it%hk_star(j+1))
            self%hk_face(j) = 0.5*(self%it%hk(j) + self%it%hk(j+1))
         end do

         j=self%npts
         self%hk_facestar(j) = self%it%hk_star(j)
         self%hk_face(j) = self%it%hk(j) 
 
         class is (ch_param)
 
         !first set the values at nodes
         call self%it%hydraulic_conductivity(param)
         call self%st%hydraulic_conductivity(param)


         do j=1,self%npts-1
            self%hk_facestar(j) = 0.5*(self%it%hk_star(j) + self%it%hk_star(j+1))
            self%hk_face(j) = 0.5*(self%it%hk(j) + self%it%hk(j+1))
         end do

         j=self%npts
         self%hk_facestar(j) = self%it%hk_star(j)
         self%hk_face(j) = self%it%hk(j) 

         end select

         end associate
 
      end subroutine         

      function pan_get_soln(self)
         class(pan_controller), intent(inout)   :: self
         real(kind=uk), dimension(self%npts) :: pan_get_soln

         real(kind=uk), dimension(self%npts) :: liq_eff,liq_r_eff

         integer :: ilev

         ilev = self%npts
         associate(state => self%st,&
         it_state => self%it, &
         param => self%sp,&
         opts => self%opts,&
         tri => self%tri)

         select type (param)
         class is (bc_param) 

            self%st%h = self%st%head_from_pt(param)
      
            call self%st%update_diags(self%sp)
      
            call self%iterate()

            pan_get_soln(:) = self%it%pt - self%st%pt      
         
         class is (ch_param)

            self%st%h = self%st%head_from_pt(param)
      
            call self%st%update_diags(self%sp)
      
            call self%iterate()

            pan_get_soln(:) = self%it%pt - self%st%pt      
         
         end select

         end associate         

         return
 
      end function


      subroutine pan_hk_pt_correction(self)
         class(pan_controller), intent(inout) :: self

         integer :: i,ilve

         associate(state => self%st,&
         it_state => self%it, &
         param => self%sp,&
         opts => self%opts,&
         tri => self%tri)

         select type (param)
         class is (bc_param) 

            state%pt = it_state%pt
            state%h = state%head_from_pt(param)
            it_state%h = it_state%head_from_pt(param)

            call state%swc(param)
            call it_state%swc(param)

         class is (ch_param)

            state%pt = it_state%pt

            state%h = state%head_from_pt(param)
            it_state%h = it_state%head_from_pt(param)

            call state%swc(param)
            call it_state%swc(param)

         end select

         end associate

      end subroutine pan_hk_pt_correction


      subroutine pan_initialize(self)
         class(pan_controller), intent(inout)   :: self
         
         integer :: ilev

         ilev = self%npts

         associate(state => self%st,&
         it_state => self%it,&
         param => self%sp,&
         opts => self%opts,&
         tri => self%tri)

         state%wtd = opts%initial_wtd

         call state%equilibrium_potential(param)

         state%h(:) = state%zq(:)
         state%se(:) = state%se_eq(:)

         state%pt = state%h / (1._uk + state%bt * state%h )

         select type (param)
         class is (bc_param)

            state%smp = min(state%h, param%smp_sat)
            state%liq_vol = param%liq_avl * state%se + param%liq_res

            it_state%h = state%h
            it_state%smp = state%smp
            it_state%se = state%se
            it_state%liq_vol = state%liq_vol
            it_state%pt = state%pt

         class is (ch_param)

            state%smp = min(state%h, param%smp_sat)
            state%liq_vol = param%liq_sat * state%se 

            it_state%h = state%h
            it_state%smp = state%smp
            it_state%se = state%se
            it_state%liq_vol = state%liq_vol
            it_state%pt = state%pt

         end select

         end associate 
 
      end subroutine pan_initialize

end module pan_control_mod


