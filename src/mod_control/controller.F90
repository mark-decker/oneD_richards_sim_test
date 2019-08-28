module mod_control_mod

   use constants_mod
   use matrix_utils
   use BrooksCorey_parameters_mod
   use ClappHorn_parameters_mod
   use mass_soil_state_mod
   use modrich_run_options_mod
   use abstract_sim_controller_mod

   implicit none

   !defines procedures soil_water_sim class must have
   !defines how the solution is found, typically created
   !for classes that extend ab_soil_state
   type, extends(ab_sim_controller) :: mod_controller
      class(mod_opts), allocatable :: opts 
      class(mass_soil_state), allocatable :: st
      class(ab_soil_param), allocatable :: sp
      class(tridiag_matrix), allocatable :: tri 

      integer :: iteraction_count, ntime
      real(kind=uk) :: dtime
      real(kind=uk), dimension(:), allocatable ::  qflx,               &
                                                  dqin_ds0,dqin_ds1,&
                                                  dqout_ds1, dqout_ds2,&
                                                  dqds1_grid, dse_dt0,&
                                                  hk_face, dhkdf_ds0,dhkf_ds1,dhkf_ds2
      integer, allocatable, dimension(:) :: sat_layers
      contains
         procedure :: setup => modrich_setup
         procedure :: alloc                     => modrich_alloc
         procedure :: update                   => modrich_update
         procedure :: initialize  => modrich_initialize
         procedure :: get_soln                 => modrich_get_soln
         procedure :: step_time => modrich_step_time
         
         procedure :: determine_timestep
         procedure :: check_correct_saturation_transitions ,&
                      modrich_internode_conductivity

   end type mod_controller


  contains 
      subroutine modrich_setup(self,npts,nml_file,opts,param,state,tri)
         class(mod_controller), intent(inout) :: self
         character(len=*), intent(in) :: nml_file
         integer, intent(in) :: npts

         class(mod_opts) , intent(inout) :: opts
         class(mass_soil_state)       , intent(inout) :: state
         class(ab_soil_param)    , intent(inout) :: param
         class(tridiag_matrix)   , intent(inout) :: tri 

         self%npts = npts

         call self%alloc(opts,param,state,tri)

         call self%opts%read_namelist(nml_file)

         call self%sp%set_params(self%opts)

      end subroutine modrich_setup
      

     subroutine modrich_alloc(self,opts,param,state,tri)
        class(mod_controller), intent(inout) :: self
        class(mod_opts) , intent(inout) :: opts
        class(mass_soil_state)       , intent(inout) :: state
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
         allocate(self%tri , source= tri)
        !call self%tri%alloc(self%npts)

        allocate(self%qflx(0:self%npts) )
        allocate(self%dqin_ds0(self%npts) )
        allocate(self%dqin_ds1(self%npts) )
        allocate(self%dqout_ds1(self%npts) ) 
        allocate(self%dqout_ds2(self%npts) ) 
        allocate(self%dqds1_grid(self%npts) )
        allocate(self%dse_dt0(self%npts)   )
        allocate(self%hk_face(0:self%npts)   )
        allocate(self%dhkdf_ds0(self%npts)  )
        allocate(self%dhkf_ds1(self%npts)   ) 
        allocate(self%dhkf_ds2(self%npts)  )
        allocate(self%sat_layers(self%npts) )

     end subroutine modrich_alloc


      subroutine modrich_update(self)
         class(mod_controller), intent(inout) :: self

         integer :: i,j,k,ilev

         ilev = self%npts
         associate(state => self%st,&
         param => self%sp,&
         opts => self%opts,&
         tri => self%tri)


         !update the dependent states (smp, hk)
         call state%inv_swc(param)

         call state%hydraulic_conductivity(param)

         !update conductivity on faces
         call self%modrich_internode_conductivity()

         ilev=self%npts

         select type(param)
         type is (bc_param)

            do i=1,ilev-1
               state%smp_tot(i) = ( state%smp(i+1) - state%smp(i) ) - &
                                 ( state%zq(i+1) - state%zq(i) )
            end do
   
            state%smp_tot(ilev) = 0._uk!self%smp(ilev) - param%zq(ilev)
   
            self%qflx(0) = opts%q_infil
            do i=1,ilev-1
               self%qflx(i) = -self%hk_face(i) * state%smp_tot(i)/(param%z(i+1)-param%z(i))
            end do
            self%qflx(ilev) = opts%q_drain
   
            self%dse_dt0 = self%qflx(0:ilev-1) - self%qflx(1:ilev)
   
            self%dqin_ds0(1) = 0._uk
            self%dqin_ds1(1) = 0._uk
            do i=2,ilev
               self%dqin_ds0(i) = -(-self%hk_face(i-1)*state%dsmpds(i-1) + &
                                     state%smp_tot(i-1)*self%dhkf_ds1(i-1))/(param%z(i)-param%z(i-1))
               self%dqin_ds1(i) = -( self%hk_face(i-1)*state%dsmpds(i)   + &
                                     state%smp_tot(i-1)*self%dhkf_ds2(i-1))/(param%z(i)-param%z(i-1))
            end do
   
            self%dqout_ds1(ilev) = 0._uk
            self%dqout_ds2(ilev) = 0._uk
            do i=1,ilev-1
               self%dqout_ds1(i) = -(-self%hk_face(i)*state%dsmpds(i)   + &
                                      state%smp_tot(i)*self%dhkf_ds1(i))/(param%z(i+1)-param%z(i))
               self%dqout_ds2(i) = -( self%hk_face(i)*state%dsmpds(i+1) + &
                                      state%smp_tot(i)*self%dhkf_ds2(i))/(param%z(i+1)-param%z(i))
            end do

        
         end select

         end associate
      end subroutine modrich_update

      function modrich_get_soln(self)
         class(mod_controller), intent(inout)   :: self
         real(kind=uk), dimension(self%npts) :: modrich_get_soln

         real(kind=uk), dimension(self%npts) :: liq_eff,liq_r_eff

         integer :: ilev

         ilev = self%npts
         associate(state => self%st,&
         param => self%sp,&
         opts => self%opts,&
         tri => self%tri)

         select type (param)
         type is (bc_param) 

            call tri%set(-self%dqin_ds0(:), &   !a
                         param%liq_avl(:)*param%dz(:)/self%dtime - self%dqin_ds1(:)+self%dqout_ds1(:), & !b
                         self%dqout_ds2(:), &  !c
                         self%qflx(0:ilev-1) - self%qflx(1:ilev) ) !r

         type is (ch_param)


            call tri%set(-self%dqin_ds0(:), &   !a
                         param%liq_sat(:)*param%dz(:)/self%dtime - self%dqin_ds1(:)+self%dqout_ds1(:), & !b
                         self%dqout_ds2(:), &  !c
                         self%qflx(0:ilev-1) - self%qflx(1:ilev) ) !r

         end select

         modrich_get_soln = tri%get_soln()

         end associate         

         return
 
      end function


      subroutine modrich_initialize(self)
         class(mod_controller), intent(inout)   :: self

         
         integer :: ilev

         ilev = self%npts

         associate(state => self%st,&
         param => self%sp,&
         opts => self%opts,&
         tri => self%tri)

         state%wtd = opts%initial_wtd

         call state%equilibrium_potential(param)

         state%smp(:) = state%zq(:)
         state%se(:) = state%se_eq(:)

         call state%se_to_volliq(param)

         end associate 
 
      end subroutine modrich_initialize

      subroutine modrich_step_time(self)
         class(mod_controller), intent(inout)   :: self

         real(kind=uk), dimension(self%npts) :: dwat
         
         integer :: ilev

         ilev = self%npts
         associate(state => self%st,&
         param => self%sp,&
         opts => self%opts,&
         tri => self%tri)

         dwat = self%get_soln()

         state%se(:) = state%se(:) + dwat(:)

         end associate

      end subroutine modrich_step_time


      !procedures unique to the modrich solution

      subroutine check_correct_saturation_transitions(self)
         class(mod_controller), intent(inout) :: self

         integer :: i,j,k
         real(kind=uk), dimension(self%npts*2) :: dt_test
         real(kind=uk) :: mass_needed, delta_flux(self%npts), delta_liq_min(self%npts)
         integer ::  dt_level

         
         
         
         
         integer :: ilev

         ilev = self%npts
         associate(state => self%st,  &
         param => self%sp,&
         opts => self%opts,&
         tri => self%tri)

         ilev=self%npts

         do i=1,ilev
            if (state%se(i) .ge. self%opts%rel_sat_min) then 
               self%sat_layers(i) = 1
            else
               self%sat_layers(i) = 0
           end if
        end do

        !if a layer will saturate  with a timestep of 
        !less than dtime_min seconds, force it to saturation
        !by 1) take from direction of explicit fluxes, unless self causes newly unsat layer
        select type (param)
        type is (bc_param)

           delta_liq_min(:) = self%opts%dtime_min/(param%dz(:)*param%liq_avl(:))&
                                     * (self%qflx(0:ilev-1)-self%qflx(1:ilev))
           do i=1,ilev
              mass_needed = zero

              if ( (state%se(i) + delta_liq_min(i) .gt. self%opts%rel_sat_min) .and.&
                      (state%se(i) .lt. self%opts%rel_sat_min) ) then

                  self%sat_layers(i) = 1  !mark to be treated as a sat layer

                  mass_needed = ((1._uk-state%se(i))*param%liq_avl(i)+param%liq_res(i))*param%dz(i)
                  state%se(i) = 1._uk
                  state%liq_vol(i) = param%liq_sat(i)

                  if ( (state%se(max(1,i-1)) .lt. self%opts%rel_sat_min) .and.&
                      (state%se(min(i+1,ilev)) .ge. self%opts%rel_sat_min) ) then
                     state%liq_vol(i-1) = state%liq_vol(i-1) - mass_needed / param%dz(i-1)
                     state%se(i-1) = (state%liq_vol(i-1)-param%liq_res(i-1))/param%liq_avl(i-1)
                     mass_needed = zero
                  elseif ( (state%se(min(i+1,ilev)) .lt. self%opts%rel_sat_min)&
                                     .and. (state%se(max(1,i-1)) .ge. self%opts%rel_sat_min) ) then
                     state%liq_vol(i+1) = state%liq_vol(i+1) - mass_needed / param%dz(i+1)
                     state%se(i+1) = (state%liq_vol(i+1)-param%liq_res(i+1))/param%liq_avl(i+1)
                     mass_needed = zero
                  else  !take from first non sat from bottom
                     do j=ilev,1,-1
                        if (i .ne. j .and.&
                             state%se(j) .lt. self%opts%rel_sat_min .and.&
                                mass_needed .gt. eps) then
                           state%liq_vol(j) = state%liq_vol(j) - mass_needed / param%dz(j)
                           state%se(j) = (state%liq_vol(j)-param%liq_res(j))/param%liq_avl(j)
                           mass_needed = zero
                        end if
                     end do
                  end if
               end if
               if (mass_needed .gt. eps) then
                  j=1
                  state%liq_vol(j) = state%liq_vol(j) - mass_needed / param%dz(j)
                  state%se(j) = (state%liq_vol(j)-param%liq_res(j))/param%liq_avl(j)
                  mass_needed = zero
               end if 
           end do

           !now mark layer to be treated as unsat if will become unsat during timestep
           do i=1,ilev
              if ( (state%se(i) + delta_liq_min(i) .lt. self%opts%rel_sat_min) .and.&
                                          (state%se(i) .gt. self%opts%rel_sat_min) ) then
                 self%sat_layers(i) = 0  !mark to be solved as unsat
              end if
           end do

        type is (ch_param)

           delta_liq_min(:) = self%opts%dtime_min/(param%dz(:)*param%liq_sat(:))&
                                     * (self%qflx(0:ilev-1)-self%qflx(1:ilev))
           do i=1,ilev
              mass_needed = zero

              if ( (state%se(i) + delta_liq_min(i) .gt. self%opts%rel_sat_min) .and.&
                      (state%se(i) .lt. self%opts%rel_sat_min) ) then

                  self%sat_layers(i) = 1  !mark to be treated as a sat layer

                  mass_needed = ((1._uk-state%se(i))*param%liq_sat(i))*param%dz(i)
                  state%se(i) = 1._uk
                  state%liq_vol(i) = param%liq_sat(i)

                  if ( (state%se(max(1,i-1)) .lt. self%opts%rel_sat_min) .and.&
                      (state%se(min(i+1,ilev)) .ge. self%opts%rel_sat_min) ) then
                     state%liq_vol(i-1) = state%liq_vol(i-1) - mass_needed / param%dz(i-1)
                     state%se(i-1) = (state%liq_vol(i-1))/param%liq_sat(i-1)
                     mass_needed = zero
                  elseif ( (state%se(min(i+1,ilev)) .lt. self%opts%rel_sat_min)&
                                     .and. (state%se(max(1,i-1)) .ge. self%opts%rel_sat_min) ) then
                     state%liq_vol(i+1) = state%liq_vol(i+1) - mass_needed / param%dz(i+1)
                     state%se(i+1) = (state%liq_vol(i+1))/param%liq_sat(i+1)
                     mass_needed = zero
                  else  !take from first non sat from bottom
                     do j=ilev,1,-1
                        if (i .ne. j .and.&
                             state%se(j) .lt. self%opts%rel_sat_min .and.&
                                mass_needed .gt. eps) then
                           state%liq_vol(j) = state%liq_vol(j) - mass_needed / param%dz(j)
                           state%se(j) = (state%liq_vol(j))/param%liq_sat(j)
                           mass_needed = zero
                        end if
                     end do
                  end if
               end if
               if (mass_needed .gt. eps) then
                  j=1
                  state%liq_vol(j) = state%liq_vol(j) - mass_needed / param%dz(j)
                  state%se(j) = (state%liq_vol(j))/param%liq_sat(j)
                  mass_needed = zero
               end if 
           end do

           !now mark layer to be treated as unsat if will become unsat during timestep
           do i=1,ilev
              if ( (state%se(i) + delta_liq_min(i) .lt. self%opts%rel_sat_min) .and.&
                                          (state%se(i) .gt. self%opts%rel_sat_min) ) then
                 self%sat_layers(i) = 0  !mark to be solved as unsat
              end if
           end do


        end select

       end associate 
           
      end subroutine check_correct_saturation_transitions

      subroutine determine_timestep(self)
         class(mod_controller), intent(inout) :: self
         integer :: i,j,k,ilev
         real(kind=uk), dimension(self%npts*2) :: dt_test
         real(kind=uk) :: mass_needed, delta_flux(self%npts), delta_liq_min(self%npts)
         integer ::  dt_level

         
         
         
         

         ilev = self%npts
         associate(state => self%st,  &
         param => self%sp,&
         opts => self%opts,&
         tri => self%tri)


         !est dt
         dt_test(:) = self%opts%dtime_max
         
        delta_flux(:) = self%qflx(0:ilev-1)-self%qflx(1:ilev)
        where (abs(delta_flux(:)) .lt. half_eps) delta_flux(:) = sign(delta_flux,half_eps)

         select type(param)
         type is (bc_param)

            do i=1,ilev
               dt_test(i) = param%dz(i)*(param%liq_sat(i)-param%liq_res(i))*&
                                    (1._uk-state%se(i)) / delta_flux(i)
            end do
            do i=1,ilev
               dt_test(i+ilev) = self%opts%delta_liq_max*param%dz(i) / delta_flux(i)
            end do

            self%dtime = minval(dt_test(:),dim=1)
            dt_level   = minloc(dt_test,dim=1)

         end select         
         end associate
      end subroutine determine_timestep

      subroutine modrich_internode_conductivity(self)
         class(mod_controller), intent(inout) :: self

         integer :: i,j,ilev,k

         ilev = self%npts

         associate(state => self%st, &
         param => self%sp,&
         opts => self%opts,&
         tri => self%tri)

         select type(param)
         type is (bc_param)

         do i=1,self%npts

            j = min(i+1,self%npts)

             !simple
             self%hk_face(i) = (param%dz(i)*state%hk(i) + &
                                param%dz(j)*state%hk(j))/(param%dz(i)+param%dz(i))

             self%dhkf_ds1(i)  = param%dz(i)*state%dhkds(i)/(param%dz(i)+param%dz(j))
             self%dhkf_ds2(i)  = param%dz(j)*state%dhkds(j)/(param%dz(i)+param%dz(j))

         end do

         end select 

         end associate
      end subroutine modrich_internode_conductivity


end module mod_control_mod

