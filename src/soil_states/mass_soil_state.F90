module mass_soil_state_mod

   use constants_mod
   use matrix_utils
   use BrooksCorey_parameters_mod
   use ClappHorn_parameters_mod
   use abstract_soil_state_mod

   !two concrete soil state classes
   !soil_state uses mass as prognostic while
   !h_soil_state uses soil pressure as prognostic
   !todo: trans_h_soil_state: uses transformed pressure as prognostic

   !Note that state procedures need to be written for each parameter
   !class they will run with (here bc and ch)

   type, extends(ab_soil_state) :: mass_soil_state
      !prognostic variable
      real(kind=uk), dimension(:), allocatable :: &
                                                  dsmpds,dhkds,&
                                                  dhkds0,dhkds1,dhkds2,smp
      real(kind=uk), dimension(:), allocatable :: se_eq,smp_tot
      real(kind=uk), allocatable, dimension(:) :: se, h, hk, zq, liq_vol
       contains
         procedure, pass :: alloc => alloc_soil_state
         procedure, pass :: swc_func
         procedure, pass :: inv_swc_func
         procedure, pass :: layer_hydraulic_conductivity
         procedure, pass :: nonlin_eq_potential
         procedure, pass :: state_se_to_volliq 
         procedure, pass :: state_volliq_to_se 

         generic :: swc => swc_func
         generic :: inv_swc => inv_swc_func
         generic :: hydraulic_conductivity => layer_hydraulic_conductivity
         generic :: equilibrium_potential => nonlin_eq_potential
         generic :: se_to_volliq => state_se_to_volliq
         generic :: volliq_to_se => state_volliq_to_se
   end type mass_soil_state


   contains

     subroutine alloc_soil_state(self,npts) 
        class(mass_soil_state), intent(inout) :: self
        integer, intent(in) :: npts

        self%npts = npts

        if (allocated(self%liq_vol )) deallocate(self%liq_vol )  
        if (allocated(self%se )     ) deallocate(self%se )         
        if (allocated(self%dsmpds ) ) deallocate(self%dsmpds )   
        if (allocated(self%dhkds )  ) deallocate(self%dhkds )    
        if (allocated(self%dhkds0 ) ) deallocate(self%dhkds0 )   
        if (allocated(self%dhkds1 ) ) deallocate(self%dhkds1 )    
        if (allocated(self%dhkds2 ) ) deallocate(self%dhkds2 )   
        if (allocated(self%zq )     ) deallocate(self%zq )       
        if (allocated(self%se_eq )  ) deallocate(self%se_eq )    
        if (allocated(self%hk )     ) deallocate(self%hk )       
        if (allocated(self%smp_tot )) deallocate(self%smp_tot )   
        if (allocated(self%smp )) deallocate(self%smp )   


        allocate(self%liq_vol(self%npts) )
        allocate(self%se(self%npts) )
        allocate(self%dsmpds(self%npts) )
        allocate(self%dhkds(self%npts) )
        allocate(self%dhkds0(self%npts) )
        allocate(self%dhkds1(self%npts) )
        allocate(self%dhkds2(self%npts) )
        allocate(self%zq(self%npts) )
        allocate(self%se_eq(self%npts) )
        allocate(self%hk(self%npts) )
        allocate(self%smp_tot(self%npts) )
        allocate(self%smp(self%npts) )

     end subroutine alloc_soil_state

      subroutine swc_func(self,param)
         class(mass_soil_state), intent(inout) :: self
         class(ab_soil_param)  , intent(inout) :: param
         integer :: i,j,k
         real(kind=uk), dimension(self%npts) :: tmp_vec

         select type (param)
         type is (bc_param)

            tmp_vec(:) = min(param%smp_sat(:), max(smp_min, self%smp(:) ) )
            self%se(:) = ( param%smp_sat(:)/tmp_vec(:) )**param%lam(:)

         type is (ch_param)

           tmp_vec(:) = min(param%smp_sat(:), max(smp_min, self%smp(:) ) )
           self%se(:) = ( param%smp_sat(:)/tmp_vec(:) )**param%lam(:)
         
         end select
     end subroutine swc_func
     subroutine inv_swc_func(self,param)
         class(mass_soil_state), intent(inout) :: self
         class(ab_soil_param)  , intent(inout) :: param 

         integer :: i,j,k
         real(kind=uk), dimension(self%npts) :: tmp_vec

         select type (param)
         type is (bc_param)

            tmp_vec(:)  = max(0.01_uk,min(1._uk,self%se(:) ) )
            self%smp(:) = param%smp_sat(:) * tmp_vec(:)**(-param%bsw(:))
            self%dsmpds(:) = -param%bsw(:) * param%liq_avl(:) * self%smp(:) / tmp_vec(:)

         type is (ch_param)

            tmp_vec(:)  = max(0.01_uk,min(1._uk,self%se(:) ) )
            self%smp(:) = param%smp_sat(:) * tmp_vec(:)**(-param%bsw(:))
            self%dsmpds(:) = -param%bsw(:) * param%liq_sat(:) * self%smp(:) / tmp_vec(:)

         end select

     end subroutine inv_swc_func

      subroutine layer_hydraulic_conductivity(self,param)
         class(mass_soil_state), intent(inout) :: self
         class(ab_soil_param)  , intent(inout) :: param
         integer :: i
         real(kind=uk),dimension(self%npts) :: se_tmp
         real(kind=uk) :: tmp_Avar,tmp_Bvar

         select type (param)
         type is (bc_param)

            se_tmp = max(0.01_uk  , self%se(:) )
            self%hk    = param%hk_sat * se_tmp**(param%hk_exp(:))
            self%dhkds = param%hk_exp(:) * self%hk(:) / se_tmp

         type is (ch_param)

            se_tmp = max(0.01_uk  , self%se(:) )
            self%hk    = param%hk_sat * se_tmp**(param%hk_exp(:))
            self%dhkds = param%hk_exp(:) * self%hk(:) / se_tmp

         end select

      end subroutine layer_hydraulic_conductivity


      subroutine nonlin_eq_potential(self,param)
         class(mass_soil_state),intent(inout) :: self
         class(ab_soil_param), intent(inout) ::  param
         !real(kind=uk), dimension(:), intent(in) :: smp_wtd
         real(kind=uk) :: temp,temp0,tempi,temp_smp
         integer :: i,j,k,inxt
         real(kind=uk), dimension(self%npts) :: vol_eq,pos_smp_sat
         real(kind=uk) :: voleq1

         select type (param)
         type is (bc_param)

            self%se_eq(:) = 1._uk
            self%zq(:)    = param%smp_sat(:)

            do i=1,self%npts
               if (self%wtd .le. param%zi(i-1)) then!elseif (self%wtd .le. param%zi(i)) then

                  self%se_eq(i) = 1._uk

               elseif (self%wtd .lt. param%zi(i)) then

                  tempi  = 1._uk
                  temp0 = ( 1._uk- (self%wtd-param%zi(i-1))/param%smp_sat(i) )**(1._uk-param%lam(i))
                  temp  = (param%smp_sat(i)*param%liq_avl(i)*(tempi-temp0)/ &
                                ( (1._uk - param%lam(i))*(self%wtd-param%zi(i-1)) ) ) + param%liq_res(i)
                  temp = (temp*(self%wtd-param%zi(i-1)) + &
                              param%liq_sat(i)*(param%zi(i)-self%wtd))/(param%zi(i)-param%zi(i-1))
                  self%se_eq(i) =max(0.001_uk, min(1._uk, (temp - param%liq_res(i))/param%liq_avl(i) ) )
               else
                  tempi = ( 1._uk -(self%wtd-param%zi(i))/param%smp_sat(i) ) ** (1._uk - param%lam(i))
                  temp0 = ( 1._uk -(self%wtd-param%zi(i-1))/param%smp_sat(i) )**(1._uk-param%lam(i))
                  temp  = (param%smp_sat(i)*param%liq_avl(i)/(1-param%lam(i))*&
                                      (tempi-temp0))/(param%zi(i)-param%zi(i-1))+&
                                        param%liq_res(i)
                  self%se_eq(i) =max(0.001_uk, min(1._uk, (temp - param%liq_res(i))/param%liq_avl(i) ) )
              end if
              self%zq(i) = param%smp_sat(i) * self%se_eq(i)**(-param%bsw(i))

            end do

         type is (ch_param)

            self%se_eq(:) = 1._uk
            self%zq(:)    = param%smp_sat(:)

            self%se_eq(:) = 1.0
            self%zq(:)    = param%smp_sat(:)

            do i=1,self%npts

               if (self%wtd .le. param%zi(i-1)) then!elseif (self%wtd .le. param%zi(i)) then

                  self%se_eq(i) = 1._uk

               elseif (self%wtd .lt. param%zi(i)) then

                  tempi  = 1._uk
                  temp0 = ( 1._uk- (self%wtd-param%zi(i-1))/param%smp_sat(i) )**(1._uk-param%lam(i))
                  temp  = (param%smp_sat(i)*param%liq_sat(i)*(tempi-temp0)/ &
                                ( (1._uk - param%lam(i))*(self%wtd-param%zi(i-1)) ) ) 
                  temp = (temp*(self%wtd-param%zi(i-1)) + &
                               param%liq_sat(i)*(param%zi(i)-self%wtd))/(param%zi(i)-param%zi(i-1))
                  self%se_eq(i) =max(0.001_uk, min(1._uk, (temp  ) ) )
               else
                  tempi = ( 1._uk -(self%wtd-param%zi(i))/param%smp_sat(i) ) ** (1._uk - param%lam(i))
                  temp0 = ( 1._uk -(self%wtd-param%zi(i-1))/param%smp_sat(i) )**(1._uk-param%lam(i))
                  temp  = (param%smp_sat(i)*param%liq_sat(i)/(1-param%lam(i))*&
                                      (tempi-temp0))/(param%zi(i)-param%zi(i-1))
                  self%se_eq(i) =max(0.001_uk, min(1._uk, (temp  ) ) )
              end if
              self%zq(i) = param%smp_sat(i) * self%se_eq(i)**(-param%bsw(i))

            end do

         end select

      end subroutine nonlin_eq_potential

      subroutine state_volliq_to_se(self,param)
         class(mass_soil_state), intent(inout) :: self
         class(ab_soil_param)  , intent(inout) :: param
         integer :: i,j,k
         real(kind=uk), dimension(self%npts) :: tmp_vec

         select type (param)
         type is (bc_param)

            self%se(:) = (self%liq_vol(:) - param%liq_res(:)) / param%liq_avl(:)

         type is (ch_param)

            self%se(:) = (self%liq_vol(:)) / param%liq_sat(:)

         end select

     end subroutine state_volliq_to_se

      subroutine state_se_to_volliq(self,param)
         class(mass_soil_state), intent(inout) :: self
         class(ab_soil_param)  , intent(inout) :: param
         integer :: i,j,k
         real(kind=uk), dimension(self%npts) :: tmp_vec

         select type (param)
         type is (bc_param)

            self%liq_vol(:) = param%liq_avl(:) * self%se(:) + param%liq_res(:)

         type is (ch_param)

            self%liq_vol(:) = param%liq_sat(:) * self%se(:)

         end select

     end subroutine state_se_to_volliq



end module mass_soil_state_mod
