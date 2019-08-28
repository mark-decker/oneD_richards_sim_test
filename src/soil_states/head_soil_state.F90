module head_soil_state_mod

   use constants_mod
   use matrix_utils
   use BrooksCorey_parameters_mod
   use ClappHorn_parameters_mod
   use abstract_soil_state_mod

   type, extends(ab_soil_state) :: h_soil_state
      !prognostic variable
      real(kind=uk), dimension(:), allocatable :: smp,dsmpds,&
                                                  dsedh,&
                                                  dhkdh0,dhkdh1,dhkdh2,dhkdh
      real(kind=uk), dimension(:), allocatable :: se_eq,smp_tot
      real(kind=uk), allocatable, dimension(:) :: se, h, hk, zq, liq_vol, pt,bt,&
                                                  hk_star
       contains
         procedure, pass :: alloc => alloc_h_soil_state
         procedure, pass :: h_layer_hydraulic_conductivity
         procedure, pass :: h_swc_func
         procedure, pass :: h_inv_swc_func
         procedure, pass :: lin_eq_potential
         procedure, pass :: h_se_to_volliq
         procedure, pass :: h_volliq_to_se

         generic :: swc  => h_swc_func
         generic :: inv_swc => h_inv_swc_func
         generic :: hydraulic_conductivity => h_layer_hydraulic_conductivity
         generic :: equilibrium_potential => lin_eq_potential
         generic :: se_to_volliq => h_se_to_volliq
         generic :: volliq_to_se  => h_volliq_to_se

         procedure :: update_diags => h_update_diags, &
                      head_from_pt

   end type h_soil_state




   type, extends(h_soil_state) :: iter_h_soil_state
      real(kind=uk), dimension(:), allocatable :: smp_a,h_a,pt_a,se_a  !half time level vars
      contains
         !procedure, pass :: iter_h_layer_hydraulic_conductivity
         procedure :: h_update_iteration_diags
         procedure :: specific_wat_cap   !nneded at half level
         !generic :: hydraulic_conductivity => iter_h_layer_hydraulic_conductivity
   end type iter_h_soil_state


   interface clone
      module procedure :: clone_iter_h_soil_state
   end interface

   contains

     function clone_iter_h_soil_state(state)
         type(h_soil_state), intent(in) :: state
         type(iter_h_soil_state) :: clone_iter_h_soil_state

         call clone_iter_h_soil_state%alloc(state%npts)        
                      
         clone_iter_h_soil_state%dsmpds =  state%dsmpds              
         clone_iter_h_soil_state%bt =      state%bt                  
         clone_iter_h_soil_state%pt =      state%pt                  
         clone_iter_h_soil_state%h =       state%h                   
         clone_iter_h_soil_state%smp =     state%smp                 
         clone_iter_h_soil_state%dsedh =   state%dsedh               
         clone_iter_h_soil_state%dhkdh0 =  state%dhkdh0              
         clone_iter_h_soil_state%dhkdh1 =  state%dhkdh1              
         clone_iter_h_soil_state%dhkdh2 =  state%dhkdh2              
         clone_iter_h_soil_state%dhkdh =   state%dhkdh               
         clone_iter_h_soil_state%se =      state%se                  
         clone_iter_h_soil_state%liq_vol = state%liq_vol  
         clone_iter_h_soil_state%zq =      state%zq                  
         clone_iter_h_soil_state%se_eq =   state%se_eq               
         clone_iter_h_soil_state%hk =      state%hk                  
         clone_iter_h_soil_state%hk_star = state%hk_star  
         clone_iter_h_soil_state%smp_tot = state%smp_tot   

         clone_iter_h_soil_state%pt_a =   state%pt       
         clone_iter_h_soil_state%h_a =    state%h                   
         clone_iter_h_soil_state%se_a =   state%se                 
         clone_iter_h_soil_state%smp_a =  state%smp                 

      end function

     subroutine alloc_h_soil_state(self,npts) 
        class(h_soil_state), intent(inout) :: self
        integer, intent(in) :: npts

        real(kind=uk), dimension(:),allocatable :: zero_vec

        self%npts = npts

        if (allocated(self%dsmpds))      deallocate(self%dsmpds)      
        if (allocated(self%bt))      deallocate(self%bt)      
        if (allocated(self%pt))      deallocate(self%pt)      
        if (allocated(self%h))      deallocate(self%h)      
        if (allocated(self%smp))    deallocate(self%smp)     
        if (allocated(self%dsedh))  deallocate(self%dsedh)      
        if (allocated(self%dhkdh0)) deallocate(self%dhkdh0)     
        if (allocated(self%dhkdh1)) deallocate(self%dhkdh1)     
        if (allocated(self%dhkdh2)) deallocate(self%dhkdh2)    
        if (allocated(self%dhkdh)) deallocate(self%dhkdh)    
        if (allocated(self%se))     deallocate(self%se)        
        if (allocated(self%liq_vol)) deallocate(self%liq_vol)
        if (allocated(self%zq))     deallocate(self%zq)      
        if (allocated(self%se_eq))  deallocate(self%se_eq)    
        if (allocated(self%hk))     deallocate(self%hk)          
        if (allocated(self%hk_star))     deallocate(self%hk_star)          
        if (allocated(self%smp_tot)) deallocate(self%smp_tot)

        allocate(self%dsmpds(self%npts)) 
        allocate(self%bt(self%npts)) 
        allocate(self%pt(self%npts)) 
        allocate(self%h(self%npts)) 
        allocate(self%smp(self%npts)) 
        allocate(self%dsedh(self%npts))
        allocate(self%dhkdh0(self%npts)) 
        allocate(self%dhkdh1(self%npts)) 
        allocate(self%dhkdh2(self%npts)) 
        allocate(self%dhkdh(self%npts)) 
        allocate(self%se(self%npts)) 
        allocate(self%liq_vol(self%npts))   
        allocate(self%zq(self%npts)) 
        allocate(self%se_eq(self%npts)) 
        allocate(self%hk(self%npts)) 
        allocate(self%hk_star(self%npts)) 
        allocate(self%smp_tot(self%npts)) 

        select type (self)
        type is (iter_h_soil_state)

        !states at time n+1/2
        if (allocated(self%pt_a))      deallocate(self%pt_a)      
        if (allocated(self%se_a))      deallocate(self%se_a)      
        if (allocated(self%smp_a))      deallocate(self%smp_a)      
        if (allocated(self%h_a))      deallocate(self%h_a)      

        allocate(self%se_a(self%npts)) 
        allocate(self%h_a(self%npts)) 
        allocate(self%pt_a(self%npts)) 
        allocate(self%smp_a(self%npts)) 

        end select

     end subroutine alloc_h_soil_state

      subroutine h_swc_func(self,param)
         class(h_soil_state), intent(inout) :: self
         class(ab_soil_param)  , intent(inout) :: param
         integer :: i,j,k
         real(kind=uk), dimension(self%npts) :: tmp_vec

         select type (param)
         type is (bc_param)

            self%smp(:) = min(param%smp_sat(:), self%h(:) )
            tmp_vec(:) = min(param%smp_sat(:), max(smp_min, self%smp(:) ) )
            self%se(:) = ( param%smp_sat(:)/tmp_vec(:) )**param%lam(:)
            self%dsedh(:) = ( param%smp_sat(:)/tmp_vec(:) )**(param%lam(:)-1._uk) / tmp_vec(:)
            where (self%h .gt. param%smp_sat(:)) self%dsedh(:) = 0._uk

         type is (ch_param)

            self%smp(:) = min(param%smp_sat(:), self%h(:) )
            tmp_vec(:) = min(param%smp_sat(:), max(smp_min, self%smp(:) ) )
            self%se(:) = ( param%smp_sat(:)/tmp_vec(:) )**param%lam(:)
            self%dsedh(:) = ( param%smp_sat(:)/tmp_vec(:) )**(param%lam(:)-1._uk) / tmp_vec(:)
            where (self%h .gt. param%smp_sat(:)) self%dsedh(:) = 0._uk

         end select

     end subroutine h_swc_func

     subroutine h_inv_swc_func(self,param)
         class(h_soil_state), intent(inout) :: self
         class(ab_soil_param)  , intent(inout) :: param 

         integer :: i,j,k
         real(kind=uk), dimension(self%npts) :: tmp_vec

         select type (param)
         type is (bc_param)

            tmp_vec(:)  = max(0.01_uk,min(1._uk,self%se(:) ) )
            self%smp(:) = param%smp_sat(:) * tmp_vec(:)**(-param%bsw(:))

         type is (ch_param)

            tmp_vec(:)  = max(0.01_uk,min(1._uk,self%se(:) ) )
            self%smp(:) = param%smp_sat(:) * tmp_vec(:)**(-param%bsw(:))

         end select
     end subroutine h_inv_swc_func

      subroutine lin_eq_potential(self,param)
         class(h_soil_state),intent(inout) :: self
         class(ab_soil_param), intent(inout) ::  param
         !real(kind=uk), dimension(:), intent(in) :: smp_wtd
         real(kind=uk) :: temp,temp0,tempi,temp_smp
         integer :: i,j,k,inxt
         real(kind=uk), dimension(self%npts) :: vol_eq,pos_smp_sat
         real(kind=uk) :: voleq1

         select type (param)
         type is (bc_param)

            self%zq(:)    = min(param%smp_sat(:), param%smp_sat(:) + self%wtd - param%z(:) )
            self%se_eq(:) = param%liq_avl(:) * (self%zq(:) / param%smp_sat(:)) ** param%lam(:) + param%liq_res(:)

         type is (ch_param)

            self%zq(:)    = min(param%smp_sat(:), param%smp_sat(:) + self%wtd - param%z(:) )
            self%se_eq(:) = param%liq_sat(:) * (self%zq(:) / param%smp_sat(:)) ** param%lam(:) 

         end select

      end subroutine lin_eq_potential

      subroutine h_layer_hydraulic_conductivity(self,param)
         class(h_soil_state), intent(inout) :: self
         class(ab_soil_param), intent(inout) :: param

         integer :: j

         select type (self)
         type is (h_soil_state)

            select type(param)
            class is (bc_param)

            do j=1,self%npts
                self%hk(j) = param%hk_sat(j)*(self%h(j)/param%smp_sat(j))**(-2.0 - 3.0/param%bsw(j))
                self%hk_star(j) = self%hk(j) * ((1.0+self%bt(j)*self%h(j))**2 ) 
            end do

            class is (ch_param)

            do j=1,self%npts
                self%hk(j) = param%hk_sat(j)*(self%h(j)/param%smp_sat(j))**(-2.0 - 3.0/param%bsw(j))
                self%hk_star(j) = self%hk(j) * ((1.0+self%bt(j)*self%h(j))**2 ) 
            end do

            end select

         type is (iter_h_soil_state)

            select type(param)
            class is (bc_param)
   
            do j=1,self%npts
                self%hk(j) = param%hk_sat(j)*(self%h_a(j)/param%smp_sat(j))**(-2.0 - 3.0/param%bsw(j))
                self%hk_star(j) = self%hk(j) * ((1.0+self%bt(j)*self%h_a(j))**2 ) 
            end do
   
            class is (ch_param)
   
            do j=1,self%npts
                self%hk(j) = param%hk_sat(j)*(self%h_a(j)/param%smp_sat(j))**(-2.0 - 3.0/param%bsw(j))
                self%hk_star(j) = self%hk(j) * ((1.0+self%bt(j)*self%h_a(j))**2 ) 
            end do
   
            end select

         end select

      end subroutine h_layer_hydraulic_conductivity

!      subroutine iter_h_layer_hydraulic_conductivity(self,param)
!         class(iter_h_soil_state), intent(inout) :: self
!         class(ab_soil_param), intent(inout) :: param
!
!         integer :: j
!
!         select type(param)
!         class is (bc_param)
!
!         do j=1,self%npts
!             self%hk(j) = param%hk_sat(j)*(self%h_a(j)/param%smp_sat(j))**(-2.0 - 3.0/param%bsw(j))
!             self%hk_star(j) = self%hk(j) * ((1.0+self%bt(j)*self%h_a(j))**2 ) 
!         end do
!
!         class is (ch_param)
!
!         do j=1,self%npts
!             self%hk(j) = param%hk_sat(j)*(self%h_a(j)/param%smp_sat(j))**(-2.0 - 3.0/param%bsw(j))
!             self%hk_star(j) = self%hk(j) * ((1.0+self%bt(j)*self%h_a(j))**2 ) 
!         end do
!
!         end select
!      end subroutine iter_h_layer_hydraulic_conductivity



      subroutine h_volliq_to_se(self,param)
         class(h_soil_state), intent(inout) :: self
         class(ab_soil_param)  , intent(inout) :: param
         integer :: i,j,k
         real(kind=uk), dimension(self%npts) :: tmp_vec

         select type (param)
         type is (bc_param)

            self%se(:) = (self%liq_vol(:) - param%liq_res(:)) / param%liq_avl(:)

         type is (ch_param)

            self%se(:) = (self%liq_vol(:) ) / param%liq_sat(:)

         end select

     end subroutine h_volliq_to_se

      subroutine h_se_to_volliq(self,param)
         class(h_soil_state), intent(inout) :: self
         class(ab_soil_param)  , intent(inout) :: param
         integer :: i,j,k
         real(kind=uk), dimension(self%npts) :: tmp_vec


         select type (param)
         type is (bc_param)

            self%liq_vol(:) = param%liq_avl(:) * self%se(:) + param%liq_res(:)

         type is (ch_param)

            self%liq_vol(:) = param%liq_sat(:) * self%se(:) 

         end select

     end subroutine h_se_to_volliq

     function head_from_pt(self,param)
        class(h_soil_state), intent(inout) :: self
        class(ab_soil_param), intent(inout) :: param
        real(kind=uk), dimension(self%npts) :: head_from_pt
     
        integer :: k
    
        select type (param)
        class is (bc_param) 
        do k=1,self%npts
         if (self%pt(k) < 0.0) then
             if (self%pt(k) == (1/self%bt(k))) then
                self%pt(k) = self%pt(k) + 0.000000001
             end if
             head_from_pt(k)= self%pt(k)/(1-self%bt(k)*self%pt(k))
          else
             head_from_pt(k) = self%pt(k)
          end if
        end do
     
        class is (ch_param)
        do k=1,self%npts
         if (self%pt(k) < 0.0) then
             if (self%pt(k) == (1/self%bt(k))) then
                self%pt(k) = self%pt(k) + 0.000000001
             end if
             head_from_pt(k)= self%pt(k)/(1-self%bt(k)*self%pt(k))
          else
             head_from_pt(k) = self%pt(k)
          end if
        end do

        end select 
     end function
     
     subroutine specific_wat_cap(self,param,htmp)
        class(iter_h_soil_state), intent(inout) :: self
        class(ab_soil_param), intent(inout) :: param
        real(kind=uk), dimension(self%npts),optional :: htmp
     
        real(kind=uk), dimension(self%npts) :: h_tmp
        integer :: k
     
        if (present(htmp)) then
           h_tmp = htmp
        else
           h_tmp = self%h
        end if
     
        select type (param)
        class is (bc_param)
     
           do j=1,self%npts
              if (h_tmp(j) .le. param%smp_sat(j)) then
                 self%dsmpds(j) = (self%smp_a(j)/param%smp_sat(j))**(-1.0-1.0/param%bsw(j)) &
                                       * (-1.0/param%bsw(j)) * (param%liq_avl(j))*&
                                         (1.0 + self%bt(j) * self%h_a(j))**2.0
              else
                 self%dsmpds(j) = 0.0
              end if
           end do
     
        class is (ch_param)
     
           do j=1,self%npts
              if (h_tmp(j) .le. param%smp_sat(j)) then
                 self%dsmpds(j) = (self%smp_a(j)/param%smp_sat(j))**(-1.0-1.0/param%bsw(j)) &
                                       * (-1.0/param%bsw(j)) * (param%liq_sat(j))*&
                                         (1.0 + self%bt(j) * self%h_a(j))**2.0
              else
                 self%dsmpds(j) = 0.0
              end if
           end do
     
        end select
     
     end subroutine specific_wat_cap


      subroutine h_update_diags(self,param)
         class(h_soil_state), intent(inout) :: self
         class(ab_soil_param), intent(inout) :: param
         integer :: i,j,k

         select type (param)
         class is (bc_param)

            self%h   = self%head_from_pt(param)
            self%smp = min(param%smp_sat, self%h)
            call self%swc(param)

         class is (ch_param)

            self%h   = self%head_from_pt(param)
            self%smp = min(param%smp_sat, self%h)
            call self%swc(param)

         end select

      end subroutine h_update_diags

      subroutine h_update_iteration_diags(self,state,param)
         class(iter_h_soil_state), intent(inout) :: self
         class(h_soil_state), intent(inout) :: state
         class(ab_soil_param), intent(inout) :: param
         integer :: i,j,k
     
         select type (param)
         class is (bc_param) 

            self%h_a = half*(state%h + self%h)
            self%smp_a = min(param%smp_sat, self%h_a)
            call self%swc(param)
      
         class is (ch_param) 

            self%h_a = half*(state%h + self%h)
            self%smp_a = min(param%smp_sat, self%h_a)
            call self%swc(param)
      
         end select

      end subroutine h_update_iteration_diags
      

end module head_soil_state_mod

