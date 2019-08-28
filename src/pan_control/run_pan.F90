program pan_sim

   use constants_mod
   use matrix_utils

   !parameter classes
   use BrooksCorey_parameters_mod, only : bc_param  
   use ClappHorn_parameters_mod, only : ch_param

   !prognostic / diagnostic classes
   use mass_soil_state_mod, only : mass_soil_state
   use head_soil_state_mod, only : h_soil_state, iter_h_soil_state

   !to use modified richards
   !use modrich_run_options_mod, only : mod_opts
   !use mod_control_mod, only : modrich_soil_sim

   !to use transformed pressure head from celia
   use pan_run_options_mod, only : pan_opts
   use pan_control_mod, only : pan_controller

implicit none

! this code uses the new richards solution to solve 1-d soil water flow through variably saturated, heterogenious soil
!this code only solves it one time then outputs the necessary data 
!  last edit 10-24-2007
! by mark decker AKA Admeral Awesome

!need to define all of the variables
!integer, parameter :: pan_ctrl%npts=10          !number of soil layers
!integer, parameter  :: ntime= 96         !number of time steps
!real(uk) ::  dtime = 1800.0                !time step (s)
!!  set the boundary conditions
!real(uk) :: qinfl = 1.0/3600.0/24.0
!real(uk) :: wtd = 4999.99999
!timing vars
integer time_array_0(8), time_array_1(8)
!bunch of temp vars to store data in calce
!bunch of integers for looping
integer :: i,j,k,n,&
            jp,&  !j previous j-1
            jn    !j next j+1
integer, dimension(8) :: timeinfob     !stores the time info at the beginning     
integer, dimension(8) :: timeinfoe     !stores the time info at the end
real(uk) :: massb                    !beginning mass balance
real(uk) :: masse                    !ending mass balance
real(uk) :: massin                    !total mass into the domain
real(uk) :: masserror                       !error in the mass balance

integer, allocatable :: isoil_het(:)

integer :: isoil = 1
type(pan_controller)  :: pan_ctrl
type(pan_opts)  :: my_opts
type(h_soil_state)   :: my_state
type(bc_param)     :: my_param
type(tridiag_matrix) :: my_trid

real(kind=uk), allocatable, dimension(:) :: dwat,layer_dz
character(len=1024) :: def_nml = 'pan.run.nml'

call my_opts%read_namelist(def_nml)

call my_param%alloc(my_opts%nlevel)
call my_state%alloc(my_opts%nlevel)
call my_trid%alloc(my_opts%nlevel)

allocate(dwat(my_opts%nlevel))
allocate(layer_dz(my_opts%nlevel))

my_state%wtd = my_opts%initial_wtd

layer_dz(:) = my_opts%dz

!create the sim object using above settings
call pan_ctrl%alloc(my_opts,my_param, my_state, my_trid)

call pan_ctrl%sp%set_params(pan_ctrl%opts)

pan_ctrl%st%wtd = my_opts%initial_wtd

!sewt up the intial conditions
call pan_ctrl%initialize()

!desired api below
!create the sim object that reads nml and makes vars
!call pan_ctrl%setup(my_param, def_nml)
!
!sewt up the intial conditions
!call pan_ctrl%initialize()

!character(len=350) :: filename_isoil
open(11,file='MR_DTlo_3.2_liq.dat',status='replace',access='sequential',action='write')
open(21,file='MR_DTlo_3.2_ic.dat',status='replace',access='sequential',action='write')


massb = sum( pan_ctrl%st%liq_vol*pan_ctrl%sp%dz(:), dim=1)

do j=1,pan_ctrl%npts
   write(21,*) pan_ctrl%st%liq_vol(j)
end do


call date_and_time(values=time_array_0)
write(*,*) (time_array_0(j),j=4,8)


do n=1,pan_ctrl%opts%ntime

   call pan_ctrl%update()

   call pan_ctrl%step_time()

end do      ! end time loop


masse = sum( pan_ctrl%st%liq_vol*pan_ctrl%sp%dz(:), dim=1)

masserror = masse - massb - pan_ctrl%opts%q_infil * pan_ctrl%opts%dtime*real(pan_ctrl%opts%ntime,kind=uk)
 
 call date_and_time(values=time_array_1)
close(11)
close(21)

write(*,*) 'mass b',massb
write(*,*) 'mass e',masse
write(*,*) 'mass eerr',masserror
write(*,*) 'mass err per dt ',masserror/real(pan_ctrl%opts%ntime,kind=uk)

    
end program pan_sim


