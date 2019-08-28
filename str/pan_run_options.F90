module pan_run_options_mod

   use constants_mod

   use abstract_run_options_mod
   use option_vars_mod

   implicit none

   type, extends(run_info) :: pan_opts
      real(kind=uk) :: iteration_tolerance
      real(kind=uk), pointer :: bt
      integer :: max_iterations
      contains
         procedure :: read_namelist => read_pan_namelist
         procedure :: construct => construct_pan_opts

   end type pan_opts

   contains

     subroutine read_pan_namelist(self)
        class(pan_opts), intent(inout) :: self

        character(len=265) :: my_iomsg
        integer :: my_iostat

        logical :: file_found
        character(len=*), parameter :: sub_name='read_namelist'
        integer :: myunit, fourvec(4), i

        namelist/pan_isoil/isoil_one,isoil_two,isoil_three,isoil_four,&
                            isbg_one,isbg_two,isbg_three,isbg_four,&
                            ntime,dtime,dz,q_infil,&
                            q_drain,initial_wtd,nlevel,nhorz,iteration_tolerance,max_iterations,&
                            bt

        inquire( file=trim(self%namelist_file), exist=file_found )
        if (.not.file_found) then
           write(*,*) 'The input file '//trim(self%namelist_file)//' cannot be found.'
           stop
        end if


        open(newunit=myunit,file=trim(self%namelist_file),status='old')
        read(myunit,nml=pan_isoil)
        close(myunit)
   
        self%nml_done = .true.

     end subroutine read_pan_namelist


      subroutine construct_pan_opts(self)
         class(pan_opts), intent(inout) :: self

         integer, dimension(:), pointer, save :: saved_local_isoils, &
                                                         saved_local_isoils_beg
         real(kind=uk), dimension(:), pointer, save :: saved_local_thickness

         integer :: i

        self%ntime       =  ntime       
        self%dtime       =  dtime          
        self%dz          =  dz                 
        self%q_infil     =  q_infil       
        self%q_drain     =  q_drain                    
        self%initial_wtd =  initial_wtd 
        self%nlevel      =  nlevel     
        self%nhorz       =  nhorz
        self%bt       =  bt
        self%iteration_tolerance =   iteration_tolerance
        self%max_iterations = max_iterations

        if (associated(self%isoils)) nullify(self%isoils)
        if (associated(self%isoils_beg)) nullify(self%isoils_beg)
        if (associated(self%layer_thickness)) nullify(self%layer_thickness)

        allocate(saved_local_isoils(max(4,nhorz)))
        allocate(saved_local_isoils_beg(max(4,nhorz)))
        allocate(saved_local_thickness(nhorz))

        do i=1,nhorz
           select case(i)
           case(1)
              saved_local_isoils(i) = isoil_one
              saved_local_isoils_beg(i) = isbg_one
           case(2)
              saved_local_isoils(i) = isoil_two
              saved_local_isoils_beg(i) = isbg_two
           case(3)
              saved_local_isoils(i) = isoil_three
              saved_local_isoils_beg(i) = isbg_three
           case(4)
              saved_local_isoils(i) = isoil_four
              saved_local_isoils_beg(i) = isbg_four
           case default
              write(*,*) 'currently only four horizons are coded for.  Add more'
           end select
        end do

        saved_local_thickness(:) = spread(dz,1,nhorz)

        self%isoils => saved_local_isoils
        self%isoils_beg => saved_local_isoils_beg
        self%layer_thickness => saved_local_thickness

     end subroutine construct_pan_opts

end module pan_run_options_mod
