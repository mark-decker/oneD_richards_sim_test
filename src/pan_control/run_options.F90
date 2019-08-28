module pan_run_options_mod

   use constants_mod

   use abstract_run_options_mod

   implicit none

   type, extends(run_info) :: pan_opts
      real(kind=uk) :: iteration_tolerance 
      integer :: max_iterations
      contains
         procedure :: read_namelist => read_pan_namelist

   end type pan_opts

   contains

        
     subroutine read_pan_namelist(self,nml_file)
        class(pan_opts), intent(inout) :: self
        character(len=*) :: nml_file

        character(len=265) :: my_iomsg
        integer :: my_iostat

        logical :: file_found
        character(len=*), parameter :: sub_name='read_namelist'
        integer :: myunit

        namelist/pan_isoil/isoil_one,isoil_two,isoil_three,isoil_four,&
                            isbg_one,isbg_two,isbg_three,isbg_four,&
                            ntime,dtime,dz,q_infil,&
                            q_drain,initial_wtd,nlevel,nhorz,iteration_tolerance,max_iterations

        self%namelist_file = nml_file

        inquire( file=trim(self%namelist_file), exist=file_found )
        if (.not.file_found) then
           write(*,*) 'The input file '//trim(self%namelist_file)//' cannot be found.'
           stop
        end if


        open(newunit=myunit,file=trim(self%namelist_file),status='old')
        read(myunit,nml=pan_isoil)
        close(myunit)

        self%ntime       =  ntime       
        self%dtime       =  dtime          
        self%dz          =  dz                 
        self%q_infil     =  q_infil       
        self%q_drain     =  q_drain                    
        self%initial_wtd =  initial_wtd 
        self%nlevel      =  nlevel     
        self%nhorz       =  nhorz
        self%isoils(:)       =  (/isoil_one,isoil_two,isoil_three,isoil_four/)
        self%isoils_beg(:)   =  (/isbg_one,isbg_two,isbg_three,isbg_four/)
        self%iteration_tolerance =   iteration_tolerance
        self%max_iterations = max_iterations

        allocate( self%layer_thickness(self%nlevel) ) 
        self%layer_thickness(:) = self%dz

     end subroutine read_pan_namelist
end module pan_run_options_mod
