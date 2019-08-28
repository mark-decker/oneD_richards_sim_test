module modrich_run_options_mod

   use constants_mod

   use abstract_run_options_mod

   implicit none

   !two types, one for modrich settings other for pan
   !read soil type using one of 11 types
   type, extends(run_info) :: mod_opts
      real(kind=uk) :: &
                       dtime_max  ,&
                       dtime_min  

      real(kind=uk) :: delta_se_max,&
                       rel_sat_max, &
                       rel_sat_min, &
                       delta_liq_max

      contains
         procedure :: read_namelist => read_modrich_namelist

   end type mod_opts


   !read soil texture (mineral) and calc params
   type, extends(mod_opts) :: mod_ped_opts
      real(kind=uk), dimension(nhorz_max) :: fsnd,fslt,fcly
      contains
         procedure :: read_namelist => read_mod_ped_namelist
   end type

   contains

        
     subroutine read_modrich_namelist(self,nml_file)
        class(mod_opts), intent(inout) :: self
        character(len=*) :: nml_file

        character(len=265) :: my_iomsg
        integer :: my_iostat

        logical :: file_found
        character(len=*), parameter :: sub_name='read_namelist'

        integer :: myunit

        namelist/modrich_isoil/isoil_one,isoil_two,isoil_three,isoil_four,&
                             isbg_one,isbg_two,isbg_three,isbg_four,&
                             dtime_max,dtime_min,ntime,dtime,&
                             dz,q_infil,q_drain,initial_wtd,&
                             nlevel,nhorz

        self%namelist_file = nml_file

        inquire( file=trim(self%namelist_file), exist=file_found )
        if (.not.file_found) then
           write(*,*) 'The input file '//trim(self%namelist_file)//' cannot be found.'
           stop
        end if


        !the nml to read and the data depends on type of self
        open(newunit=myunit,file=trim(self%namelist_file),status='old')
        read(myunit,nml=modrich_isoil)
        close(myunit)
   
        self%dtime_max   =  dtime_max              
        self%dtime_min   =  dtime_min
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
        self%delta_liq_max =   delta_liq_max 
        self%delta_se_max  =   delta_se_max  
        self%rel_sat_max   =   rel_sat_max   
        self%rel_sat_min   =   rel_sat_min   

        allocate( self%layer_thickness(self%nlevel) ) 
        self%layer_thickness(:) = self%dz
     end subroutine read_modrich_namelist

        
     subroutine read_mod_ped_namelist(self,nml_file)
        class(mod_ped_opts), intent(inout) :: self
        character(len=*) :: nml_file

        character(len=265) :: my_iomsg
        integer :: my_iostat

        logical :: file_found
        character(len=*), parameter :: sub_name='read_namelist'


        integer :: myunit

        namelist/modrich_pedt/ isbg_one,isbg_two,isbg_three,isbg_four,&
                             snd_one,snd_two,snd_three,snd_four,&
                             slt_one,slt_two,slt_three,slt_four,&
                             cly_one,cly_two,cly_three,cly_four,&
                             dtime_max,dtime_min,ntime,dtime,&
                             dz,q_infil,q_drain,initial_wtd,&
                             nlevel,nhorz,rel_sat_max,rel_sat_min

        self%namelist_file = nml_file

        inquire( file=trim(self%namelist_file), exist=file_found )
        if (.not.file_found) then
           write(*,*) 'The input file '//trim(self%namelist_file)//' cannot be found.'
           stop
        end if

        open(newunit=myunit,file=trim(self%namelist_file),status='old')
        read(myunit,nml=modrich_pedt,iostat=my_iostat,iomsg=my_iomsg)
        close(myunit)

        self%dtime_max   =  dtime_max              
        self%dtime_min   =  dtime_min
        self%ntime       =  ntime       
        self%dtime       =  dtime          
        self%dz          =  dz                 
        self%q_infil     =  q_infil       
        self%q_drain     =  q_drain                    
        self%initial_wtd =  initial_wtd 
        self%nlevel      =  nlevel     
        self%nhorz       =  nhorz
        self%delta_liq_max =   delta_liq_max 
        self%delta_se_max  =   delta_se_max  
        self%rel_sat_max   =   rel_sat_max   
        self%rel_sat_min   =   rel_sat_min   
        self%isoils(:)       =  (/-1,-1,-1,-1/)
        self%isoils_beg(:)   =  (/isbg_one,isbg_two,isbg_three,isbg_four/)
        self%fsnd(:)        = (/snd_one,snd_two,snd_three,snd_four/)
        self%fslt(:)        = (/slt_one,slt_two,slt_three,slt_four/)
        self%fcly(:)        = (/cly_one,cly_two,cly_three,cly_four/)

        allocate( self%layer_thickness(self%nlevel) ) 
        self%layer_thickness(:) = self%dz
     end subroutine read_mod_ped_namelist

end module modrich_run_options_mod
