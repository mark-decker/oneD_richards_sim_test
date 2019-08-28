module modrich_run_options_mod

   use constants_mod
   use abstract_run_options_mod
   use option_vars_mod

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

      real(kind=uk) :: crnk=1.0

      contains
         procedure :: read_namelist => read_modrich_namelist
         procedure :: construct => construct_mod_opts

   end type mod_opts


   !read soil texture (mineral) and calc params
   type, extends(mod_opts) :: mod_ped_opts
      real(kind=uk), dimension(:), pointer :: fsnd,fslt,fcly
      contains
         procedure :: read_namelist => read_mod_ped_namelist
         procedure :: construct => construct_mod_ped_opts
   end type

   contains

        
     subroutine read_modrich_namelist(self)
        class(mod_opts), intent(inout) :: self

        character(len=265) :: my_iomsg
        integer :: my_iostat

        logical :: file_found
        character(len=*), parameter :: sub_name='read_namelist'

        integer :: myunit

        namelist/modrich_isoil/isoil_one,isoil_two,isoil_three,isoil_four,&
                             isbg_one,isbg_two,isbg_three,isbg_four,&
                             dtime_max,dtime_min,ntime,dtime,&
                             dz,q_infil,q_drain,initial_wtd,&
                             nlevel,nhorz,rel_sat_max,rel_sat_min

        inquire( file=trim(self%namelist_file), exist=file_found )
        if (.not.file_found) then
           write(*,*) 'The input file '//trim(self%namelist_file)//' cannot be found.'
           stop
        end if


        !the nml to read and the data depends on type of self
        open(newunit=myunit,file=trim(self%namelist_file),status='old')
        read(myunit,nml=modrich_isoil)
        close(myunit)
   

        self%nml_done = .true.

     end subroutine read_modrich_namelist


      subroutine construct_mod_opts(self)
         class(mod_opts), intent(inout) :: self

         real(kind=uk), dimension(:), pointer, save :: saved_local_thickness 
         integer, dimension(:), pointer, save :: saved_local_isoils, &
                                                         saved_local_isoils_beg
         integer :: i

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
        self%rel_sat_max = rel_sat_max
        self%rel_sat_min = rel_sat_min
        self%delta_liq_max = delta_liq_max
        self%delta_se_max = delta_se_max

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

     end subroutine construct_mod_opts


     subroutine read_mod_ped_namelist(self)
        class(mod_ped_opts), intent(inout) :: self

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

        inquire( file=trim(self%namelist_file), exist=file_found )
        if (.not.file_found) then
           write(*,*) 'The input file '//trim(self%namelist_file)//' cannot be found.'
           stop
        end if

        open(newunit=myunit,file=trim(self%namelist_file),status='old')
        read(myunit,nml=modrich_pedt,iostat=my_iostat,iomsg=my_iomsg)
        close(myunit)


        self%nml_done = .true.

     end subroutine read_mod_ped_namelist

      subroutine construct_mod_ped_opts(self)
         class(mod_ped_opts), intent(inout) :: self

         real(kind=uk), dimension(:), pointer, save :: saved_local_thickness , &
                                                         saved_local_sand,saved_local_clay,saved_local_silt
         integer, dimension(:), pointer, save :: saved_local_isoils, &
                                                         saved_local_isoils_beg
         integer :: i


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

        if (associated(saved_local_isoils)) nullify(saved_local_isoils)
        allocate(saved_local_isoils(4))

        if (associated(saved_local_isoils_beg)) nullify(saved_local_isoils_beg)
        allocate(saved_local_isoils_beg(4))

        if (associated(saved_local_thickness)) nullify(saved_local_thickness)
        allocate(saved_local_thickness(4))

        allocate(saved_local_thickness(nhorz))

        self%isoils       => saved_local_isoils
        self%isoils_beg   => saved_local_isoils_beg

        if (associated(self%fcly)) nullify(self%fcly)
        if (associated(self%fslt)) nullify(self%fslt)
        if (associated(self%fsnd)) nullify(self%fsnd)

        allocate(saved_local_sand(max(4,nhorz)))
        allocate(saved_local_silt(max(4,nhorz)))
        allocate(saved_local_clay(max(4,nhorz)))

        do i=1,nhorz
           select case(i)
           case(1)
              saved_local_sand(i) = snd_one
              saved_local_silt(i) = slt_one
              saved_local_clay(i) = cly_one
           case(2)
              saved_local_sand(i) = snd_two
              saved_local_silt(i) = slt_two
              saved_local_clay(i) = cly_two
           case(3)
              saved_local_sand(i) = snd_three
              saved_local_silt(i) = slt_three
              saved_local_clay(i) = cly_three
           case(4)
              saved_local_sand(i) = snd_four
              saved_local_silt(i) = slt_four
              saved_local_clay(i) = cly_four
           case default
              write(*,*) 'currently only four horizons are coded for.  Add more'
           end select
        end do

        saved_local_thickness(:) = spread(dz,1,nhorz)

        self%fsnd => saved_local_sand
        self%fslt => saved_local_silt
        self%fcly => saved_local_clay

        self%layer_thickness => saved_local_thickness



     end subroutine construct_mod_ped_opts


end module modrich_run_options_mod
