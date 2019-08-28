module constants_mod


  implicit none

  integer, parameter :: uk  = SELECTED_REAL_KIND(12, 50) 
  integer, parameter :: nhorz_max = 4
  real(kind=uk), parameter :: smp_min   = -1.0e8                                  ,&
                              eps       = epsilon(1._uk)                          ,&
                              zero      = 0._uk                                   ,&
                              one       = 1._uk                                   ,&
                              two       = 2._uk                                   ,&
                              three     = 3._uk                                   ,&
                              ten       = 10._uk                                  ,&
                              half      = 0.5_uk                                  ,&
                              half_eps  = 10._uk ** (-real(uk,kind=uk) )          ,&
                              close_val = 1.0e-5_uk                           
end module constants_mod


