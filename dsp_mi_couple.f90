MODULE dsp_mi_couple

  use dsp_header
  use dsp_slv
  use dsp_fbi
  use dsp_mag

  implicit none

  real(kind=DP), parameter :: ell   = 1000._DP
  real(kind=DP)            :: theta = 0._DP

  private

  public mi_couple_find, mi_couple_func


CONTAINS


  SUBROUTINE mi_couple_find

    complex(kind=DP) :: omega, k_para

    complex(kind=DP) :: zin, zout
    real(kind=DP)    :: k_perp

    complex(kind=DP) :: phi_jpara, psi_bottom, psi_m, phi_m


    real(kind=DP)    :: dk, dth 
    integer          :: i, j, k, m

    real(kind=DP)    :: pi = 3.141592653589793_DP

    integer, parameter :: iloop=500, jloop=4, nth = 1
    complex(kind=DP), dimension (0:jloop,0:nth-1,-1:1) :: zsave

    integer :: iconv = 0

      zin    = cmplx( pi/ell, 0.0_DP )

      dk     = pi * 100._DP / real( iloop, kind=DP )
      dth    = pi * 0.5_DP  / real( nth, kind=DP )
 
      k_perp = dk


      call fbi_set_k_perp( k_perp )

      call fbi_init
      call mag_init
      call mi_couple_init


      do m = 0, nth-1

!!!      do k = 1, -1, -2
      do k = 1, 1

      do j = 0, jloop

!!!        print *, "# j = ", j

        if( m == 0 ) then
!!!          zin = cmplx( pi*(j+0.5_DP)*k/ell, 0._DP, kind=DP )
          zin = cmplx( pi*(j+1)*k/ell, 1.d-5, kind=DP )
        else
          zin = zsave(j,m-1,k)
        end if

        theta = dth * real( m, kind=DP ) + pi*0.5_DP

        call fbi_set_theta( theta )
        call mag_set_theta( theta )


        do i = iloop, 2, -1
!!!        do i = 1, iloop

!!!          print *, "# i = ", i

          k_perp = dk*real(i,kind=DP)
          call fbi_set_k_perp( k_perp )
          call mag_set_k_perp( k_perp )

          call slv_newton( mi_couple_func, zin, zout, iconv )

          if( iconv /= 0 ) exit

          k_para = zout
!          call fbi_func( omega )
          omega  = zout

          write(olog,fmt="(1p,6e15.7)") k_para, k_perp, theta, omega

          zin = zout

          if( i == iloop ) then
            zsave(j,m,k) = zout
          end if

          if( i == 20 ) then
!          if( i == 50 ) then
!          if( i == 450 ) then
            call fbi_set_isw( 1 )
            call mag_set_isw( 1 )
            call mag_response( omega, phi_jpara, psi_m, phi_m )
!!!            call fbi_func( omega )
            call fbi_func_inhomo( omega, psi_m, phi_m, psi_bottom )
            call fbi_set_isw( 0 )
            call mag_set_isw( 0 )
          end if


        end do

          write(olog,fmt="(1p,5e15.7)")

      end do

      end do

!!!          write(olog,fmt="(1p,5e15.7)")

      end do


  END SUBROUTINE mi_couple_find


  SUBROUTINE mi_couple_init

      write(olog,fmt="(a,1p,e15.7)")
      write(olog,fmt="(a,1p,e15.7)") "# ell   = ", ell

  END SUBROUTINE mi_couple_init


  FUNCTION mi_couple_func( zin )

    complex(kind=DP) :: mi_couple_func
    complex(kind=DP), intent(in) :: zin

    complex(kind=DP) :: k_para, omega_m, omega_i, kpl
    complex(kind=DP) :: phi_jpara
    complex(kind=DP) :: psi_bottom, psi_m, phi_m

    complex(kind=DP) :: ui = ( 0._DP, 1._DP )

    real(kind=DP)    :: aomg, delt

    complex(kind=DP) :: k_para_sq


!      k_para  = zin
!
!! test
!      omega_m = k_para
!! test
!
!      k_para_sq  = k_para **2 - kgp * cos(theta)**2
!
!      aomg    = abs( k_para_sq )
!      delt    = atan2( aimag( k_para_sq ), real( k_para_sq ) )
!      omega_m = sqrt( aomg ) * exp( ui*delt*0.5_DP ) &
!              * sign( 1._DP, real( k_para ) )
!
!!!      print *, "# k_para, omega_m = ", k_para, omega_m
!
!      kpl     = k_para * ell
!      phi_jpara = - ( exp(ui*kpl) + exp(-ui*kpl) ) &
!                  / ( exp(ui*kpl) - exp(-ui*kpl) )

      omega_m = zin
      call mag_response( omega_m, phi_jpara, psi_m, phi_m )

      call fbi_set_phi_jpara( phi_jpara )

!! for conventional dispersion relation
!      call fbi_func( omega_i )
!
!      mi_couple_func = omega_i / omega_m - 1._DP


! for inhomogeneous ionosphere
      call fbi_func_inhomo( omega_m, psi_m, phi_m, psi_bottom )

      mi_couple_func = psi_bottom

       

  END FUNCTION mi_couple_func


END MODULE dsp_mi_couple
