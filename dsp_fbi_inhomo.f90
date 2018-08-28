MODULE dsp_fbi

  use dsp_header

  implicit none

  real(kind=DP), parameter :: pi = 3.141592653589793d0

  real(kind=DP),    save :: k_perp, e0, dfs, alp, sigma_p
  real(kind=DP),    save :: theta, mp
  real(kind=DP),    save :: height, nu_in, et_in
  complex(kind=DP), save :: phi_jpara

  integer, save :: isw

! initial values !
    data phi_jpara / ( 0.0_DP, 0._DP ) /
    data k_perp  / 1.0_DP    /
!!!    data e0      / 0.001_DP  / !!! Feb06
    data e0      / 0.0032467_DP  /
!    data dfs     / 2.d-5     /
    data dfs     / 0.d+0     /
!    data alp     / 7.d-4     /
    data alp     / 0.d+0     /
    data sigma_p / 5._DP     /

    data theta   / 1.5707963267949d0 /
!!!    data mp      / 0.5_DP    / !!! Feb06
    data mp      / 0.154_DP    /
    data height  / 1.0_DP    /
    data nu_in   / 10._DP    /
!    data nu_in   / 100._DP    /
    data et_in   / 0._DP    /

    data isw / 0 /

  private

  public fbi_init, fbi_func, fbi_set_k_perp, fbi_set_phi_jpara, fbi_set_isw
  public fbi_set_theta
  public fbi_func_inhomo


CONTAINS


  SUBROUTINE fbi_init


      write(olog,fmt="(a,1p,1e15.7)")
      write(olog,fmt="(a,1p,2e15.7)") "# phi_jpara = ", phi_jpara
      write(olog,fmt="(a,1p,1e15.7)") "# k_perp = ", k_perp
      write(olog,fmt="(a,1p,1e15.7)") "# e0     = ", e0
      write(olog,fmt="(a,1p,1e15.7)") "# dfs    = ", dfs
      write(olog,fmt="(a,1p,1e15.7)") "# alp    = ", alp
      write(olog,fmt="(a,1p,1e15.7)") "# sigma_p= ", sigma_p
      write(olog,fmt="(a,1p,1e15.7)") "# theta  = ", theta
      write(olog,fmt="(a,1p,1e15.7)") "# mp     = ", mp
      write(olog,fmt="(a,1p,1e15.7)") "# height = ", height
      write(olog,fmt="(a,1p,1e15.7)") "# nu_in  = ", nu_in 
      write(olog,fmt="(a,1p,1e15.7)") "# et_in  = ", et_in 
      write(olog,fmt="(a,1p,1e15.7)")


  END SUBROUTINE fbi_init


  SUBROUTINE fbi_set_k_perp( k_perp_new )

    real(kind=DP) :: k_perp_new

      k_perp = k_perp_new

  END SUBROUTINE fbi_set_k_perp


  SUBROUTINE fbi_set_theta( theta_new )

    real(kind=DP) :: theta_new

!      theta = theta_new
!
      print *, "# theta is unchanged, and fixed to pi/2"


  END SUBROUTINE fbi_set_theta


  SUBROUTINE fbi_set_phi_jpara( phi_jpara_new )

    complex(kind=DP) :: phi_jpara_new

      phi_jpara = phi_jpara_new

  END SUBROUTINE fbi_set_phi_jpara


  SUBROUTINE fbi_set_isw( isw_new )

    integer :: isw_new

      isw = isw_new

  END SUBROUTINE fbi_set_isw


  SUBROUTINE fbi_func( omega )

    complex(kind=DP), intent(out) :: omega

    complex(kind=DP) :: ui = ( 0._DP, 1._DP )

!!!      omega = ( k_perp * e0 - ui * k_perp**2 * dfs ) &

!!!    write(6,*) "theta = ", theta

      omega = ( k_perp * e0 * ( mp * sin(theta) - cos(theta) ) &
              - ui * k_perp**2 * dfs ) &
            / ( 1._DP + phi_jpara * sigma_p )        &
            - ui * alp * 2._DP


  END SUBROUTINE fbi_func


  SUBROUTINE fbi_func_inhomo( z, psi_m, phi_m, psi_bottom )

    complex(kind=DP), intent(in)  :: z, psi_m, phi_m
    complex(kind=DP), intent(out) :: psi_bottom

    complex(kind=DP) :: ui = ( 0._DP, 1._DP )

!!!      omega = ( k_perp * e0 - ui * k_perp**2 * dfs ) &

!!!    write(6,*) "theta = ", theta

!! height integrated ionosphere model
!      omega = ( k_perp * e0 * ( mp * sin(theta) - cos(theta) ) &
!              - ui * k_perp**2 * dfs ) &
!            / ( 1._DP + phi_jpara * sigma_p )        &
!            - ui * alp * 2._DP


    integer, parameter :: n_zeta = 4000 

    complex(kind=DP) :: psi, chi, phi, omg, bph


    complex(kind=DP), dimension(0:n_zeta) :: zpsi, zphi, zdns
    real(kind=DP),    dimension(0:n_zeta) :: zeta, h_nu, h_ph, valp, sss, b0, dzds

    real(kind=DP),    dimension(0:n_zeta) :: kgp, kperp_sq, kphi_sq
    real(kind=DP),    dimension(0:n_zeta) :: mp_z, nu_z, dn_z, ff_z, al_z, et_z

    real(kind=DP)    :: dzeta, al, tau
    real(kind=DP),    dimension(4) :: tau_k, tau_w


    complex(kind=DP), dimension(4) :: chi_k, chi_w, bph_k, bph_w
    real(kind=DP),    dimension(4) :: ck, ci

    real(kind=DP)    :: dzeta2
    complex(kind=DP) :: fc1, fc2

    integer       :: i, j

!    integer, save :: isw
!    data isw / 1 /


      dzeta   = height / real( n_zeta, kind=DP )

      zeta(0) = - height

! slab geometry
        do i = 0, n_zeta
          zeta(i)     = zeta(0) + dzeta*dble(i)
          valp(i)     = 1._DP
          b0(i)       = 1._DP
          h_nu(i)     = 1._DP
          h_ph(i)     = 1._DP
          dzds(i)     = 1._DP
          kperp_sq(i) = k_perp**2
        end do

! inhomogeneous ionospher
        do i = 0, n_zeta
          nu_z(i) = nu_in * ( 1._DP - cos( pi * zeta(i) / height ) )
          mp_z(i) = nu_z(i) / ( 1 + nu_z(i)**2 )
!!!          mp_z(i) = 2._DP * mp * ( sin( pi * zeta(i) / height ) )**2 - 1.d-14
!!!          nu_z(i) = ( 1._DP - sqrt( 1._DP - 4._DP*mp_z(i)**2 ) ) / ( 2._DP*mp_z(i) )
          dn_z(i) = sigma_p / ( mp * height )
          ff_z(i) = ( 1._DP - nu_z(i)**2 ) / ( 1._DP + nu_z(i)**2 )**2
          al_z(i) = alp
          et_z(i) = et_in * ( 1._DP - cos( pi * zeta(i) / height ) )
        end do


! coefficients used in the RK method
          ck(1)   = 1._DP/6._DP
          ck(2)   = 1._DP/3._DP
          ck(3)   = 1._DP/3._DP
          ck(4)   = 1._DP/6._DP

          ci(1)   = 0._DP
          ci(2)   = 0.5_DP
          ci(3)   = 0.5_DP
          ci(4)   = 1._DP


        chi_k = ( 0._DP, 0._DP )
        chi_w = ( 0._DP, 0._DP )
        bph_k = ( 0._DP, 0._DP )
        bph_w = ( 0._DP, 0._DP )

! boundary codition given by the magnetosphere
!        bph   = ( 1._DP, 0._DP )
!        psi   = - bph / phi_jpara

        bph   = phi_m
        psi   = psi_m

        zpsi(:) = ( 0._DP, 0._DP )
        zphi(:) = ( 0._DP, 0._DP )
        zdns(:) = ( 0._DP, 0._DP )


!!!        kperp_sq = k_perp**2                           &
!!!                 * ( ( cos(theta) / h_ph(n_zeta) )**2  &
!!!                   + ( sin(theta) / h_nu(n_zeta) )**2 )

        chi   = - kperp_sq(n_zeta) * psi
        omg   = - kperp_sq(n_zeta) * bph / b0(n_zeta)

        zpsi(n_zeta) = psi
        zphi(n_zeta) = bph / b0(n_zeta)
        zdns(n_zeta) = ( 0._DP, 0._DP )

        dzeta2 = dzeta * 2._DP


        tau      = 0._DP
        tau_k(:) = 0._DP
        tau_w(:) = 0._DP

        do i = n_zeta, 2, -2

          do j = 1, 4
            if( j == 1 ) then
              tau_w(j) = tau
            else
              tau_w(j) = tau + ci(j) * tau_k(j-1) * dzeta2
            end if

            if( j == 1 ) then
              tau_k(j) = 1._DP / ( valp(i  ) * dzds(i  ) )
            else if ( j == 4 ) then
              tau_k(j) = 1._DP / ( valp(i-2) * dzds(i-2) )
            else
              tau_k(j) = 1._DP / ( valp(i-1) * dzds(i-1) )
            end if

          end do

          do j = 1, 4
            tau   = tau + dzeta2 * ck(j)*tau_k(j)
          end do

        end do

!!!        write(unit=olog,fmt="(a,1p,e15.7)") "# tau = ", tau

! re-scaling of valp
        valp(:)   = valp(:) * tau / height
        


        do i = n_zeta, 2, -2

! Runge-Kutta starts

          do j = 1, 4
            if( j .eq. 1 ) then
              chi_w(j) = chi
              bph_w(j) = bph
            else
              chi_w(j) = chi + ci(j) * chi_k(j-1) * dzeta2
              bph_w(j) = bph + ci(j) * bph_k(j-1) * dzeta2
            end if

            if( j == 1 ) then
!!!              fc1  = - kperp_sq(i  ) /(valp(i  )**2 * dzds(i  ) * b0(i  ))

              fc1  = + kperp_sq(i  ) / ( 1._DP - k_perp * mp_z(i  ) * e0 / ( z + ui*2._DP*al_z(i  ) ) ) &
                   * ( - ui * z * ff_z(i  ) / valp(i  )**2  + mp_z(i  ) * dn_z(i  ) ) &
                   / ( dzds(i  ) * b0(i  ) )
!!!              fc2  = - b0(i  ) / (dzds(i  )*kperp_sq(i  ))
              fc2  = - b0(i  ) / (dzds(i  )*kperp_sq(i  )) * ( ui * z + et_z(i  ) * kperp_sq(i  ) )

            else if ( j == 4 ) then
!!!              fc1  = - kperp_sq(i-2) /(valp(i-2)**2 * dzds(i-2) * b0(i-2))
              fc1  = + kperp_sq(i-2) / ( 1._DP - k_perp * mp_z(i-2) * e0 / ( z + ui*2._DP*al_z(i-2) ) ) &
                   * ( - ui * z * ff_z(i-2) / valp(i-2)**2  + mp_z(i-2) * dn_z(i-2) ) &
                   / ( dzds(i-2) * b0(i-2) )
!!!              fc2  = - b0(i-2) / (dzds(i-2)*kperp_sq(i-2))
              fc2  = - b0(i-2) / (dzds(i-2)*kperp_sq(i-2)) * ( ui * z + et_z(i-2) * kperp_sq(i-2) )

            else
!!!              fc1  = - kperp_sq(i-1) /(valp(i-1)**2 * dzds(i-1) * b0(i-1))
              fc1  = + kperp_sq(i-1) / ( 1._DP - k_perp * mp_z(i-1) * e0 / ( z + ui*2._DP*al_z(i-1) ) ) &
                   * ( - ui * z * ff_z(i-1) / valp(i-1)**2  + mp_z(i-1) * dn_z(i-1) ) &
                   / ( dzds(i-1) * b0(i-1) )
!!!              fc2  = - b0(i-1) / (dzds(i-1)*kperp_sq(i-1))
              fc2  = - b0(i-1) / (dzds(i-1)*kperp_sq(i-1)) * ( ui * z + et_z(i-1) * kperp_sq(i-1) )

            end if

!!!            chi_k(j) = ui * bph_w(j) *   fc1 * z
            chi_k(j) =      bph_w(j) *   fc1

!!!            bph_k(j) = ui * chi_w(j) *   fc2 * z
            bph_k(j) =      chi_w(j) *   fc2
          end do

          do j = 1, 4
            chi   = chi + dzeta2 * ck(j)*chi_k(j)
            bph   = bph + dzeta2 * ck(j)*bph_k(j)
          end do

          psi  = - chi / kperp_sq(i-2)
          omg  = - kperp_sq(i-2) * bph / b0(i-2)

          zpsi(i-2) = psi
          zphi(i-2) = bph / b0(i-2)
!!!          zdns(i-2) = bph / b0(i-2) * fc1 * ui / z
!!!          zdns(i-2) = bph / b0(i-2) * fc1 * ui / ( z + ui*2._DP*al_z(i-2) )
!
! bug fixed (Jul 17, 2018)
          zdns(i-2) = bph / b0(i-2) * fc1 * ui / ( z + ui*2._DP*al_z(i-2) ) / ( -dn_z(i-2) )

        end do

        phi   = - omg / k_perp**2

!
! The following omega is to find zpsi(0) = 0 by minimizing 
! omega_i / omega_m - 1 in mi_couple_func
!
        psi_bottom = zpsi(0)


        if( isw == 1 ) then
            write(unit=30,fmt="(a, 1p, 3e15.7)") "# k_perp, z = ", k_perp, z
          do i = 0, n_zeta, 2
            write(unit=30,fmt="(1p, 13e15.7)") &
!!!              zeta(i), valp(i), mp_z(i), nu_z(i), dn_z(i), ff_z(i), zphi(i), zpsi(i)
              zeta(i), zpsi(i), zphi(i), zdns(i), valp(i), mp_z(i), nu_z(i), dn_z(i), ff_z(i), al_z(i)
          end do
            write(30,*)
            isw = 0
        end if


      return

  END SUBROUTINE fbi_func_inhomo


END MODULE dsp_fbi
