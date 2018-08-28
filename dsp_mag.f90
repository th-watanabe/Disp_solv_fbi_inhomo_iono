MODULE dsp_mag
!------------------------------------------------------
!
!   Compute magnetospheric response, k_perp^2 Phi/J
!
!------------------------------------------------------

  use dsp_header

  implicit none

!  real(kind=DP), parameter :: ell   = 1000._DP
!  real(kind=DP)            :: kgp   = 0._DP / ell**2


    complex(kind=DP), parameter :: ui = ( 0._DP, 1._DP )
    real(kind=DP),    parameter :: pi = 3.141592653589793_DP
    integer,          parameter :: n_zeta = 1000

    real(kind=DP),    parameter :: zl = 1000._DP, zeta_0 = pi / 9._DP
    real(kind=DP),    parameter :: beta = 0.d+0, dpdl =-1._DP

    real(kind=DP),    save :: k_perp, theta

    integer, save :: isw


    data  k_perp / 1.0_DP /
    data  theta  / 0.0_DP /

    data isw / 0 /

    private

    public mag_response, mag_set_k_perp, mag_set_theta, mag_init, mag_set_isw


CONTAINS

  SUBROUTINE mag_init

      write(olog,fmt="(a,1p,1e15.7)")
      write(olog,fmt="(a,1p,1e15.7)") "# zl     = ", zl
      write(olog,fmt="(a,1p,1e15.7)") "# zeta_0 = ", zeta_0
      write(olog,fmt="(a,1p,2e15.7)") "# beta   = ", beta
      write(olog,fmt="(a,1p,1e15.7)") "# dpdl   = ", dpdl

  END SUBROUTINE mag_init


  SUBROUTINE mag_set_k_perp( k_perp_new )

    real(kind=DP) :: k_perp_new

      k_perp = k_perp_new

  END SUBROUTINE mag_set_k_perp


  SUBROUTINE mag_set_theta( theta_new )

    real(kind=DP) :: theta_new

      theta = theta_new

  END SUBROUTINE mag_set_theta


  SUBROUTINE mag_set_isw( isw_new )

    integer :: isw_new

      isw = isw_new

  END SUBROUTINE mag_set_isw



  SUBROUTINE mag_response( z, zres, psi_m, phi_m )

    complex(kind=DP), intent(in)  :: z
    complex(kind=DP), intent(out) :: zres, psi_m, phi_m


    complex(kind=DP) :: psi, chi, phi, omg, bph


    complex(kind=DP), dimension(0:n_zeta) :: zpsi, zphi
    real(kind=DP),    dimension(0:n_zeta) :: zeta, h_nu, h_ph, valp, sss, b0, dzds

    real(kind=DP),    dimension(0:n_zeta) :: kgp, kperp_sq, kphi_sq

    real(kind=DP)    :: dzeta, al, tau
    real(kind=DP),    dimension(4) :: tau_k, tau_w


    complex(kind=DP), dimension(4) :: chi_k, chi_w, bph_k, bph_w
    real(kind=DP),    dimension(4) :: ck, ci

    real(kind=DP)    :: w1, w2, w10, w20, dzeta2, fc1, fc2, fc3

    integer       :: i, j

!    integer, save :: isw
!    data isw / 1 /


!      dzeta   = ( pi*0.5_DP - zeta_0 ) / real( n_zeta, kind=DP )
!
!      zeta(0) = zeta_0
!      h_nu(0) = sin( zeta(0) )**3 / sqrt( 1._DP + 3._DP*cos( zeta(0) )**2 )
!      h_ph(0) = sin( zeta(0) )**3
!
!      w10     = sqrt( 3._DP )*cos( zeta(0) )
!      w20     = sqrt( 1._DP + w10**2 )
!      al      = zl*2._DP*sqrt( 3._DP ) / ( w10*w20 + log( abs( w10+w20 ) ) )
!
!        do i = 1, n_zeta
!          zeta(i) = zeta_0 + dzeta*dble(i)
!          h_nu(i) = sin( zeta(i) )**3                     &
!                  / sqrt( 1._DP + 3._DP*cos( zeta(i) )**2 ) &
!                  / h_nu(0)
!          h_ph(i) = sin( zeta(i) )**3 / h_ph(0)
!        end do
!
!        h_nu(0) = 1._DP
!        h_ph(0) = 1._DP
!
!
!        do i = 0, n_zeta
!          w1      = sqrt( 3._DP ) * cos( zeta(i) )
!          w2      = sqrt( 1._DP + w1**2 )
!          sss(i)  = al / ( 2._DP*sqrt(3._DP) )              &
!                  * ( w10 * w20 + log( abs( w10 + w20 ) )   &
!                    - w1  * w2  - log( abs( w1  + w2  ) ) )
!!!!          valp(i) = 1._DP
!!!!          valp(i) = 1._DP / ( h_nu(i) * h_ph(i) )
!!
!! rho is proportional to B
!          valp(i) = 1._DP / sqrt( h_nu(i) * h_ph(i) )
!
!!!!          valp(i) = 1._DP / ( h_nu(i) * h_ph(i) )**(1.d0/3.d0)
!!!!          valp(i) = 1._DP / sqrt( sqrt( h_nu(i) * h_ph(i) ) )
!          b0(i)   = 1._DP / ( h_nu(i) * h_ph(i) )
!          dzds(i) = 1._DP / ( al * sin( zeta(i) ) )         &
!                  / sqrt( 1._DP + 3._DP*cos( zeta(i) )**2 )
!
!          kperp_sq(i) = k_perp**2                       &
!                      * ( ( cos(theta) / h_ph(i) )**2   &
!                        + ( sin(theta) / h_nu(i) )**2 )
!          kphi_sq(i) = k_perp**2                        &
!                      *   ( cos(theta) / h_ph(i) )**2
!
!          kgp(i)  = 1._DP / ( al**2 * b0(i)**3 )        &
!                  / sin( zeta(i) )**4                   &
!                  * ( 1._DP + cos( zeta(i) )**2 )       &
!                  / ( 1._DP + 3._DP*cos( zeta(i) )**2 ) &
!                  *  kphi_sq(i)
!!!!                  *  kphi_sq(i) / kperp_sq(i)        ! fixed (Mar 30, 2016)
!        end do
!          sss(0)       = 0._DP
!          sss(n_zeta)  = zl
!
! for comparison with the slab case
          dzeta   = zl / dble( n_zeta )
        do i = 0, n_zeta
          zeta(i) = dzeta*dble(i)
          b0(i)   = 1._DP
          h_nu(i) = 1._DP
          h_ph(i) = 1._DP
          dzds(i) = 1._DP
          valp(i) = 1._DP

          sss(i)  = zeta(i)

          kperp_sq(i) = k_perp**2                       &
                      * ( ( cos(theta) / h_ph(i) )**2   &
                        + ( sin(theta) / h_nu(i) )**2 )
          kphi_sq(i) = k_perp**2                        &
                      *   ( cos(theta) / h_ph(i) )**2
          kgp(i)  = 0._DP
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

! boundary condition on the equator
! symmetric for phi
        psi   = ( 0._DP, 0._DP )
        bph   = ( 1._DP, 0._DP )
!!!! anti-symmetric for phi
!!!        psi   = cmplx( 0._DP, 1._DP/b0(n_zeta) )
!!!        bph   = ( 0._DP, 0._DP )

        zpsi(:) = ( 0._DP, 0._DP )
        zphi(:) = ( 0._DP, 0._DP )


!!!        kperp_sq = k_perp**2                           &
!!!                 * ( ( cos(theta) / h_ph(n_zeta) )**2  &
!!!                   + ( sin(theta) / h_nu(n_zeta) )**2 )

        chi   = - kperp_sq(n_zeta) * psi
        omg   = - kperp_sq(n_zeta) * bph / b0(n_zeta)

        zpsi(n_zeta) = psi
        zphi(n_zeta) = bph / b0(n_zeta)

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
        valp(:)   = valp(:) * tau / zl
        


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
!!!              kperp_sq = k_perp**2                     &
!!!                   * ( ( cos(theta) / h_ph(i  ) )**2   &
!!!                     + ( sin(theta) / h_nu(i  ) )**2 )
              fc1  = - kperp_sq(i  ) /(valp(i  )**2 * dzds(i  ) * b0(i  ))
              fc2  = - b0(i  ) / (dzds(i  )*kperp_sq(i  ))
              fc3  =   kgp(i  ) * beta * dpdl
            else if ( j == 4 ) then
!!!              kperp_sq = k_perp**2                     &
!!!                   * ( ( cos(theta) / h_ph(i-2) )**2   &
!!!                     + ( sin(theta) / h_nu(i-2) )**2 )
              fc1  = - kperp_sq(i-2) /(valp(i-2)**2 * dzds(i-2) * b0(i-2))
              fc2  = - b0(i-2) / (dzds(i-2)*kperp_sq(i-2))
              fc3  =   kgp(i-2) * beta * dpdl
            else
!!!              kperp_sq = k_perp**2                     &
!!!                   * ( ( cos(theta) / h_ph(i-1) )**2   &
!!!                     + ( sin(theta) / h_nu(i-1) )**2 )
              fc1  = - kperp_sq(i-1) /(valp(i-1)**2 * dzds(i-1) * b0(i-1))
              fc2  = - b0(i-1) / (dzds(i-1)*kperp_sq(i-1))
              fc3  =   kgp(i-1) * beta * dpdl
            end if

            chi_k(j) = ui * bph_w(j) * ( fc1 * z + fc3 / z )
            bph_k(j) = ui * chi_w(j) *   fc2 * z
          end do

          do j = 1, 4
            chi   = chi + dzeta2 * ck(j)*chi_k(j)
            bph   = bph + dzeta2 * ck(j)*bph_k(j)
          end do

!!!          kperp_sq = k_perp**2                           &
!!!                     * ( ( cos(theta) / h_ph(i-2) )**2   &
!!!                       + ( sin(theta) / h_nu(i-2) )**2 )

          psi  = - chi / kperp_sq(i-2)
          omg  = - kperp_sq(i-2) * bph / b0(i-2)

          zpsi(i-2) = psi
          zphi(i-2) = bph / b0(i-2)

        end do

        phi   = - omg / k_perp**2

        zres = k_perp**2 * phi / ( b0(0) * chi )

        psi_m = psi
        phi_m = phi


        if( isw == 1 ) then
            write(unit=20,fmt="(a, 1p, 3e15.7)") "# k_perp, z = ", k_perp, z
          do i = 0, n_zeta, 2
            write(unit=20,fmt="(1p, 11e15.7)") &
!!!              sss(i), zeta(i), b0(i), valp(i), (b0(i)/valp(i))**2, kgp(i)
              zeta(i), zpsi(i), zphi(i), b0(i), h_ph(i), h_nu(i), sss(i), valp(i), kgp(i)
          end do
            write(20,*)
            isw = 0
        end if


      return

  END SUBROUTINE mag_response


END MODULE dsp_mag
