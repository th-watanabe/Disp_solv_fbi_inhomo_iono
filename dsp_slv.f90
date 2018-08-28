MODULE dsp_slv

  use dsp_header

  implicit none

  private

  public slv_newton


CONTAINS


  SUBROUTINE slv_newton( fun, zin, zout, iconv )

    complex(kind=DP), external :: fun
    complex(kind=DP), intent(in)  :: zin
    complex(kind=DP), intent(out) :: zout
    integer, intent(inout) :: iconv

    complex(kind=DP) :: z0, z1, z2, dz
    complex(kind=DP) :: f0, f1, f2
    real(kind=DP)    :: eps, eps0 = 0.000001_DP
    real(kind=DP)    :: eps_min = 1.d-99
    integer          :: ic, nc = 20


      z0   = zin

      iconv = 0

      eps  = 1._DP
      ic   = 0

      do while ( eps > eps0  .AND.  ic < nc )

        ic = ic + 1

        dz = z0 * 0.00000001_DP

        z1 = z0 - dz
        z2 = z0 + dz

        f0 = fun( z0 )
        f1 = fun( z1 )
        f2 = fun( z2 )

        zout = z0 - f0 * ( z2 - z1 ) / ( f2 - f1 )
        eps  = abs( zout - z0 ) / abs( z0 )
        z0   = zout

!!!        print *, "# ic, zout, f0 = ", ic, zout, f0

      end do

        f0 = fun( z0 )

      if( ic >= nc ) then
        print *, "# no convergence in iteration: ic = ", ic
        iconv = 1
!!!        stop   !!! July 11
      end if

      if( abs( real (zout) ) < eps_min ) then
        zout = cmplx( 0._DP, aimag(zout), kind=DP )
      end if

      if( abs( aimag(zout) ) < eps_min ) then
        zout = cmplx( real (zout), 0._DP, kind=DP )
      end if


  END SUBROUTINE slv_newton


END MODULE dsp_slv
