function zero ( a, b, machep, t, f )

!*****************************************************************************80
!
!! ZERO seeks the root of a function F(X) in an interval [A,B].
!
!  Discussion:
!
!    The interval [A,B] must be a change of sign interval for F.
!    That is, F(A) and F(B) must be of opposite signs.  Then
!    assuming that F is continuous implies the existence of at least
!    one value C between A and B for which F(C) = 0.
!
!    The location of the zero is determined to within an accuracy
!    of 6 * MACHEPS * abs ( C ) + 2 * T.
!
!    Thanks to Thomas Secretin for pointing out a transcription error in the
!    setting of the value of P, 11 February 2013.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 February 2013
!
!  Author:
!
!    Original FORTRAN77 version by Richard Brent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the change of 
!    sign interval.
!
!    Input, real ( kind = 8 ) MACHEP, an estimate for the relative machine
!    precision.
!
!    Input, real ( kind = 8 ) T, a positive error tolerance.
!
!    Input, external real ( kind = 8 ) F, the name of a user-supplied
!    function, of the form "FUNCTION F ( X )", which evaluates the
!    function whose zero is being sought.
!
!    Output, real ( kind = 8 ) ZERO, the estimated value of a zero of
!    the function F.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) e
  real ( kind = 8 ) f
  real ( kind = 8 ) fa
  real ( kind = 8 ) fb
  real ( kind = 8 ) fc
  real ( kind = 8 ) m
  real ( kind = 8 ) machep
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) sa
  real ( kind = 8 ) sb
  real ( kind = 8 ) t
  real ( kind = 8 ) tol
  real ( kind = 8 ) zero
!
!  Make local copies of A and B.
!
  sa = a
  sb = b
  fa = f ( sa )
  fb = f ( sb )

  c = sa
  fc = fa
  e = sb - sa
  d = e

  do

    if ( abs ( fc ) < abs ( fb ) ) then

      sa = sb
      sb = c
      c = sa
      fa = fb
      fb = fc
      fc = fa

    end if

    tol = 2.0D+00 * machep * abs ( sb ) + t
    m = 0.5D+00 * ( c - sb )

    if ( abs ( m ) <= tol .or. fb == 0.0D+00 ) then
      exit
    end if

    if ( abs ( e ) < tol .or. abs ( fa ) <= abs ( fb ) ) then

      e = m
      d = e

    else

      s = fb / fa

      if ( sa == c ) then

        p = 2.0D+00 * m * s
        q = 1.0D+00 - s

      else

        q = fa / fc
        r = fb / fc
        p = s * ( 2.0D+00 * m * q * ( q - r ) - ( sb - sa ) * ( r - 1.0D+00 ) )
        q = ( q - 1.0D+00 ) * ( r - 1.0D+00 ) * ( s - 1.0D+00 )

      end if

      if ( 0.0D+00 < p ) then
        q = - q
      else
        p = - p
      end if

      s = e
      e = d

      if ( 2.0D+00 * p < 3.0D+00 * m * q - abs ( tol * q ) .and. &
        p < abs ( 0.5D+00 * s * q ) ) then
        d = p / q
      else
        e = m
        d = e
      end if

    end if

    sa = sb
    fa = fb

    if ( tol < abs ( d ) ) then
      sb = sb + d
    else if ( 0.0D+00 < m ) then
      sb = sb + tol
    else
      sb = sb - tol
    end if

    fb = f ( sb )

    if ( ( 0.0D+00 < fb .and. 0.0D+00 < fc ) .or. &
         ( fb <= 0.0D+00 .and. fc <= 0.0D+00 ) ) then
      c = sa
      fc = fa
      e = sb - sa
      d = e
    end if

  end do

  zero = sb

  return
end
