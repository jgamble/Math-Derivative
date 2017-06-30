package Math::Derivative;

use 5.010001;
use Exporter;

our @ISA = qw(Exporter);
our %EXPORT_TAGS = (all => [qw(
	Derivative1
	Derivative2
	centraldiff
	forwarddiff
	seconddx
)]);

our @EXPORT_OK = (@{$EXPORT_TAGS{all}});

our $VERSION = 0.04;

use strict;
use warnings;
use Carp;

=head1 NAME

Math::Derivative - Numeric 1st and 2nd order differentiation

=head1 SYNOPSIS

    use Math::Derivative qw(:all);

    @dydx = forwarddiff(\@x, \@y);

    @dydx = centraldiff(\@x, \@y);

    @d2ydx2 = seconddx(\@x, \@y);

    @dydx = Derivative1(\@x, \@y);     # A synonym for centraldiff()
    @d2ydx2 = Derivative2(\@x, \@y);   # A synonym for seconddx()

=head1 DESCRIPTION

This Perl package exports functions that numerically approximate first
and second order differentiation on vectors of data. The accuracy of
the approximation will depend upon the differences between the
successive values in the X array.

=head2 FUNCTIONS

The functions may be imported by name or by using the tag ":all".

=head3 forwarddiff()

    @dydx = forwarddiff(\@x, \@y);

Take the references to two arrays containing the x and y ordinates of
the data, and return an array of approximate first derivatives at the
given x ordinates, using the forward difference approximation.

=cut

sub forwarddiff
{
	my($x, $y) = @_;
	my @y2;
	my $n = $#{$x};

	croak "X and Y array lengths don't match." unless ($n == $#{$y});

	$y2[$n] = ($y->[$n] - $y->[$n-1])/($x->[$n] - $x->[$n-1]);

	for my $i (0 .. $n-1)
	{
		$y2[$i] = ($y->[$i+1] - $y->[$i])/($x->[$i+1] - $x->[$i]);
	}

	return @y2;
}

=head3 centraldiff()

    @dydx = centraldiff(\@x, \@y);

Take the references to two arrays containing the x and y ordinates of
the data, and return an array of approximate first derivatives at the
given x ordinates.

The algorithm used three data points to calculate the derivative, except
at the end points, where by necessity the forward difference algorithm
is used instead.

=cut

sub centraldiff
{
	my($x, $y) = @_;
	my @y2;
	my $n = $#{$x};

	croak "X and Y array lengths don't match." unless ($n == $#{$y});

	$y2[0] = ($y->[1] - $y->[0])/($x->[1] - $x->[0]);
	$y2[$n] = ($y->[$n] - $y->[$n-1])/($x->[$n] - $x->[$n-1]);

	for my $i (1 .. $n-1)
	{
		$y2[$i] = ($y->[$i+1] - $y->[$i-1])/($x->[$i+1] - $x->[$i-1]);
	}

	return @y2;
}

=head3 seconddx()

    @d2ydx2 = seconddx(\@x, \@y);

or

    @d2ydx2 = seconddx(\@x, \@y, $yp0, $ypn);

Take references to two arrays containing the x and y ordinates of the
data and return an array of approximate second derivatives at the given x ordinates.

You may optionally give values to use as the first derivatives at the
start and end points of the data. If you don't, first derivative values
will be calculated from the arrays for you.

=cut

sub seconddx
{
	my($x, $y, $yp1, $ypn) = @_;
	my(@y2, @u);
	my($qn, $un);
	my $n = $#{$x};

	croak "X and Y array lengths don't match." unless ($n == $#{$y});

	if (defined $yp1)
	{
		$y2[0] = -0.5;
		$u[0] = (3/($x->[1] - $x->[0])) *
			(($y->[1] - $y->[0])/($x->[1] - $x->[0]) - $yp1);
	}
	else
	{
		$y2[0] = 0;
		$u[0] = 0;
	}

	for my $i (1 .. $n-1)
	{
		my $sig = ($x->[$i] - $x->[$i-1])/($x->[$i+1] - $x->[$i-1]);
		my $p = $sig * $y2[$i-1] + 2.0;

		$y2[$i] = ($sig - 1.0)/$p;
		$u[$i] = (6.0 * (
			($y->[$i+1] - $y->[$i])/($x->[$i+1] - $x->[$i]) -
			($y->[$i] - $y->[$i-1])/($x->[$i] - $x->[$i-1]))/
			($x->[$i+1] - $x->[$i-1]) - $sig * $u[$i-1])/$p;
	}

	if (defined $ypn)
	{
		$qn = 0.5;
		$un = (3.0/($x->[$n]-$x->[$n-1])) *
			($ypn - ($y->[$n] - $y->[$n-1])/($x->[$n] - $x->[$n-1]));
	}
	else
	{
		$qn = 0;
		$un = 0;
	}

	$y2[$n] = ($un - $qn * $u[$n-1])/($qn * $y2[$n-1] + 1.0);

	for my $i (reverse 0 .. $n-1)
	{
		$y2[$i] = $y2[$i] * $y2[$i+1] + $u[$i];
	}

	return @y2;
}

=head3 Derivative1()

A synonym for forwarddiff().

=head3 Derivative2()

A synonym for seconddx().

=cut

#
# Alias Derivative1() to centraldiff(), and Derivative2() to
# seconddx(), preserving the old names.
#
*Derivative1 = \&centraldiff;
*Derivative2 = \&seconddx;

=head1 COMPARISON OF FUNCTIONS

It is best to compare with a function for which we have exact derivatives. We'll use
the polynomial C<2*x**4 - 7*x**3 - 2*x**2 - x + 1>:

    use Math::Utils qw(:polynomial);
    use Math::Derivative qw(:all);

    #
    # The polynomial, its derivative, and its second derivative.
    #
    my @coef = (1, -1, -2, -7, 2);
    my @dx_coef = @{ pl_derivative(\@coef) };
    my @d2x_coef = @{ pl_derivative(\@dx_coef) };

    #
    # The X values in steps of 0.10,
    # the Y values from the polynomial and its derivatives.
    #
    my @xvals = map{$_ / 10} (20 .. 39);

    my @yvals = pl_evaluate(\@coef, \@xvals);
    my @dx_yvals = pl_evaluate(\@dx_coef, \@xvals);

    #
    # Get the first derivative approximations,
    # and compare them.
    #
    my @cen_dxdy = centraldiff(\@xvals, \@yvals);
    my @for_dxdy = forwarddiff(\@xvals, \@yvals);

    print " X      Y            dy       center   forward\n";
    for my $j (0 .. $#xvals)
    {
        printf("%5.2f, %5.2f:     %5.2f   %5.2f   %5.2f\n",
                $xvals[$j],
                $yvals[$j],
                $dx_yvals[$j],
                $cen_dxdy[$j],
                $for_dxdy[$j]);
    }

In this example C<centraldiff()> has better approximations
than C<forwarddiff()>.

    #
    # Now compare the values from the second derivative
    # with the second derivative approximations.
    #
    my @d2x_yvals = pl_evaluate(\@d2x_coef, \@xvals);
    my @d2xdy2 = seconddx(\@xvals, \@yvals);


=head1 REFERENCES

L<http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/central-differences/>

L<http://www.robots.ox.ac.uk/~sjrob/Teaching/EngComp/ecl6.pdf>

=head1 AUTHOR

John A.R. Williams B<J.A.R.Williams@aston.ac.uk>

John M. Gamble B<jgamble@cpan.org> (current maintainer)

=cut

1;
