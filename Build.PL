use 5.006;
use strict;
use warnings FATAL => 'all';
use Module::Build;

my $builder = Module::Build->new(
	module_name => 'Math::Derivative',

	dist_abstract => q(Numeric 1st and 2nd order differentiation.),
	dist_author => [q{J. A. R. Williams},
			q{John M. Gamble <jgamble@cpan.org>}],
	dist_version => '1.01',

	requires => {
		perl => '5.10.1',
	},

	configure_requires => {
		'Module::Build' => 0,
	},

	build_requires => {
		'Test::More' => 0,
		'Math::Utils' => '1.10',
	},

	license => 'perl',
	create_license => 1,
	create_readme => 0,
	create_makefile_pl => 'traditional',
	dynamic_config => 0,
	release_status => 'stable',
	add_to_cleanup     => [ 'Math-Derivative-*' ],
	meta_merge => {
            resources => {
                repository => 'git://github.com/jgamble/Math-Derivative.git'
            },
	}
);

$builder->create_build_script();
