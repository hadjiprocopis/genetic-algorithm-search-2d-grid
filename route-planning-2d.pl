#!/usr/bin/env perl

####################################################################
# demonstrate the use of genetic algorithms in searching a 2D grid
# heavily based on Algorithm::Evolutionary by J. J. Merelo-Guervos
# author: bliako (bliako//at//cpan.org)
# github: https://github.com/hadjiprocopis/genetic-algorithm-search-2d-grid
# date: 13/03/2020
# basic usage: $0 --grid-width 32 --grid-height 32
####################################################################

use strict;
use warnings;

use Getopt::Long;

use Algorithm::Evolutionary::Experiment;
use Algorithm::Evolutionary::Op::Easy;
use Algorithm::Evolutionary::Op::Bitflip;
use Algorithm::Evolutionary::Op::Crossover;

##############################
# begin the Grid etc. packages
##############################
{
package Path;
use strict;
use warnings;

use overload
	'""' => 'toString';

sub new {
	my ($class, $grid) = @_;
	my $W = $grid->{'w'};
	my $self = {
		'w' => $W,
		'start' => $grid->{'start'},
		'end' => $grid->{'end'},
		# the first and last part of the path are on the same col as the start/end square of the spec
		# that's ok, it gives extra flexibility.
		'Y' => [ (0)x$W ],
		# however, we need 1 extra slide-bit from the move from last path square to the end-square of the spec
		# the slide-bit, i.e. if 1 then do not collect score from passing squares
		'slide' => [ (0)x($W+1) ],
		'score' => undef,
		'distance' => undef,
	};
	bless $self => $class;
	print "Path::new() : called for width ".$self->{'w'}.".\n";
	return $self;
}
sub Y {
	my ($self, $x, $value) = @_;
	$self->{'Y'}->[$x] = $value if defined $value;
	return $self->{'Y'}->[$x];
}
sub slide {
	my ($self, $x, $value) = @_;
	$self->{'slide'}->[$x] = $value if defined $value;
	return $self->{'slide'}->[$x];
}
sub toString {
	my $self = $_[0];
	my $W = $self->{'w'};

	my ($y);
	my $y0 = $self->{'start'}->[1];
	my $x0 = 0;
	my $ret = "[";
	for(my $x=0;$x<$W;$x++){
		$y = $self->{'Y'}->[$x];
		$ret .= 
			 ($self->{'slide'}->[$x]==1?'S':'N')
			.":($x0,$y0)->($x,$y)"
			.','
		;
		$x0 = $x; $y0 = $y;
	}
	return $ret
		. ($self->{'slide'}->[$W]==1?'S':'N')
		. ":($x0,$y0)->(".($W-1).",".$self->{'end'}->[1].")"
		.']'
}
} # end package Path
{
package Cell;
use strict;
use warnings;

use overload
	'""' => 'toString';

our $VISITED_MARKS = {0=>' ', 1=>'+', 2=>'*'};

sub new {
	my ($class, $x, $y, $score) = @_;
	my $self = {
		'x' => $x,
		'y' => $y,
		'score' => $score,
		# zero means we have not visited, 1 means we visited but sliding (did not collect score)
		# 2 means visited and collected score (no slide)
		'visited' => 0,
	};
	bless $self => $class;
	#print "Cell::new() : called: for ($x, $y).\n";
	return $self;
}
sub visit { $_[0]->{'visited'} = $_[1] }
sub visited { return $_[0]->{'visited'} }
sub visited_mark { return $Cell::VISITED_MARKS->{$_[0]->visited()} }
sub score {
	my $self = $_[0];
	my $m = $_[1];
	if( defined $m ){ $self->{'score'} = $m; }
	return $_[0]->{'score'}
}
sub toString {
	my $self = $_[0];
	# a visited of 0 means we did not visit and no mark is there (mark is ' ')
	# visited of 1 means we visited and collected score, mark is '+'
	# visited of 2 means we visited but sliding (no score collected), mark is '*'
	return
		$self->{'x'}.",".$self->{'y'}
		. '='
		. $self->{'score'}
		. $self->visited_mark()
}
} # end package Cell

{
package Grid;
use strict;
use warnings;

use overload
	'""' => 'toString';

sub new {
	my ($class, $w, $h, $start, $end, $maxscore) = @_;
	my $self = {
		'w' => $w,
		'h' => $h,
		'start' => $start,
		'end' => $end,
		'cells' => undef,
		'path' => undef,
		'cache' => {},
		'cache-usage' => {'read' => 0, 'write' => 0},
		'maxscore' => $maxscore // 10
	};
	bless $self => $class;
	print "Grid::new() : called for dims: ($w x $h), start: @{$start}, end: @{$end}.\n";
	$self->populate_randomly($self->{'maxscore'});
	return $self;
}
sub populate_randomly {
	# put random scores to each cell for testing purposes
	my ($self, $maxscore) = @_;
	my $x = 0;
	$self->{'cells'} = [];
	$self->{'cache'} = {};
	for(my $i=0;$i<$self->{'w'};$i++){
		my @acol = ();
		for(my $j=0;$j<$self->{'h'};$j++){
			$acol[$j] = Cell->new($i, $j, $maxscore - int(rand(2*$maxscore+1)));
		}
		$self->{'cells'}->[$i] = \@acol;
	}
}
sub populate_test_path_randomly {
	# create a chain of cells (a path) with high scores which we will use to test our algorithm
	my $self = $_[0];
	my $pathscore = $_[1] // $self->{'maxscore'}*2;

	my $C = $self->{'cells'};
	my $w = $self->{'w'};
	my $h = $self->{'h'};
	for(my $i=0;$i<$w;$i++){
		$C->[$i]->[int(rand($h))]->score($pathscore);
	}
}
sub calculate_path {
	my ($self, $path, $setpath) = @_;
	#$self->unpath();
	my $C = $self->{'cells'};

	#print "Grid::calculate_path() : called.\n";
	my ($acell, $y, $j, $slide, $A, $B, $ascore, $adistance, $X);

	# create these shortcuts so we don't seek them inside the loop
	my $PW = $path->{'w'};
	my $PY = $path->{'Y'};
	my $PS = $path->{'slide'};
	my $cache = $self->{'cache'};

	# start from the starting square, do not collect score because this will happen inside the loop
	my $y0 = $self->{'start'}->[1];
	my $score = 0;
	my $distance = 0;
	for(my $x=0;$x<=$PW;$x++){
		$slide = $PS->[$x];
		if( $x==$PW ){
			$y = $self->{'end'}->[1];
			if( $y == $y0 ){ last; }
			$X = $x-1;
		} else { $y = $PY->[$x]; $X = $x; }

		my $segment = join(':', $X, $y0, $y, $slide);
		#print "entering ($x,$y0)->($x->$y) slide=$slide\n";
		if( $setpath==0 && exists $cache->{$segment} ){
			($ascore, $adistance, $y0) = @{$cache->{$segment}};
			#print "Grid::calculate_path() : using cache for '$segment' : $ascore, $adistance, $y0\n";
			$self->{'cache-usage'}->{'read'}++;
			$score += $ascore;
			$distance += $adistance;
		} else {
			_calcpath($X, $y, $y0, $slide, $C, $setpath, \$ascore, \$adistance);
			$y0 = $y;
			$cache->{$segment} = [$ascore, $adistance, $y0];
			#print "Grid::calculate_path() : set cache for '$segment' : $ascore, $adistance, $y0\n";
			$self->{'cache-usage'}->{'write'}++;
			$score += $ascore;
			$distance += $adistance;
		}
		#print "ascore=$ascore, total=$score\n";
		#print $self;
	}
	# score of end square is already collected

	$path->{'score'} = $score;
	$path->{'distance'} = $distance;
	return [$distance, $score];
}
sub _calcpath {
	my ($X, $y, $y0, $slide, $C, $setpath, $ascoreref, $adistanceref) = @_;
	my ($A, $B);
	my $DEBUG=0;
	if( $y0 < $y ){ $A = $y0; $B = $y; } else { $A = $y; $B = $y0; }

	# collect the score of the starting square even if we slide but only the START square
	my $acell = $C->[$X]->[$A];
	$acell->visit(1) if $setpath==1;
	my $ascore = $acell->{'score'};
	print "collecting1 ($X,$A) = ".$acell->{'score'}."\n" if $DEBUG;

	for(my $j=$A+1;$j<$B;$j++){
		$acell = $C->[$X]->[$j];
		if( $slide == 0 ){
			$ascore += $acell->{'score'};
			print "collecting2 ($X,$j) = ".$acell->{'score'}."\n" if $DEBUG;
		}# else { print "visiting2 ($X,$j) = ".$acell->{'score'}."\n"; }
		# visited of 1 means we visited and collected score, mark is '+'
		# visited of 2 means we visited but sliding (no score collected), mark is '*'
		$acell->visit($slide+1) if $setpath==1; # 1 means not collected score, 2 means score collected
	}
	# collect the score of the end square even if we slide but only the END square
	if( $y0 != $y ){
		$acell = $C->[$X]->[$B];
		print "collecting3 ($X,$B) = ".$acell->{'score'}."\n" if $DEBUG;
		$ascore += $acell->{'score'};
		$acell->visit(1) if $setpath==1;
	}
	$$ascoreref = $ascore;
	$$adistanceref = $B-$A+1;
}

sub unpath {
	my $self = $_[0];
	my $C = $self->{'cells'};
	my $w = $self->{'w'};
	my $h = $self->{'h'};
	my ($i, $j);
	for($i=0;$i<$w;$i++){
		for($j=0;$j<$h;$j++){
			$C->[$i]->[$j]->visit(0);
		}
	}
}
sub toString {
	my $self = $_[0];
	my $ret = "";

	my $C = $self->{'cells'};
	my $w = $self->{'w'};
	my $h = $self->{'h'};

	my ($i, $j, $acell);
	$ret .= sprintf("%4s|", ""); for($i=0;$i<$w;$i++){ $ret .= sprintf("%4s|", "$i"); } $ret .= "\n";
	for($i=0;$i<($w+1);$i++){ for($j=0;$j<5;$j++){ $ret .= "-" } } $ret .= "\n";
	for($j=0;$j<$h;$j++){
		$ret .= sprintf("%4s|", $j);
		for($i=0;$i<$w;$i++){
			$acell = $C->[$i]->[$j];
			$ret .= sprintf("%s%3s|", $acell->visited_mark(), "".$C->[$i]->[$j]->score())
		}
		$ret .= "\n";
	}
	$ret .= "cache usage: read: ".$self->{'cache-usage'}->{'read'}.", write: ".$self->{'cache-usage'}->{'write'}."\n";
	return $ret
}
} # end package Grid
#############################
# end the Grid etc. packages
#############################

########################
# begin the main program
########################

# grid size
my $W = 16;
my $H = 16; # <<< Height must be multiple of 2
my $seed_create = 123; # seed to create a test path, fix it so that all the same
my $seed_search = time; # seed for searching with GA, can fix it with options or leave it dangling
my $Y0 = -1;
my $YN = -1;
# search for those maximum iterations or when fitness changes less by the amount specified
my $MAXITERS=100;
my $BREAK_WHEN_FITNESS_DOES_NOT_CHANGE = 10E-03;
my $STOPEARLY = 1;
if( ! Getopt::Long::GetOptions(
	"grid-width|w=i" => \$W,
	"grid-height|h=i" => \$H,
	"seed-create=i" => \$seed_create,
	"seed-search=i" => \$seed_search,
	"start-at-y=i" => \$Y0,
	"end-at-y=i" => \$YN,
	"max-iters=i" => \$MAXITERS,
	"stop-early!" => \$STOPEARLY,
	"help|h" => sub {
		print "Usage : $0 [--grid-width W] [--grid-height H] [--seed-create S] [--seed-search S] [--start-at-y Y] [--end-at-y Y] [--max-iters M] [--(no-)stop-early]\n";
		exit(0);
	},
) ){ die "error in command line arguments.\n"; }

my $bits_per_height_gene = int(log($H)/log(2));
die "Grid height must be a power of 2 but it was $H != 2^${bits_per_height_gene}" unless 2**(int(log($H)/log(2))) == $H;

# Y coord of the starting and ending squares
$Y0 = int($H/2) unless $Y0 >= 0; 
$YN = int($H/4) unless $YN >= 0;

if( ! $STOPEARLY ){ $BREAK_WHEN_FITNESS_DOES_NOT_CHANGE = -1 }

# score in each square ranges between +- MAXSCORE
my $MAXSCORE=10;

print "$0 : creating the grid and a test-path (seed $seed_create) ...\n";

srand($seed_create);
# nothing to change below
my $tstarted = time;
my $G = Grid->new($W, $H, [0,$Y0], [$W-1, $YN], $MAXSCORE);
$G->populate_test_path_randomly($MAXSCORE*2+1);
my $myPath = Path->new($G);

print "$0 : searching (seed $seed_search) ...\n";
srand($seed_search);

###### make genetic
# our chromosome consists of bits
# the first bit is whether we slide from start square (which is fixed and not part of the chromosome)
# to the next square, call it B
# the second bit is whether we slide from B to next square, C
# the next N bits represent the y-coordinate of the B square
# the max y-coordinate is $H, so N = log_2($H)
# this continues up to the last square of the path
# then we have one extra bit representing the slide from last square of path to our
# ending square (which is not part of the path)
# So the path consists of $W-1 y-coordinates plus their slide-bit (of each)
# plus 2 extra slide-bits
my $chromosomeSizeInBits = $W*(1 + $bits_per_height_gene) + 1;

print "$0 : number of bits in the chromosome: $chromosomeSizeInBits\n";

sub calculate_path_fitness {
	my $individual = $_[0];
	my $chromosome = $individual->Chrom();
	chromosome2genes($chromosome, $myPath); # fills $directions and $Y
	my ($distance, $score) = @{ $G->calculate_path($myPath, 0) };
	return $score / $distance;
}
sub set_path_to_grid {
	my $individual = $_[0];
	my $chromosome = $individual->Chrom();
	chromosome2genes($chromosome, $myPath); # fills $directions and $Y
	my ($distance, $score) = @{ $G->calculate_path($myPath, 1) };
	return $score / $distance;
}

sub chromosome2genes {
	# convert a chromosome which consists of genes which consist of bits(alleles)
	# into a set of numbers to be applied to our problem.
	# that is: 1bit for slide, and as many bits required for the y-coordinate of the target square
	# plus 1bit for the ending slide

	# chromosome bit string containing all genes as 10101	
	my ($achromosome, $apath) = @_;
	#print "chromosome2genes() : $achromosome\n";
	my $x = 0;
	while( $achromosome =~ /([01])([01]{$bits_per_height_gene})/g ){
		# it means we move to Y with this slide-bit
		$apath->{'slide'}->[$x] = $1;
		$apath->{'Y'}->[$x] = bin2dec($2);
		$x++;
	}
	# now we have one extra slide-bit left at the end for moving to the target square given by the spec
	$achromosome =~ /([01])$/;
	$apath->{'slide'}->[$x] = $1;
}
sub bin2dec {
	# MSB is the last one, e.g. 011 = 6
	my $in = $_[0];
	my $g = 0;
	my $j = 1;
	map { $g += $_*$j; $j*=2; } split(//, $in);
	return $g;
}

my $m = Algorithm::Evolutionary::Op::Bitflip->new(3); # flip this number of bits randomly
my $c = Algorithm::Evolutionary::Op::Crossover->new(2, 4); # crossover with 2 points
my $ez = new Algorithm::Evolutionary::Op::Easy \&calculate_path_fitness, 0.8, [$m,$c];
#my $ez = new Algorithm::Evolutionary::Op::CanonicalGA \&fitness, 0.8, [$m,$c];
my $popSize = 500; # population size, each individual in this pop has a chromosome which consists of 2 genes
my $chromosomeType = 'BitString'; # the chromosome is a sequence of bits as a string
my $e = new Algorithm::Evolutionary::Experiment $popSize, $chromosomeType, $chromosomeSizeInBits, $ez;

my $populationRef;
my $previous_fitness = 0;
my $current_fitness = 0;
my ($best_solution, $best_fitness);
my $iter = 0;
my $stale=0;
while( (++$iter<$MAXITERS) && ($stale<10) ){
	# create a new generation of solutions, this is one iteration in the genetic algorithm:
	$populationRef = $e->go();
	# the first in the population of solutions is the best for this generation:
	$best_solution = $populationRef->[0];
	$best_fitness = $best_solution->Fitness();
	print "$iter / $MAXITERS) : fitness: ($previous_fitness -> $current_fitness ->) $best_fitness\n";
	if( ($BREAK_WHEN_FITNESS_DOES_NOT_CHANGE>0) && (($current_fitness - $previous_fitness) < $BREAK_WHEN_FITNESS_DOES_NOT_CHANGE) ){ $stale++ } else { $stale = 0; }
	$previous_fitness = $current_fitness;
	$current_fitness = $best_fitness;
}
set_path_to_grid($best_solution);
print "Best solution at iteration $iter:\n$G\n";
print "$myPath\n";
print "score=".$myPath->{'score'}.", distance=".$myPath->{'distance'}.", fitness=$current_fitness, iteration=$iter/$MAXITERS\n";
print "time taken: ".(time-$tstarted)." seconds.\n";
print "$0 : done.\n";
