#!/usr/bin/env perl
# Allocate reactivity to emitted species
# Version 0: Jane Coates 30/8/2014

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use PDL::Bad;
use Statistics::R;

my $model_run = "/local/home/coates/MECCA/MCM_3.2_tagged";
my $boxmodel = "$model_run/boxmodel";
my $mecca = MECCA->new($boxmodel);
my $eqn = "$model_run/gas.eqn";
my $kpp = KPP->new($eqn);

my $reactant = 'OH';
my @inorganic = qw(OH O3 H2O2 NO NO2 HO2 HNO3 HONO HO2NO2 H2 NO3);
my %inorganic; 
$inorganic{$_} = 1 foreach (@inorganic);
my $cair = $mecca->cair;
my $consumers = $kpp->consuming($reactant);
die "No reactions found for $reactant\n" if (@$consumers == 0) ;

my $total_reactivity;
foreach my $reaction (@$consumers) {
    my $reactants = $kpp->reactants($reaction);
    my ($other_reactant) = grep { $_ ne $reactant } @$reactants;
    next if $other_reactant eq 'hv';
    my $rnum = $kpp->reaction_number($reaction);
    my $rconst = $mecca->rconst($rnum);
    next if isbad($rconst(-1));
    my $other_reactant_conc = $mecca->tracer($other_reactant) * $cair;
    my $reactivity = $rconst * $other_reactant_conc;
    print "$reactivity\n";
    #$total_reactivity += $reactivity;
}
