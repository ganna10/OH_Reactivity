#! /usr/bin/env perl
# plot total OH reactivity and OH reactivity just of emitted initial VOC OH reactivity
# Version 0: Jane Coates 20/9/2014
# Version 1: Jane Coates 1/10/2014 including specific inorganic section

use strict; 
use diagnostics;
use MECCA;
use KPP;
use PDL;
use PDL::NiceSlice;
use PDL::Bad;
use Statistics::R;

my $model_run = "/local/home/coates/Documents/OH_Reactivity";
my $boxmodel = "$model_run/boxmodel";
my $mecca = MECCA->new($boxmodel);
my $eqn = "$model_run/gas.eqn";
my $kpp = KPP->new($eqn);
my $NTIME = $mecca->time->nelem;
my @organic = qw( CO CH4 C2H6 C3H8 NC4H10 IC4H10 NC5H12 IC5H12 NC6H14 NC7H16 NC8H18 C2H4 C3H6 BUT1ENE MEPROPENE C5H8 BENZENE TOLUENE MXYL OXYL PXYL EBENZ );
my @inorganic = qw( NO NO2 O3  );

my $reactant = 'OH';
my $cair = $mecca->cair;
my $consumers = $kpp->consuming($reactant);
die "No reactions found for $reactant\n" if (@$consumers == 0) ;

my %total_reactivity;
foreach my $reaction (@$consumers) {
    my $reactants = $kpp->reactants($reaction);
    my ($other_reactant) = grep { $_ ne $reactant } @$reactants;
    next if $other_reactant eq 'hv';
    my $rnum = $kpp->reaction_number($reaction);
    my $rconst = $mecca->rconst($rnum);
    next if (isbad($rconst(-1)));
    my $other_reactant_conc = $mecca->tracer($other_reactant) * $cair;
    my $reactivity = $rconst * $other_reactant_conc;
    next if ($reactivity->sum == 0);
    my ($number, $parent) = split /_/, $reaction;
    $total_reactivity{'Total'} += $reactivity(1:$NTIME-2); 
    if ($other_reactant ~~ @organic) {
        $total_reactivity{'VOCs'} += $reactivity(1:$NTIME-2);
    } elsif ($other_reactant ~~ @inorganic) {
        $total_reactivity{'Inorganic'} += $reactivity(1:$NTIME-2);
    }
}

my @plot_array = (
    {'Inorganic'    => $total_reactivity{'Inorganic'}},
    { 'VOCs'        => $total_reactivity{'VOCs'} },
    { 'Total'       => $total_reactivity{'Total'} } ,
);

my $times = $mecca->time; #time axis
$times -= $times->at(0);
$times = $times(1:$NTIME-2);
$times /= 86400;
my @time_axis = map { $_ } $times->dog; 

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(reshape2) `,
        q` library(plyr) `,
        q` library(Cairo) `,
        q` library(grid) `,
);

$R->set('Time', [@time_axis]);
$R->run(q` data = data.frame(Time) `);
foreach my $ref (@plot_array) {
    foreach my $item (sort keys %$ref) {
        $R->set('run', $item);
        $R->set('reactivity', [map { $_ } $ref->{$item}->dog]);
        $R->run(q` data[run] = reactivity `);
    }
}

$R->run(q` data$Total = data$Total - data$VOCs - data$Inorganic `);
$R->run(q` data = melt(data, id.vars = c("Time"), variable.name = "Run", value.name = "OH.Reactivity") `);
$R->run(q` data$Run = factor(data$Run, levels = c("Total", "VOCs", "Inorganic")) `);
#my $p = $R->run(q` print(data) `);
#print $p, "\n";

$R->run(q` plot = ggplot(data, aes(x = Time, y = OH.Reactivity, fill = Run)) `,
        q` plot = plot + geom_area(position = "stack") `,
        q` plot = plot + geom_area(position = "stack", colour = "black", show_guide = FALSE) `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + xlab("Time (days)") `,
        q` plot = plot + ylab(expression(bold(paste("OH Reactivity (", s^-1, ")")))) `,
        q` plot = plot + ggtitle("Total OH Reactivity and Emitted VOC OH Reactivity") `,
        q` plot = plot + scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, 1)) `,
        q` plot = plot + theme(panel.grid.major = element_blank()) `,
        q` plot = plot + theme(panel.grid.minor = element_blank()) `,
        q` plot = plot + theme(panel.border = element_blank()) `,
        q` plot = plot + theme(axis.line = element_line(colour = "black")) `,
        q` plot = plot + theme(axis.title.y = element_text(size = 25)) `,
        q` plot = plot + theme(axis.title.x = element_text(size = 25, face = "bold")) `,
        q` plot = plot + theme(axis.text.y = element_text(size = 20)) `,
        q` plot = plot + theme(axis.text.x = element_text(size = 20)) `,
        q` plot = plot + theme(plot.title = element_text(size = 32, face = "bold")) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(legend.text = element_text(size = 25)) `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + theme(legend.position = c(0.99, 0.6)) `,
        q` plot = plot + theme(legend.justification = c(0.99, 0.6)) `,
        q` plot = plot + theme(legend.key.size = unit(1.3, "cm")) `,
        q` plot = plot + scale_fill_manual(values = c("#2b9eb3", "#f9c500", "#8c1531"), labels = c("Total OH Reactivity", "Emitted VOC OH Reactivity", "Inorganic OH Reactivity")) `,
);

$R->run(q` CairoPDF(file = "Total_vs_parents_OH_reactivity.pdf", width = 20, height = 14) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();
