#!/usr/bin/env perl
# Allocate HCHO OH reactivity to emitted species
# Version 0: Jane Coates 20/9/2014

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use PDL::Bad;
use Statistics::R;

my $model_run = "/local/home/coates/Documents/OH_Reactivity";
my $boxmodel = "$model_run/boxmodel";
my $mecca = MECCA->new($boxmodel);
my $eqn = "$model_run/gas.eqn";
my $kpp = KPP->new($eqn);
my $ntime = $mecca->time->nelem;

my $reactant = 'OH';
my $cair = $mecca->cair;
my $consumers = $kpp->consuming($reactant);
die "No reactions found for $reactant\n" if (@$consumers == 0) ;

my %total_reactivity;
foreach my $reaction (@$consumers) {
    my $reactants = $kpp->reactants($reaction);
    my ($other_reactant) = grep { $_ ne $reactant } @$reactants;
    next unless ($other_reactant =~ /\bHCHO/);
    my $rnum = $kpp->reaction_number($reaction);
    my $rconst = $mecca->rconst($rnum);
    next if (isbad($rconst(-1)));
    my $other_reactant_conc = $mecca->tracer($other_reactant) * $cair;
    my $reactivity = $rconst * $other_reactant_conc;
    next if ($reactivity->sum == 0);
    my ($number, $parent) = split /_/, $reaction;
    $total_reactivity{$parent} += $reactivity(1:$ntime-2);
}

my $others_max = 2.5;
foreach my $process (keys %total_reactivity) {
    if ($total_reactivity{$process}->sum <= $others_max) {
        $total_reactivity{'Others'} += $total_reactivity{$process};
        delete $total_reactivity{$process};
    }
}

my $sort_function = sub { $_[0]->sum };
my @sorted_data = sort { &$sort_function($total_reactivity{$b}) <=> &$sort_function($total_reactivity{$a}) } keys %total_reactivity;

my @final_sorted_data;
foreach (@sorted_data) { 
    next if ($_ eq 'Others');
    push @final_sorted_data, { $_ => $total_reactivity{$_} };
} 
push @final_sorted_data, { 'Others' => $total_reactivity{'Others'} } if (defined $total_reactivity{'Others'}); 

my @plot_data;
foreach my $ref (@final_sorted_data) {#extract reaction and rates for each plot
    foreach my $item (keys %$ref) {
        my @rate_array = map { $_ } $ref->{$item}->dog;
        push @plot_data, { $item => \@rate_array };
    }
} 

my $times = $mecca->time; #time axis
$times -= $times->at(0);
$times = $times(1:$ntime-2);
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
foreach my $ref (@plot_data) {
    foreach my $process (sort keys %$ref) {
        $R->set('Process', $process);
        $R->set('Reactivity', [@{$ref->{$process}}]);
        $R->run(q` data[Process] = Reactivity `);
    }
}

$R->run(q` data = melt(data, id.vars = c("Time"), variable.name = "Process", value.name = "Reactivity") `,
        q` data$Process = factor(data$Process, levels = rev(levels(data$Process))) `,
);

$R->run(q` my.colours = c(  "CO" = "#2b9eb3" ,
                            "NO2" = "#b569b3" ,
                            "CH4" = "#1b695b" ,
                            "IC5H12" = "#ae4901" ,
                            "C3H8" = "#e7e85e" ,
                            "TOLUENE" = "#0e5c28" ,
                            "NC4H10" = "#f3aa7f" ,
                            "C2H4" = "#898989" ,
                            "O3" = "#1c3e3d" ,
                            "NC5H12" = "#f9c500" ,
                            "C5H8" = "#8c1531" ,
                            "C2H6" = "#86b650" ,
                            "MXYL" = "#ef6638" ,
                            "C3H6" = "#0352cb" ,
                            "NO" = "#c9a415" ,
                            "NC6H14" = "#9bb18d" ,
                            "IC4H10" = "#a67c52" ,
                            "PXYL" = "#dc3522" ,
                            "BENZENE" = "#f7c56c" ,
                            "OXYL" = "#4c9383" ,
                            "EBENZ" = "#ba8b01" ,
                            "Others" = "#58691b" ) `,
        q` my.names = c(    "CH4" = "Methane" ,
                            "IC5H12" = "2-Methyl Butane" ,
                            "C3H8" = "Propane" ,
                            "TOLUENE" = "Toluene" ,
                            "NC4H10" = "Butane" ,
                            "C2H4" = "Ethene" ,
                            "NC5H12" = "Pentane" ,
                            "C5H8" = "Isoprene" ,
                            "C2H6" = "Ethane" ,
                            "MXYL" = "m-Xylene" ,
                            "C3H6" = "Propene" ,
                            "NC6H14" = "Hexane" ,
                            "IC4H10" = "2-Methyl Propane" ,
                            "PXYL" = "p-Xylene" ,
                            "BENZENE" = "Benzene" ,
                            "OXYL" = "o-Xylene" ,
                            "EBENZ" = "Ethylbenzene") `,
);

$R->run(q` plot = ggplot(data = data, aes(x = Time, y = Reactivity, fill = Process)) `, #plot data
        q` plot = plot + geom_area(position = "stack") `,
        q` plot = plot + geom_area(position = "stack", colour = "black", show_guide = FALSE) `,
        q` plot = plot + ggtitle("Formaldehyde OH Reactivity Distributed by VOC Source\n") `,
        q` plot = plot + ylab(expression(bold(paste("OH Reactivity (", s^-1, ")")))) `,
        q` plot = plot + xlab("Time (days)") `,
        q` plot = plot + scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, 1)) `,
        q` plot = plot + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + theme(plot.title = element_text(size = 32, face = "bold")) `,
        q` plot = plot + theme(axis.title.x = element_text(size = 25, face = "bold")) `,
        q` plot = plot + theme(axis.title.y = element_text(size = 25, face = "bold")) `,
        q` plot = plot + theme(axis.text.x = element_text(size = 20)) `,
        q` plot = plot + theme(axis.text.y = element_text(size = 20)) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(legend.text = element_text(size = 20)) `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + theme(legend.key.size = unit(1.3, "cm")) `,
        q` plot = plot + theme(panel.grid.major = element_blank()) `,
        q` plot = plot + theme(panel.grid.minor = element_blank()) `,
        q` plot = plot + theme(panel.border = element_blank()) `,
        q` plot = plot + theme(axis.line = element_line(colour = "black")) `,
        q` plot = plot + scale_fill_manual(values = my.colours, labels = my.names) `,
);

$R->run(q` CairoPDF(file = "HCHO_OH_reactivity_VOC_allocation.pdf", width = 20, height = 14) `, #save plot to file
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop;
