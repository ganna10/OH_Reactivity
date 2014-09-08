#!/usr/bin/env perl
# Allocate NO3 reactivity to emitted species
# Version 0: Jane Coates 4/9/2014

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

my $reactant = 'NO3';
my $cair = $mecca->cair;
my $consumers = $kpp->consuming($reactant);
die "No reactions found for $reactant\n" if (@$consumers == 0) ;

my %total_reactivity;
foreach my $reaction (@$consumers) {
    my $reactants = $kpp->reactants($reaction);
    next if (scalar @$reactants == 1);
    my ($other_reactant) = grep { $_ ne $reactant } @$reactants;
    next if $other_reactant eq 'hv';
    my $rnum = $kpp->reaction_number($reaction);
    my $rconst = $mecca->rconst($rnum);
    next if (isbad($rconst(-1)));
    my $other_reactant_conc = $mecca->tracer($other_reactant) * $cair;
    my $reactivity = $rconst * $other_reactant_conc;
    next if ($reactivity->sum == 0);
    my ($number, $parent) = split /_/, $reaction;
    my $key_string = $kpp->reaction_string($reaction);
    if ($key_string =~ /_/) {
        $total_reactivity{$parent} += $reactivity(1:$ntime-2);
    } elsif ($key_string =~ /NO3\s\+\sOH/) {
        $total_reactivity{'OH'} += $reactivity(1:$ntime-2);
    } elsif ($key_string =~ /NO2\s\+\sNO3/) {
        $total_reactivity{'NO2'} += $reactivity(1:$ntime-2);
    } elsif ($key_string =~ /NO\s\+\sNO3/) {
        $total_reactivity{'NO'} += $reactivity(1:$ntime-2);
    } elsif ($key_string =~ /HO2\s\+\sNO3/) {
        $total_reactivity{'HO2'} += $reactivity(1:$ntime-2);
    } else {
        $total_reactivity{$key_string} += $reactivity(1:$ntime-2); 
    }
}

my $others_max = 2e-2;
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
    next if ($_ eq 'Others' or $_ eq 'NO' or $_ eq 'NO2' or $_ eq 'OH' or $_ eq 'HO2') ;
    push @final_sorted_data, { $_ => $total_reactivity{$_} };
} 
push @final_sorted_data, { 'Others' => $total_reactivity{'Others'} } if (defined $total_reactivity{'Others'}); 
unshift @final_sorted_data, { 'HO2' => $total_reactivity{'HO2'} } if (defined $total_reactivity{'HO2'}); 
unshift @final_sorted_data, { 'OH' => $total_reactivity{'OH'} } if (defined $total_reactivity{'OH'}); 

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
$times /= 3600;
my @time_axis = map { $_ } $times->dog; 
my @time_blocks;
foreach my $time (@time_axis) {#map to day and night
    if ($time <= 12) {
        push @time_blocks, "Day 1";
    } elsif ($time > 12 and $time <= 24) {
        push @time_blocks, "Night 1";
    } elsif ($time > 24 and $time <= 36) {
        push @time_blocks, "Day 2";
    } elsif ($time > 36 and $time <= 48) {
        push @time_blocks, "Night 2";
    } elsif ($time > 48 and $time <= 60) {
        push @time_blocks, "Day 3",
    } elsif ($time > 60 and $time <= 72) {
        push @time_blocks, "Night 3";
    } elsif ($time > 72 and $time <= 84) {
        push @time_blocks, "Day 4";
    } elsif ($time > 84 and $time <= 96) {
        push @time_blocks, "Night 4";
    } elsif ($time > 96 and $time <= 108) {
        push @time_blocks, "Day 5";
    } elsif ($time > 108 and $time <= 120) {
        push @time_blocks, "Night 5";
    } elsif ($time > 120 and $time <= 132) {
        push @time_blocks, "Day 6";
    } elsif ($time > 132 and $time <= 144) {
        push @time_blocks, "Night 6";
    } elsif ($time > 144 and $time <= 156) {
        push @time_blocks, "Day 7";
    } else {
        push @time_blocks, "Night 7";
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(reshape2) `,
        q` library(plyr) `,
        q` library(Cairo) `,
        q` library(grid) `,
);

$R->set('Time', [@time_blocks]);
$R->run(q` data = data.frame(Time) `);
foreach my $ref (@plot_data) {
    foreach my $process (sort keys %$ref) {
        $R->set('Process', $process);
        $R->set('Reactivity', [@{$ref->{$process}}]);
        $R->run(q` data[Process] = Reactivity `);
    }
}

$R->run(q` data = ddply(data, .(Time), colwise(sum)) `, #arrange data
        q` data = melt(data, id.vars = c("Time"), variable.name = "Process", value.name = "Reactivity") `,
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
                            "HO2" = "#2c9daf" ,
                            "NC5H12" = "#f9c500" ,
                            "C5H8" = "#8c1531" ,
                            "C2H6" = "#86b650" ,
                            "MXYL" = "#ef6638" ,
                            "C3H6" = "#0352cb" ,
                            "NO" = "#dc3522" ,
                            "NC6H14" = "#9bb18d" ,
                            "IC4H10" = "#a67c52" ,
                            "PXYL" = "#77aecc" ,
                            "BENZENE" = "#f7c56c" ,
                            "OXYL" = "#4c9383" ,
                            "EBENZ" = "#ba8b01" ,
                            "O" = "#8d1435" ,
                            "OH" = "#6db875" ,
                            "BUT1ENE" = "#77afc8" ,
                            "MEPROPENE" = "#cc6329" ,
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
                            "BUT1ENE" = "1-Butene",
                            "MEPROPENE" = "2-Methyl Propene",
                            "EBENZ" = "Ethylbenzene") `,
);

$R->run(q` plot = ggplot(data = data, aes(x = Time, y = Reactivity, fill = Process)) `, #plot data
        q` plot = plot + geom_bar(stat = "identity", width = 0.7) `,
        q` plot = plot + scale_x_discrete(limits = c("Day 1", "Night 1", "Day 2", "Night 2", "Day 3", "Night 3", "Day 4", "Night 4", "Day 5", "Night 5", "Day 6", "Night 6", "Day 7", "Night 7")) `,
        q` plot = plot + scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, 0.5)) `,
        q` plot = plot + ggtitle("NO3 Reactivity") `,
        q` plot = plot + ylab(expression(paste("Reactivity (", s^-1, ")"))) `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + theme(plot.title = element_text(size = 22, face = "bold")) `,
        q` plot = plot + theme(axis.title.x = element_blank()) `,
        q` plot = plot + theme(axis.title.y = element_text(size = 20, face = "bold")) `,
        q` plot = plot + theme(axis.text.x = element_text(size = 18)) `,
        q` plot = plot + theme(axis.text.y = element_text(size = 18)) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(legend.text = element_text(size = 18)) `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + theme(legend.key.size = unit(1.3, "cm")) `,
        q` plot = plot + scale_fill_manual(values = my.colours, labels = my.names) `,
);

$R->run(q` CairoPDF(file = "NO3_reactivity_allocation_without_NOx.pdf", width = 20, height = 14) `, #save plot to file
        q` print(plot) `,
        q` dev.off() `,
);

#my $p = $R->run(q` print(data) `);
#print "$p\n";

$R->stop;
