#!/usr/bin/env perl
# Allocate reactivity to Pentane's degradation products
# Version 0: Jane Coates 8/9/2014

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use PDL::Bad;
use Statistics::R;

my $VOC = "NC5H12";
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
    my ($number, $parent) = split /_/, $reaction;
    next unless (defined $parent and $parent eq $VOC);
    my $reactants = $kpp->reactants($reaction);
    my ($other_reactant) = grep { $_ ne $reactant } @$reactants;
    next if $other_reactant eq 'hv';
    my $rnum = $kpp->reaction_number($reaction);
    my $rconst = $mecca->rconst($rnum);
    next if (isbad($rconst(-1)));
    my $other_reactant_conc = $mecca->tracer($other_reactant) * $cair;
    my $reactivity = $rconst * $other_reactant_conc;
    next if ($reactivity->sum == 0);
    $other_reactant =~ s/_(.*?)$//g;
    $total_reactivity{$other_reactant} += $reactivity(1:$ntime-2);
}
my $others_max = 1;
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
    next if ($_ eq 'Others') ;
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
    foreach my $reactant (sort keys %$ref) {
        $R->set('Reactant', $reactant);
        $R->set('Reactivity', [@{$ref->{$reactant}}]);
        $R->run(q` data[Reactant] = Reactivity `);
    }
}

$R->run(q` data = ddply(data, .(Time), colwise(sum)) `, #arrange data
        q` data = melt(data, id.vars = c("Time"), variable.name = "Reactant", value.name = "Reactivity") `,
        q` data$Reactant = factor(data$Reactant, levels = rev(levels(data$Reactant))) `,
);

$R->run(q` my.colours = c(  "CO" = "#2b9eb3" ,
                            "HCHO" = "#b569b3" ,
                            "CH3CHO" = "#1b695b" ,
                            "CO2C5OH" = "#ae4901" ,
                            "CO2C3CHO" = "#ef6638" ,
                            "DIEK" = "#e7e85e" ,
                            "MPRK" = "#0e5c28" ,
                            "CO2C4CHO" = "#f3aa7f" ,
                            "C2H5CHO" = "#898989" ,
                            "NC5H12" = "#f9c500" ,
                            "Others" = "#58691b" ) `,
        q` my.names = c(    "NC5H12" = "Pentane" ) `,
);

$R->set('title', "OH Reactivity during Pentane Degradation");
$R->run(q` plot = ggplot(data = data, aes(x = Time, y = Reactivity, fill = Reactant)) `, #plot data
        q` plot = plot + geom_bar(stat = "identity", width = 0.7) `,
        q` plot = plot + scale_x_discrete(limits = c("Day 1", "Night 1", "Day 2", "Night 2", "Day 3", "Night 3", "Day 4", "Night 4", "Day 5", "Night 5", "Day 6", "Night 6", "Day 7", "Night 7")) `,
        q` plot = plot + scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, 2)) `,
        q` plot = plot + ggtitle(title) `,
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

$R->set('filename', "${VOC}_OH_allocated_reactivity.pdf");
$R->run(q` CairoPDF(file = filename, width = 20, height = 14) `, #save plot to file
        q` print(plot) `,
        q` dev.off() `,
);

#my $p = $R->run(q` print(data) `);
#print "$p\n";

$R->stop;
