#!/usr/bin/env perl
# Compare total Isoprene and toluene OH reactivity 
# Version 0: Jane Coates 8/9/2014

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
my @parents = qw( C5H8 TOLUENE );

my %plot_data;
foreach my $VOC (@parents) {
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
        $other_reactant =~ s/_(.*?)\b$//;
        $total_reactivity{$other_reactant} += $reactivity(1:$ntime-2);
    }

    my $others_max = 1.35;
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
    foreach my $ref (@final_sorted_data) {
        foreach my $item (keys %$ref) {
            my @rate_array = map { $_ } $ref->{$item}->dog;
            push @plot_data, { $item => \@rate_array };
        }
    } 
    $plot_data{$VOC} = \@plot_data;
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
        q` library(gridExtra) `,
);

$R->set('Time', [@time_axis]);
$R->run(q` plots = list() `);

$R->run(q` my.colours = c(  "CO" = "#2b9eb3" , 
                            "HCHO" = "#e7e85e" , 
                            "TOLUENE" = "#0e5c28" ,
                            "C5H8" = "#8c1531" ,
                            "GLYOX" = "#77aecc" ,
                            "MGLYOX" = "#f7c56c" ,
                            "ACCOMECHO" = "#4c9383" ,
                            "HOCH2CHO" = "#ba8b01" ,
                            "MACR" = "#b569b3",
                            "ACETOL" = "#8ed6d2",
                            "MVK" = "#6db875",
                            "HCOCO2H" = "#f7c56c",
                            "MMALANHY" = "#0c3f78",
                            "C5COO2NO2" = "#cc6329",
                            "MALANHY" = "#1b695b",
                            "NC4MDCO2H" = "#c9a415",
                            "Others" = "#58691b" ) `,
        q` my.names = c( "TOLUENE" = "Toluene" , "C5H8" = "Isoprene" ) `,
);

$R->run(q` plotting = function (data) {   plot = ggplot(data = data, aes(x = Time, y = Reactivity, fill = Reactant)) ;
                                            plot = plot + geom_area(position = "stack") ;
                                            plot = plot + geom_area(position = "stack", colour = "black", show_guide = FALSE) ;
                                            plot = plot + scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7,1)) ;
                                            plot = plot + scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, 0.5)) ;
                                            plot = plot + theme_bw() ;
                                            plot = plot + theme(axis.title.x = element_blank()) ;
                                            plot = plot + theme(axis.title.y = element_blank()) ;
                                            plot = plot + theme(panel.border = element_blank()) ;
                                            plot = plot + theme(axis.line = element_line(colour = "black")) ;
                                            plot = plot + theme(panel.grid.major = element_blank()) ;
                                            plot = plot + theme(panel.grid.minor = element_blank()) ;
                                            plot = plot + theme(axis.text.x = element_text(size = 20)) ;
                                            plot = plot + theme(axis.text.y = element_text(size = 20)) ;
                                            plot = plot + theme(legend.title = element_blank()) ;
                                            plot = plot + theme(legend.text = element_text(size = 20)) ;
                                            plot = plot + theme(legend.key = element_blank()) ;
                                            plot = plot + theme(legend.key.size = unit(1.3, "cm")) ;
                                            plot = plot + theme(legend.justification = c(0.99, 0.7)) ;
                                            plot = plot + theme(legend.position = c(0.99, 0.7)) ;
                                            plot = plot + scale_fill_manual(values = my.colours, labels = my.names) ;
                                            return(plot) } `,
);

foreach my $VOC (sort keys %plot_data) {
    $R->run(q` data = data.frame(Time) `);
    foreach my $ref (@{$plot_data{$VOC}}) {
        foreach my $reactant (sort keys %$ref) {
            $R->set('Reactant', $reactant);
            $R->set('Reactivity', [@{$ref->{$reactant}}]);
            $R->run(q` data[Reactant] = Reactivity `);
        }
    }
    $R->run(q` data = ddply(data, .(Time), colwise(sum)) `, #arrange data
            q` data = melt(data, id.vars = c("Time"), variable.name = "Reactant", value.name = "Reactivity") `,
            q` data$Reactant = factor(data$Reactant, levels = rev(levels(data$Reactant))) `,
            q` plot = plotting(data) `,
            q` plots = c(plots, list(plot)) `,
    );
}

$R->run(q` CairoPDF(file = "Isoprene_Toluene_OH_reactivity.pdf", width = 20, height = 14) `, #save plot to file
        q` multiplot = grid.arrange(    arrangeGrob(plots[[1]] + annotate("text", x = 4, y = 1.4, label = "Isoprene", size = 15, fontface = "bold"),
                                                    plots[[2]] + annotate("text", x = 4, y = 1.4, label = "Toluene", size = 15, fontface = "bold")+ theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                                                    nrow = 1), 
                                       nrow = 1, ncol = 1,
                                       main = textGrob("OH Reactivity from Isoprene and Toluene Degradation", gp = gpar(fontsize = 32, fontface = "bold")),
                                       sub = textGrob("Time (days)\n", gp = gpar(fontsize = 25, fontface = "bold")) ,
                                       left = textGrob(expression(bold(paste("Reactivity (", s ^-1, ")"))), rot = 90, gp = gpar(fontsize = 25), vjust = 0.5) ) `,
        q` print(multiplot) `,
        q` dev.off() `,
);

#my $p = $R->run(q` print(data) `);
#print "$p\n";

$R->stop();
