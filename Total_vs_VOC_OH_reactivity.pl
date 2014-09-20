#! /usr/bin/env perl
# plot total OH reactivity and OH reactivity just of emitted initial VOC OH reactivity
# Version 0: Jane Coates 20/9/2014

use strict; 
use diagnostics;
use MECCA;
use KPP;
use PDL;
use PDL::NiceSlice;
use Statistics::R;
