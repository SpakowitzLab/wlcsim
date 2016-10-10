#!/usr/bin/perl -w

#use strict;
#require 'calcVecs.pl';

my $npara=41;       # Number of simulations
my $paraind=4;    # Index of parameter to vary
my $val1=0.0;       # First value of parameter
my $val2=1.0;       # Second value of parameter
my $logspace=0;   # Are parameter values log spaced?
my $intval=0;     # Is the parameter an integer?

my $ii=1;
my $val;

system("cp input/input input/input-save");

while ($ii <= $npara) {

# Evaluate the value (either linear or log spaced)
    if ($logspace == 1) {
        $val=$val1*exp(($ii-1)/($npara-1)*log($val2/$val1));
    } else {
        $val=$val1+($ii-1)/($npara-1)*($val2-$val1);
    }

# Reset value to an integer (if $intval==1)

    if ($intval == 1) {
        $val=int($val+$val/abs($val*2));
    }

# Alter the value in the input file

    my $filein="input/input-save";
    my $fileout="input/input";
    open my $in,  '<',  $filein      or die "Can't read old file: $!";
    open my $out, '>',  $fileout     or die "Can't write new file: $!";

    while( <$in> )   # print the lines before the change
    {
        print $out $_;
        last if $. == 2+4*$paraind-1; # line number before change
    }

    my $line = <$in>;
    $line = "  $val\n";
    print $out $line;

    while( <$in> )   # print the rest of the lines
    {
        print $out $_;
    }

    close $in;
    close $out;

# Run simulation and move data to savedata directory in numbered folder

    unless(-e "savedata" or mkdir "savedata") { die "Unable to stat savedata dir"; }
    system("./runsim.sh");
    system("cp -r data savedata/data$ii");
    system("cp input/input savedata/data$ii/input$ii");

    $ii++;
}
