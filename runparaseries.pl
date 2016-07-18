#!/usr/bin/perl -w

#use strict;
#require 'calcVecs.pl';

my $npara=20;       # Number of simulations
my $paraind=6;    # Index of parameter to vary
my $val1=5;       # First value of parameter
my $val2=500;       # Second value of parameter
my $logspace=1;   # Are parameter values log spaced?
my $intval=1;     # Is the parameter an integer?
my $nsim=50;     # Total number of duplicates

my $ii;
my $val;
my $itot;

system("cp input/input input/input-save");

$itot=1;
while ($itot <= $nsim) {
    $ii=1;
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

    unless(-e savedata or mkdir savedata) { die "Unable to stat savedata dir"; }
    system("./runsim");
    if ($itot == 1) {
        system("mkdir savedata/data$ii");
    }
    system("cp data/out1 savedata/data$ii/out1-$itot");
    system("cp data/out2 savedata/data$ii/out2-$itot");
    system("cp data/out3 savedata/data$ii/out3-$itot");
    system("cp input/input savedata/data$ii/input$ii");

    $ii++;
    }
    $itot++;
}
