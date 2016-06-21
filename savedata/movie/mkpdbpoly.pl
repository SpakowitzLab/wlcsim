#!/usr/bin/perl -w

#use strict;
#require 'calcVecs.pl';

# Converts coordinate into properly formatted pdb files
# Written by Andy Spakowitz 5-20-05

my $total=2000;         # Total number of snapshots
my $skip=1;            # Skip number of snapshots
my $np=1;               # Number of polymers
my $n0=2000;

my $radius=1;           # Radius of chain
my $ratio=1;            # Resize ratio
my $snapcount=1;        # File count for output
my $filecount=1;        # File count for input

my $nbpbead=1;         # Number of bp per bead
my $nbpturn=10;         # Number of bp per turn
my $phidna=135;         # Rotation angle for the double helix
my $pi=3.14159265359;   # Value of pi
my $gamma=2*$pi*$nbpbead/$nbpturn; # Twist angle per bead

my $filein1;            # File with bead coordinates
my $fileout1;           # Output file for pdb

while($filecount<=$total){

   if ($snapcount < 10)
   {
       $filein1="xyz/r$filecount";
       $fileout1=">pdb/snap00$snapcount.pdb";
   }
   elsif (($snapcount >= 10) && ($snapcount < 100))
   {
       $filein1="xyz/r$filecount";
       $fileout1=">pdb/snap0$snapcount.pdb";
   }
   else
   {
       $filein1="xyz/r$filecount";
       $fileout1=">pdb/snap$snapcount.pdb";
   }
   open(COORD1, $filein1) || die('cannot open file:'. $!);
   open(PDB1, $fileout1) || die('cannot open file:'. $!);

   my $count=1;

   # Array to hold positions of backbone
   my @atomx1;
   my @atomy1;
   my @atomz1;
   my @atomx4;
   my @atomy4;
   my @atomz4;

   while(<COORD1>){
      my @info = split;
      $atomx1[$count] = $ratio*$info[0];
      $atomy1[$count] = $ratio*$info[1];
      $atomz1[$count] = $ratio*$info[2];
   
      $count++;

   }
   my $nbead = $count-1;          # Index of last element in atomx1
   my $n=$nbead/$np;

   ##############################################################
   # Assemble single PDB file
   ##############################################################

   my $atomname1 = "A1";           # Chain atom type
   my $atomname2 = "A2";           # Ribbon atom type
   my $atomname3 = "A3";           # Extra atom type
   my $atomname4 = "A4";           # Extra atom type
   my $atomname5 = "A5";           # Extra atom type
   my $atomname6 = "A6";           # Extra atom type
   my $atomname7 = "A7";           # Extra atom type
   my $atomname8 = "A8";           # Extra atom type
   my $resname = "SSN";           # Type of residue (UNKnown/Single Stranded Nucleotide)
   my $chain = "A";               # Chain identifier
   my $resnum = "1";
   my $numresidues = $nbead;
   my $descrip = "Pseudo atom representation of DNA";
   my $chemicalname = "Body and ribbon spatial coordinates";
   
   # Het Header info
   printf PDB1 "HET    %3s  %1s%4d   %5d     %-38s\n",$resname,$chain,$resnum, $numresidues,$descrip ;
   printf PDB1 "HETNAM     %3s %-50s\n",$resname, $chemicalname;
   printf PDB1 "FORMUL  1   %3s    C20 N20 P21\n",$resname;

   $count=1;
   $ii=1;
   while ($ii <= $np) {
       $jj=1;
       $rat=($n-5)/($n0-5);
       while ($jj <= $n0) {
	   if ($jj <= $n) {
	       $ind=$jj+$n*($ii-1);
	       printf PDB1 "ATOM%7d %4s %3s %1s        %8.3f%8.3f%8.3f%6.2f%6.2f           C\n",$count,$atomname1,$resname,1,$atomx1[$ind],$atomy1[$ind],$atomz1[$ind],$rat,$rat;
	   } else {
	       $ind=$n*$ii;
	       printf PDB1 "ATOM%7d %4s %3s %1s        %8.3f%8.3f%8.3f%6.2f%6.2f           C\n",$count,$atomname1,$resname,1,$atomx1[$ind],$atomy1[$ind],$atomz1[$ind],$rat,$rat;   
	   }
	   $count++;
	   $jj++;
       }
       $ii++;
   }
   
# Connect up the ribbon

   $jj=1;
   $con=1;
   while ($jj <= $np) {
       $ii=1;
       while ($ii <= $n0) {
	   if ($ii == 1) {
	       printf PDB1 "CONECT%5d%5d\n",$con,$con + 1;
	   } elsif (($ii > 1) && ($ii < $n0)) {
	       printf PDB1 "CONECT%5d%5d%5d\n",$con,$con - 1, $con + 1;
	   } else {
	       printf PDB1 "CONECT%5d%5d\n",$con,$con - 1;
	   }
	   $con++;
	   $ii++;
       }
       $jj++;
   }

   printf PDB1 "END";
   # Clean up and close files

   close(PDB1);
   close(COORD1);
   $snapcount++;
   $filecount=$filecount+$skip;
}

###########################################################################
sub arccos
{
    atan2( sqrt(1.0 - $_[0] * $_[0]), $_[0] );
}

###########################################################################

