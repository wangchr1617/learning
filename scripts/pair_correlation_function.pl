#!/usr/bin/perl

use strict;
use warnings;

#configuration for which ensemble average is to be calculated
my $confmin=1;            #starting index of configurations in XDATCAR file for pair correlation function
my $confmax=20000;           #last index of configurations in XDATCAR file for pair correlation function
my $confskip=1;           #stepsize for configuration loop
my $species_1 = 1;        #species 1 for which pair correlation function is going to be calculated
my $species_2 = 1;        #species 2 for which pair correlation function is going to be calculated
#setting radial grid 
my $rmin=0.0;             #minimal value of radial grid
my $rmax=10.0;            #maximum value of radial grid
my $nr=300;                #number of equidistant steps in radial grid
my $dr=($rmax-$rmin)/$nr; #stepsize in radial grid
my $tol=0.0000000001;     #tolerance limit for r->0

my $z=0;                  #counter
my $numelem;              #number of elements
my @elements;             #number of atoms per element saved in list/array
my $lattscale;            #scaling factor for lattice
my @b;                    #Bravais matrix
my $nconf=0;              #number of configurations in XDATCAR file
my @cart;                 #Cartesian coordinates for each atom and configuration
my $atmin_1=0;            #first index of species one
my $atmax_1;              #last index of species one
my $atmin_2=0;            #first index of species two
my $atmax_2;              #last index of species two
my @vol;                  #volume of cell (determinant of Bravais matrix)
my $pi=4*atan2(1, 1);     #constant pi
my $natom=0;              #total number of atoms in cell
my @pcf;                  #pair correlation function (list/array)
my $mult_x=1;             #periodic repetition of cells in x dimension
my $mult_y=1;             #periodic repetition of cells in y dimension
my $mult_z=1;             #periodic repetition of cells in z dimension
my @cart_super;           #Cartesian cells over multiple cells
my @vec_len;              #Length of lattice vectors in 3 spatial coordinates
#my $ensemble_type="NpT";  #Set Npt or NVT. Needs to be set since both have different XDATCAR file.
my $ensemble_type="NVT";  #Set Npt or NVT. Needs to be set since both have different XDATCAR file.
my $av_vol=0;               #Average volume in cell


#reading in XDATCAR file
while (<>)
{
   chomp;
   $_=~s/^/ /;
   my @help=split(/[\t,\s]+/);
   $z++;
   if ($z==2)
   {
      $lattscale = $help[1];
   }
   if ($z==3)
   {
      $b[$nconf+1][1][1]=$help[1]*$lattscale;
      $b[$nconf+1][1][2]=$help[2]*$lattscale;
      $b[$nconf+1][1][3]=$help[3]*$lattscale;
   }
   if ($z==4)
   {
      $b[$nconf+1][2][1]=$help[1]*$lattscale;
      $b[$nconf+1][2][2]=$help[2]*$lattscale;
      $b[$nconf+1][2][3]=$help[3]*$lattscale;
   }
   if ($z==5)
   {
      $b[$nconf+1][3][1]=$help[1]*$lattscale;
      $b[$nconf+1][3][2]=$help[2]*$lattscale;
      $b[$nconf+1][3][3]=$help[3]*$lattscale;
   }
   if ($z==7)
   {
      if ($nconf==0)
      {
         $numelem=@help-1;
         for (my $i=1;$i<=$numelem;$i++)
         {
            $elements[$i]=$help[$i];
            $natom=$natom+$help[$i];
         }
      }
   }
   if ($_=~m/Direct/) 
   {
      $nconf=$nconf+1;
      #for NVT ensemble only one Bravais matrix exists, so it has to be copied
      if ($ensemble_type eq "NVT")
      {
         for (my $i=1;$i<=3;$i++)
         {
            for (my $j=1;$j<=3;$j++)
            {
               $b[$nconf][$i][$j]=$b[1][$i][$j];
            }
         }
      }
      for (my $i=1;$i<=$natom;$i++)
      {
         $_=<>;
         chomp;
         $_=~s/^/ /;
         my @helpat=split(/[\t,\s]+/);
         $cart[$nconf][$i][1]=$b[1][1][1]*$helpat[1]+$b[1][1][2]*$helpat[2]+$b[1][1][3]*$helpat[3];
         $cart[$nconf][$i][2]=$b[1][2][1]*$helpat[1]+$b[1][2][2]*$helpat[2]+$b[1][2][3]*$helpat[3];
         $cart[$nconf][$i][3]=$b[1][3][1]*$helpat[1]+$b[1][3][2]*$helpat[2]+$b[1][3][3]*$helpat[3];
      } 
      if ($ensemble_type eq "NpT")
      {
         $z=0;
      }
   }
   last if eof;
}

if ($confmin>$nconf)
{
   print "Error, confmin larger than number of configurations. Exiting...\n";
   exit;
}
if ($confmax>$nconf)
{
   $confmax=$nconf;
}

for (my $i=1;$i<=$nconf;$i++)
{
   #calculate lattice vector lengths
   $vec_len[$i][1]=($b[$i][1][1]*$b[$i][1][1]+$b[$i][1][2]*$b[$i][1][2]+$b[$i][1][3]*$b[$i][1][3])**0.5;
   $vec_len[$i][2]=($b[$i][2][1]*$b[$i][2][1]+$b[$i][2][2]*$b[$i][2][2]+$b[$i][2][3]*$b[$i][2][3])**0.5;
   $vec_len[$i][3]=($b[$i][3][1]*$b[$i][3][1]+$b[$i][3][2]*$b[$i][3][2]+$b[$i][3][3]*$b[$i][3][3])**0.5;
   #calculate volume of cell
   $vol[$i]=$b[$i][1][1]*$b[$i][2][2]*$b[$i][3][3]+$b[$i][1][2]*$b[$i][2][3]*$b[$i][3][1]+$b[$i][1][3]*$b[$i][2][1]*$b[$i][3][2]-$b[$i][3][1]*$b[$i][2][2]*$b[$i][1][3]-$b[$i][3][2]*$b[$i][2][3]*$b[$i][1][1]-$b[$i][3][3]*$b[$i][2][1]*$b[$i][1][2];
   $av_vol=$av_vol+$vol[$i];
}
$av_vol=$av_vol/$nconf;

#choose species 1 for which pair correlation function is going to be calculated
$atmin_1=1;
if ($species_1>1)
{
   for (my $i=1;$i<$species_1;$i++)
   {
     $atmin_1=$atmin_1+$elements[$i];
   }
}
$atmax_1=$atmin_1+$elements[$species_1]-1;
#choose species 2 to which paircorrelation function is calculated to
$atmin_2=1;
if ($species_2>1)
{
   for (my $i=1;$i<$species_2;$i++)
   {
     $atmin_2=$atmin_2+$elements[$i];
   }
}
$atmax_2=$atmin_2+$elements[$species_2]-1;
#initialize pair correlation function
for (my $i=0;$i<=($nr-1);$i++)
{
   $pcf[$i]=0.0;
}
# loop over configurations, make histogram of pair correlation function 
for (my $j=$confmin;$j<=$confmax;$j=$j+$confskip)
{
   for (my $k=$atmin_1;$k<=$atmax_1;$k++)
   {
       for (my $l=$atmin_2;$l<=$atmax_2;$l++)
       {
          if ($k==$l) {next};
          for (my $g_x=-$mult_x;$g_x<=$mult_x;$g_x++)
          {
             for (my $g_y=-$mult_y;$g_y<=$mult_y;$g_y++)
             {
                for (my $g_z=-$mult_y;$g_z<=$mult_z;$g_z++)
                {
                   my $at2_x=$cart[$j][$l][1]+$vec_len[$j][1]*$g_x;
                   my $at2_y=$cart[$j][$l][2]+$vec_len[$j][2]*$g_y;
                   my $at2_z=$cart[$j][$l][3]+$vec_len[$j][3]*$g_z;
                   my $dist=($cart[$j][$k][1]-$at2_x)**2.0+($cart[$j][$k][2]-$at2_y)**2.0+($cart[$j][$k][3]-$at2_z)**2.0;
                   $dist=$dist**0.5;
                   #determine integer multiple 
                   my $zz=int(($dist-$rmin)/$dr+0.5);
                   if ($zz<$nr)
                   {
                      $pcf[$zz]=$pcf[$zz]+1.0;
                   }
                }
             }
          }
       }
   }
}


#make ensemble average, rescale functions and print
for (my $i=0;$i<=($nr-1);$i++)
{
   my $r=$rmin+$i*$dr;
   if ($r<$tol)
   {
      $pcf[$i]=0.0;
   }
   else
   {
      $pcf[$i]=$pcf[$i]*$av_vol/(4*$pi*$r*$r*$dr*(($confmax-$confmin)/$confskip)*($atmax_2-$atmin_2+1)*($atmax_1-$atmin_1+1));#*((2.0*$mult_x+1.0)*(2.0*$mult_y+1.0)*(2.0*$mult_z+1.0)));
   }
   print $r," ",$pcf[$i],"\n";
}
