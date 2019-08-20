#!/bin/perl

sub change 
{
  $a=$_[0];

$a=~s/alfa/lapse/g;
$a=~s/alfa/lapse/g;
$a=~s/\(2\*M\+r\)\^\(3\/2\)/pow(2*M+r,1.5)/g;
$a=~s/\(r\+2\*M\)\^\(3\/2\)/pow(2*M+r,1.5)/g;
$a=~s/\(2\*M\+r\)\^\(3\/2\)/pow(2*M+r,1.5)/g;
$a=~s/\(r\+2\*M\)\^\(3\/2\)/pow(2*M+r,1.5)/g;
$a=~s/Vt1\^2/pow(Vt1,2.0)/g;
$a=~s/Vt2\^2/pow(Vt2,2.0)/g;
$a=~s/Vt3\^2/pow(Vt3,2.0)/g;
$a=~s/([A-Za-z])\^(\d\d)/pow($1,$2$3.0)/g;
$a=~s/([A-Za-z])\^(\d)/pow($1,$2.0)/g;
$a=~s/([A-Za-z])\^\((\d)\/(\d)\)/pow($1,$2.0\/$3.0)/g;
$a=~s/([A-Za-z])\^\((\d)(\d)\/(\d)\)/pow($1,$2$3.0\/$4.0)/g;
$a=~s/abs/fabs/g;
$a=~s/'diff\(([A-Za-z]),([a-z]),1\)/d$1$2/g;
$a=~s/r/x1/g;
$a=~s/theta/x2/g;
$a=~s/phi/x3/g;
$a=~s/sqx1t/sqrt/g;
$a=~s/x2t11/yt11/g;
$a=~s/x2t22/yt22/g;
$a=~s/x2t33/yt33/g;
$a=~s/M/MM/g;
$a=~s/k/K/g;
$a=~s/(sin\(x[1-3]\))\^(\d)/pow($1,$2.0)/g;
$a=~s/(sin\(x[1-3]\))\^\((\d)\/(\d)\)/pow($1,$2.0\/$3.0)/g;
$a=~s/(sin\(x[1-3]\))\^\((\d)(\d)\/(\d)\)/pow($1,$2$3.0\/$4.0)/g;
$a=~s/(cos\(x[1-3]\))\^(\d)/pow($1,$2.0)/g;
$a=~s/(cos\(x[1-3]\))\^\((\d)\/(\d)\)/pow($1,$2.0\/$3.0)/g;
$a=~s/(cos\(x[1-3]\))\^\((\d)(\d)\/(\d)\)/pow($1,$2$3.0\/$4.0)/g;
$a=~s/(tan\(x[1-3]\))\^(\d)/pow($1,$2.0)/g;
$a=~s/(tan\(x[1-3]\))\^\((\d)\/(\d)\)/pow($1,$2.0\/$3.0)/g;
$a=~s/(tan\(x[1-3]\))\^\((\d)(\d)\/(\d)\)/pow($1,$2$3.0\/$4.0)/g;
return $a;
}

########################################################
##Headers##
########################################################

open(Matrix, ">../Headers/vector.h");
print Matrix "int funct_A(double *a, double *uu);\n";
print Matrix "int funct_Dm(double *a, double *uu);\n";
print Matrix "int funct_Dn(double *a, double *uu);\n";
print Matrix "int funct_Do(double *a, double *uu);\n";
print Matrix "int funct_U2Q(double *a, double *uu);\n";
print Matrix "int funct_Q2U(double *a, double *uu);\n";
print Matrix "int funct_Q(double *a, double *uu);\n";
print Matrix "int funct_F(double *a, double *uu);\n";
print Matrix "int funct_G(double *a, double *uu);\n";
print Matrix "int funct_H(double *a, double *uu);\n";
print Matrix "int funct_S(double *a, double *uu);\n";
print Matrix "int GAUGE(double *a, double g1, double g2, double g3);\n";
close(Matrix);

########################################################
##MATRIX A##
########################################################

open(Matrix, 'DIFF');
@lines = <Matrix>;
my $size = @lines;
my $s = $size - 1;
my @b = (0) x @lines;
#print @lines;

foreach my $i (1..$size){
$b[$i-1] = $lines[$i];
}

foreach my $i (1..$size){
$b[$i-1] = change $b[$i-1];
}

open(Matrix, 'PARAM');
@lines = <Matrix>;
my $size = @lines;
my $s = $size - 1;
my @c = (0) x @lines;
#print @lines;

foreach my $i (1..$size){
$c[$i-1] = $lines[$i];
}

foreach my $i (1..$size){
$c[$i-1] = change $c[$i-1];
}

open(Matrix, 'MATRIX_A');
@lines = <Matrix>;
my $size = @lines;
my $s = $size - 1;
my @a = (0) x @lines;
#print @lines;

foreach my $i (1..$size){
$a[$i-1] = $lines[$i];
}

foreach my $i (1..$size){
$a[$i-1] = change $a[$i-1];
}

open(Matrix, ">amatrix.c");
print Matrix "#include<stdio.h>\n";
print Matrix "#include<math.h>\n";
print Matrix "#include\"../Headers/vector.h\"\n";
print Matrix "#include\"../Headers/main.h\"\n";
print Matrix "\n";

print Matrix "int funct_A(double *a, double *uu)\n";
print Matrix "{\n";                                                 
print Matrix "   int i, j;\n";
print Matrix "   double n, p, u=0, v=0, w=0;\n";
print Matrix "   double R, W, h;\n";
print Matrix "   double dWu, dWv, dWw;\n";
print Matrix "   double dhn, dhp;\n";
print Matrix "   n = uu[0];\n";
print Matrix "   p = uu[1];\n";
print Matrix "   u = uu[2];\n";
print Matrix "   if(dim >= 2){v = uu[3];}\n";
print Matrix "   if(dim == 3){w = uu[4];}\n";
print Matrix "\n";
print Matrix "   R = $c[0]";
print Matrix "   W = $c[1]";
print Matrix "   h = $c[4]";
print Matrix "\n";
print Matrix "   dWu = $b[0]";
print Matrix "   dWv = $b[1]";
print Matrix "   dWw = $b[2]";
print Matrix "   dhn = $b[3]";
print Matrix "   dhp = $b[4]";
print Matrix "\n";
foreach my $i (0..4){
   foreach my $j (0..4){   
      print Matrix "   a[$i*eq + $j] = $a[$i*5 + $j]";     
   }   
}
print Matrix "\n";
print Matrix "   return 0;\n";                                       
print Matrix "}\n";                                                 
close(Matrix);

########################################################
##MATRIX Dm##
########################################################

open(Matrix, 'MATRIX_Dm');
@lines = <Matrix>;
my $size = @lines;
my $s = $size - 1;
my @a = (0) x @lines;
#print @lines;

foreach my $i (1..$size){
$a[$i-1] = $lines[$i];
}

foreach my $i (1..$size){
$a[$i-1] = change $a[$i-1];
}

open(Matrix, ">dmmatrix.c");
print Matrix "#include<stdio.h>\n";
print Matrix "#include<math.h>\n";
print Matrix "#include\"../Headers/vector.h\"\n";
print Matrix "#include\"../Headers/main.h\"\n";
print Matrix "\n";

print Matrix "int funct_Dm(double *a, double *uu)\n";
print Matrix "{\n";                                                 
print Matrix "   int i;\n";  
print Matrix "   double n, p, u=0, v=0, w=0;\n";
print Matrix "   double R, V;\n";
print Matrix "   double yt11, Bt1, Vt1, c;\n";
print Matrix "   n = uu[0];\n";
print Matrix "   p = uu[1];\n";
print Matrix "   u = uu[2];\n";
print Matrix "   if(dim >= 2){v = uu[3];}\n";
print Matrix "   if(dim == 3){w = uu[4];}\n";
print Matrix "\n";
print Matrix "   R    = $c[0]";
print Matrix "   V    = $c[2]";
print Matrix "   yt11 = $c[6]";
print Matrix "   Bt1  = $c[9]";
print Matrix "   Vt1  = $c[12]";
print Matrix "   c    = $c[5]";
print Matrix "\n";
foreach my $i (0..$#a-1){
print Matrix "   a[$i] = $a[$i]";
}
print Matrix "\n";
print Matrix "   return 0;\n";                                       
print Matrix "}\n";                                                 
close(Matrix);

########################################################
##MATRIX Dn##
########################################################

open(Matrix, 'MATRIX_Dn');
@lines = <Matrix>;
my $size = @lines;
my $s = $size - 1;
my @a = (0) x @lines;
#print @lines;

foreach my $i (1..$size){
$a[$i-1] = $lines[$i];
}

foreach my $i (1..$size){
$a[$i-1] = change $a[$i-1];
}

open(Matrix, ">dnmatrix.c");
print Matrix "#include<stdio.h>\n";
print Matrix "#include<math.h>\n";
print Matrix "#include\"../Headers/vector.h\"\n";
print Matrix "#include\"../Headers/main.h\"\n";
print Matrix "\n";

print Matrix "int funct_Dn(double *a, double *uu)\n";
print Matrix "{\n";                                                 
print Matrix "   int i;\n";  
print Matrix "   double n, p, u=0, v=0, w=0;\n";
print Matrix "   double R, V;\n";
print Matrix "   double yt22, Bt2, Vt2, c;\n";
print Matrix "   n = uu[0];\n";
print Matrix "   p = uu[1];\n";
print Matrix "   u = uu[2];\n";
print Matrix "   if(dim >= 2){v = uu[3];}\n";
print Matrix "   if(dim == 3){w = uu[4];}\n";
print Matrix "\n";
print Matrix "   R     = $c[0]";
print Matrix "   V     = $c[2]";
print Matrix "   yt22  = $c[7]";
print Matrix "   Bt2   = $c[10]";
print Matrix "   Vt2   = $c[13]";
print Matrix "   c     = $c[5]";
print Matrix "\n";
foreach my $i (0..$#a-1){
print Matrix "   a[$i] = $a[$i]";
}  
print Matrix "\n";
print Matrix "   return 0;\n";                                       
print Matrix "}\n";                                                 
close(Matrix);

########################################################
##MATRIX Do##
########################################################

open(Matrix, 'MATRIX_Do');
@lines = <Matrix>;
my $size = @lines;
my $s = $size - 1;
my @a = (0) x @lines;
#print @lines;

foreach my $i (1..$size){
$a[$i-1] = $lines[$i];
}

foreach my $i (1..$size){
$a[$i-1] = change $a[$i-1];
}

open(Matrix, ">domatrix.c");
print Matrix "#include<stdio.h>\n";
print Matrix "#include<math.h>\n";
print Matrix "#include\"../Headers/vector.h\"\n";
print Matrix "#include\"../Headers/main.h\"\n";
print Matrix "\n";

print Matrix "int funct_Do(double *a, double *uu)\n";
print Matrix "{\n";                                                 
print Matrix "   int i;\n";  
print Matrix "   double n, p, u=0, v=0, w=0;\n";
print Matrix "   double R, V;\n";
print Matrix "   double yt33, Bt3, Vt3, c;\n";
print Matrix "   double dWu, dWv, dWw;\n";
print Matrix "   double dhn, dhp;\n";
print Matrix "   n = uu[0];\n";
print Matrix "   p = uu[1];\n";
print Matrix "   u = uu[2];\n";
print Matrix "   if(dim >= 2){v = uu[3];}\n";
print Matrix "   if(dim == 3){w = uu[4];}\n";
print Matrix "\n";
print Matrix "   R    = $c[0]";
print Matrix "   V    = $c[2]";
print Matrix "   yt33 = $c[8]";
print Matrix "   Bt3  = $c[11]";
print Matrix "   Vt3  = $c[14]";
print Matrix "   c    = $c[5]";
print Matrix "\n";
foreach my $i (0..$#a-1){
print Matrix "   a[$i] = $a[$i]";
}  
print Matrix "\n";
print Matrix "   return 0;\n";                                       
print Matrix "}\n";                                                 
close(Matrix);

########################################################
##VECTOR U2Q##
########################################################

open(Matrix, 'VECTOR_Q');
@lines = <Matrix>;
my $size = @lines;
my $s = $size - 1;
my @a = (0) x @lines;
#print @lines;

foreach my $i (1..$size){
$a[$i-1] = $lines[$i];
}

foreach my $i (1..$size){
$a[$i-1] = change $a[$i-1];
}

open(Matrix, ">u2qvector.c");
print Matrix "#include<stdio.h>\n";
print Matrix "#include<math.h>\n";
print Matrix "#include\"../Headers/vector.h\"\n";
print Matrix "#include\"../Headers/main.h\"\n";
print Matrix "\n";

print Matrix "int funct_U2Q(double *a, double *uu)\n";
print Matrix "{\n";                                                 
print Matrix "   int i, j, k;\n";  
print Matrix "   double n, p, u=0, v=0, w=0;\n";
print Matrix "   double R, W, h;\n";
print Matrix "   double dWu, dWv, dWw;\n";
print Matrix "   double dhn, dhp;\n";
print Matrix "   double lapse, dety;\n";
print Matrix "\n";
print Matrix "   if(dim == 1)\n";
print Matrix "   {\n";
print Matrix "      for(i = 0; i <= Nx1; i++)\n";
print Matrix "      {\n";
print Matrix "         n = uu[c1(0,i)];\n";
print Matrix "         p = uu[c1(1,i)];\n";
print Matrix "         u = uu[c1(2,i)];\n";
print Matrix "         v = 0;\n";
print Matrix "         w = 0;\n";
print Matrix "\n";
print Matrix "         x1 = X1[i];\n";
print Matrix "         x2 = 0.0;\n";
print Matrix "         x3 = 0.0;\n";
print Matrix "\n";
print Matrix "         R = $c[0]";
print Matrix "         W = $c[1]";
print Matrix "         h = $c[4]";
print Matrix "\n";
print Matrix "         lapse = $c[18]";
print Matrix "         dety  = $c[19]";
print Matrix "\n";
foreach my $i (0..$#a-1){                                   
print Matrix "         a[c1($i,i)] = $a[$i]";
}
print Matrix "      }\n";
print Matrix "   }\n";
print Matrix "   else if(dim == 2)\n";
print Matrix "   {\n";
print Matrix "      for(i = 0; i <= Nx1; i++)\n";
print Matrix "      {\n";
print Matrix "         for(j = 0; j <= Nx2; j++)\n";
print Matrix "         {\n";
print Matrix "            n = uu[c2(0,i,j)];\n";
print Matrix "            p = uu[c2(1,i,j)];\n";
print Matrix "            u = uu[c2(2,i,j)];\n";
print Matrix "            v = uu[c2(3,i,j)];\n";
print Matrix "            w = 0;\n";
print Matrix "\n";
print Matrix "            x1 = X1[i];\n";
print Matrix "            x2 = X2[j];\n";
print Matrix "            x3 = 0.0;\n";
print Matrix "\n";
print Matrix "            R = $c[0]";
print Matrix "            W = $c[1]";
print Matrix "            h = $c[4]";
print Matrix "\n";
print Matrix "            lapse = $c[18]";
print Matrix "            dety  = $c[19]";
print Matrix "\n";
foreach my $i (0..$#a-1){                                   
print Matrix "            a[c2($i,i,j)] = $a[$i]";
}                                    
print Matrix "         }\n";
print Matrix "      }\n";
print Matrix "   }\n";
print Matrix "   if(dim == 3)\n";
print Matrix "   {\n";
print Matrix "      for(i = 0; i <= Nx1; i++)\n";
print Matrix "      {\n";
print Matrix "         for(j = 0; j <= Nx2; j++)\n";
print Matrix "         {\n";
print Matrix "            for(k = 0; k <= Nx3; k++)\n";
print Matrix "            {\n";
print Matrix "               n = uu[c3(0,i,j,k)];\n";
print Matrix "               p = uu[c3(1,i,j,k)];\n";
print Matrix "               u = uu[c3(2,i,j,k)];\n";
print Matrix "               v = uu[c3(3,i,j,k)];\n";
print Matrix "               w = uu[c3(4,i,j,k)];\n";
print Matrix "\n";
print Matrix "               x1 = X1[i];\n";
print Matrix "               x2 = X2[j];\n";
print Matrix "               x3 = X3[k];\n";
print Matrix "\n";
print Matrix "               R = $c[0]";
print Matrix "               W = $c[1]";
print Matrix "               h = $c[4]";
print Matrix "\n";
print Matrix "               lapse = $c[18]";
print Matrix "               dety  = $c[18]";
print Matrix "\n";
foreach my $i (0..$#a-1){                                   
print Matrix "               a[c3($i,i,j,k)] = $a[$i]";
}                                    
print Matrix "            }\n";
print Matrix "         }\n";
print Matrix "      }\n";
print Matrix "   }\n";
print Matrix "\n";
print Matrix "   return 0;\n";                                       
print Matrix "}\n";                                                 
close(Matrix);

########################################################
##VECTOR Q##
########################################################

open(Matrix, 'VECTOR_Q');
@lines = <Matrix>;
my $size = @lines;
my $s = $size - 1;
my @a = (0) x @lines;
#print @lines;

foreach my $i (1..$size){
$a[$i-1] = $lines[$i];
}

foreach my $i (1..$size){
$a[$i-1] = change $a[$i-1];
}

open(Matrix, ">qvector.c");
print Matrix "#include<stdio.h>\n";
print Matrix "#include<math.h>\n";
print Matrix "#include\"../Headers/vector.h\"\n";
print Matrix "#include\"../Headers/main.h\"\n";
print Matrix "\n";

print Matrix "int funct_Q(double *a, double *uu)\n";
print Matrix "{\n";                                                 
print Matrix "   int i;\n";  
print Matrix "   double n, p, u=0, v=0, w=0;\n";
print Matrix "   double R, W, h;\n";
print Matrix "   double dWu, dWv, dWw;\n";
print Matrix "   double dhn, dhp;\n";
print Matrix "   n = uu[0];\n";
print Matrix "   p = uu[1];\n";
print Matrix "   u = uu[2];\n";
print Matrix "   if(dim >= 2){v = uu[3];}\n";
print Matrix "   if(dim == 3){w = uu[4];}\n";
print Matrix "\n";
print Matrix "   R = $c[0]";
print Matrix "   W = $c[1]";
print Matrix "   h = $c[4]";
print Matrix "\n";
foreach my $i (0..$#a-1){                                   
print Matrix "   a[$i] = $a[$i]";
}                                    
print Matrix "\n";
print Matrix "   return 0;\n";                                       
print Matrix "}\n";                                                 
close(Matrix);

########################################################
##VECTOR F##
########################################################

open(Matrix, 'VECTOR_F');
@lines = <Matrix>;
my $size = @lines;
my $s = $size - 1;
my @a = (0) x @lines;
#print @lines;

foreach my $i (1..$size){
$a[$i-1] = $lines[$i];
}

foreach my $i (1..$size){
$a[$i-1] = change $a[$i-1];
}

open(Matrix, ">fvector.c");
print Matrix "#include<stdio.h>\n";
print Matrix "#include<math.h>\n";
print Matrix "#include\"../Headers/vector.h\"\n";
print Matrix "#include\"../Headers/main.h\"\n";
print Matrix "\n";

print Matrix "int funct_F(double *a, double *uu)\n";
print Matrix "{\n";                                                 
print Matrix "   int i;\n";  
print Matrix "   double n, p, u=0, v=0, w=0;\n";
print Matrix "   double R, W, h;\n";
print Matrix "   double dWu, dWv, dWw;\n";
print Matrix "   double dhn, dhp;\n";
print Matrix "   n = uu[0];\n";
print Matrix "   p = uu[1];\n";
print Matrix "   u = uu[2];\n";
print Matrix "   if(dim >= 2){v = uu[3];}\n";
print Matrix "   if(dim == 3){w = uu[4];}\n";
print Matrix "\n";
print Matrix "   R = $c[0]";
print Matrix "   W = $c[1]";
print Matrix "   h = $c[4]";
print Matrix "\n";
foreach my $i (0..$#a-1){                                   
print Matrix "   a[$i] = $a[$i]";
}                                    
print Matrix "\n";
print Matrix "   return 0;\n";                                       
print Matrix "}\n";                                                 
close(Matrix);

########################################################
##VECTOR G##
########################################################

open(Matrix, 'VECTOR_G');
@lines = <Matrix>;
my $size = @lines;
my $s = $size - 1;
my @a = (0) x @lines;
#print @lines;

foreach my $i (1..$size){
$a[$i-1] = $lines[$i];
}

foreach my $i (1..$size){
$a[$i-1] = change $a[$i-1];
}

open(Matrix, ">gvector.c");
print Matrix "#include<stdio.h>\n";
print Matrix "#include<math.h>\n";
print Matrix "#include\"../Headers/vector.h\"\n";
print Matrix "#include\"../Headers/main.h\"\n";
print Matrix "\n";

print Matrix "int funct_G(double *a, double *uu)\n";
print Matrix "{\n";                                                 
print Matrix "   int i;\n";  
print Matrix "   double n, p, u=0, v=0, w=0;\n";
print Matrix "   double R, W, h;\n";
print Matrix "   double dWu, dWv, dWw;\n";
print Matrix "   double dhn, dhp;\n";
print Matrix "   n = uu[0];\n";
print Matrix "   p = uu[1];\n";
print Matrix "   u = uu[2];\n";
print Matrix "   if(dim >= 2){v = uu[3];}\n";
print Matrix "   if(dim == 3){w = uu[4];}\n";
print Matrix "\n";
print Matrix "   R = $c[0]";
print Matrix "   W = $c[1]";
print Matrix "   h = $c[4]";
print Matrix "\n";
foreach my $i (0..$#a-1){                                   
print Matrix "   a[$i] = $a[$i]";
}
print Matrix "\n";
print Matrix "   return 0;\n";                                       
print Matrix "}\n";                                                 
close(Matrix);

########################################################
##VECTOR H##
########################################################

open(Matrix, 'VECTOR_H');
@lines = <Matrix>;
my $size = @lines;
my $s = $size - 1;
my @a = (0) x @lines;
#print @lines;

foreach my $i (1..$size){
$a[$i-1] = $lines[$i];
}

foreach my $i (1..$size){
$a[$i-1] = change $a[$i-1];
}

open(Matrix, ">hvector.c");
print Matrix "#include<stdio.h>\n";
print Matrix "#include<math.h>\n";
print Matrix "#include\"../Headers/vector.h\"\n";
print Matrix "#include\"../Headers/main.h\"\n";
print Matrix "\n";

print Matrix "int funct_H(double *a, double *uu)\n";
print Matrix "{\n";                                                 
print Matrix "   int i;\n";  
print Matrix "   double n, p, u=0, v=0, w=0;\n";
print Matrix "   double R, W, h;\n";
print Matrix "   double dWu, dWv, dWw;\n";
print Matrix "   double dhn, dhp;\n";
print Matrix "   n = uu[0];\n";
print Matrix "   p = uu[1];\n";
print Matrix "   u = uu[2];\n";
print Matrix "   if(dim >= 2){v = uu[3];}\n";
print Matrix "   if(dim == 3){w = uu[4];}\n";
print Matrix "\n";
print Matrix "   R = $c[0]";
print Matrix "   W = $c[1]";
print Matrix "   h = $c[4]";
print Matrix "\n";
foreach my $i (0..$#a-1){                                   
print Matrix "   a[$i] = $a[$i]";
}                                    
print Matrix "\n";
print Matrix "   return 0;\n";                                       
print Matrix "}\n";                                                 
close(Matrix);

########################################################
##VECTOR S##
########################################################

open(Matrix, 'VECTOR_S');
@lines = <Matrix>;
my $size = @lines;
my $s = $size - 1;
my @a = (0) x @lines;
#print @lines;

foreach my $i (1..$size){
$a[$i-1] = $lines[$i];
}

foreach my $i (1..$size){
$a[$i-1] = change $a[$i-1];
}

open(Matrix, ">svector.c");
print Matrix "#include<stdio.h>\n";
print Matrix "#include<math.h>\n";
print Matrix "#include\"../Headers/vector.h\"\n";
print Matrix "#include\"../Headers/main.h\"\n";
print Matrix "\n";

print Matrix "int funct_S(double *a, double *uu)\n";
print Matrix "{\n";                                                 
print Matrix "   int i;\n";  
print Matrix "   double n, p, u=0, v=0, w=0;\n";
print Matrix "   double R, W, h;\n";
print Matrix "   double dWu, dWv, dWw;\n";
print Matrix "   double dhn, dhp;\n";
print Matrix "   n = uu[0];\n";
print Matrix "   p = uu[1];\n";
print Matrix "   u = uu[2];\n";
print Matrix "   if(dim >= 2){v = uu[3];}\n";
print Matrix "   if(dim == 3){w = uu[4];}\n";
print Matrix "\n";
print Matrix "   R = $c[0]";
print Matrix "   W = $c[1]";
print Matrix "   h = $c[4]";
print Matrix "\n";
foreach my $i (0..$#a-1){                                   
print Matrix "   a[$i] = $a[$i]";
}
print Matrix "\n";
print Matrix "   return 0;\n";                                       
print Matrix "}\n";                                                 
close(Matrix);

########################################################
##VECTOR Q2U##
########################################################

open(Matrix, ">q2uvector.c");
print Matrix "#include<stdio.h>\n";
print Matrix "#include<math.h>\n";
print Matrix "#include\"../Headers/vector.h\"\n";
print Matrix "#include\"../Headers/main.h\"\n";
print Matrix "\n";
print Matrix "int funct_Q2U(double *a, double *uu)\n";
print Matrix "{\n";
print Matrix "   int i, j, k;\n";
print Matrix "   double D, t, m=0, n=0, o=0;\n";
print Matrix "   double R, SS;\n";
print Matrix "   double dWu, dWv, dWw;\n";
print Matrix "   double dhn, dhp;\n";
print Matrix "   double theta, theta_0;\n";
print Matrix "   double f, derf, h, derh, lor;\n";
print Matrix "\n";
print Matrix "   R  = sqrt(pow(x2,2.0)+pow(x1,2.0));\n";
print Matrix "\n";
print Matrix "   if(dim == 1)\n";
print Matrix "   {\n";
print Matrix "      for(i = 0; i <= Nx1; i++)\n";
print Matrix "      {\n";
print Matrix "         D = uu[c1(0,i)];\n";
print Matrix "         t = uu[c1(1,i)];\n";
print Matrix "         m = uu[c1(2,i)];\n";
print Matrix "         n = 0;\n";
print Matrix "         o = 0;\n";
print Matrix "\n";
print Matrix "         x1 = X1[i];\n";
print Matrix "         x2 = 0;\n";
print Matrix "         x3 = 0;\n";
print Matrix "\n";
print Matrix "         SS = $c[3]";
print Matrix "\n";
print Matrix "         theta_0 = 0.0;\n";
print Matrix "         f = 2.0;\n";
print Matrix "\n";
print Matrix "         h   = 1.0;\n";
print Matrix "         lor  = sqrt(1.0 + SS / pow(D*h,2.0));\n";
print Matrix "\n";
print Matrix "         a[c1(0,i)] = D / lor;\n";
print Matrix "         a[c1(1,i)] = 0.0;\n";
print Matrix "         a[c1(2,i)] = m / (D * h * lor);\n";
print Matrix "      }\n";
print Matrix "   }\n";
print Matrix "   else if(dim == 2)\n";
print Matrix "   {\n";
print Matrix "      for(i = 0; i <= Nx1; i++)\n";
print Matrix "      {\n";
print Matrix "         for(j = 0; j <= Nx2; j++)\n";
print Matrix "         {\n";
print Matrix "            D = uu[c2(0,i,j)];\n";
print Matrix "            t = uu[c2(1,i,j)];\n";
print Matrix "            m = uu[c2(2,i,j)];\n";
print Matrix "            n = uu[c2(3,i,j)];\n";
print Matrix "            o = 0;\n";
print Matrix "\n";
print Matrix "            x1 = X1[i];\n";
print Matrix "            x2 = X2[j];\n";
print Matrix "            x3 = 0;\n";
print Matrix "\n";
print Matrix "            SS = $c[3]";
print Matrix "\n";
print Matrix "            theta_0 = 0.0;\n";
print Matrix "            f = 2.0;\n";
print Matrix "\n"; 
print Matrix "            h   = 1.0;\n";
print Matrix "            lor  = sqrt(1.0 + SS / pow(D*h,2.0));\n";
print Matrix "\n";
print Matrix "            a[c2(0,i,j)] = D / lor;\n";
print Matrix "            a[c2(1,i,j)] = 0.0;\n";
print Matrix "            a[c2(2,i,j)] = m / (D * h * lor);\n";
print Matrix "            a[c2(3,i,j)] = n / (D * h * lor);\n";
print Matrix "         }\n";
print Matrix "      }\n";
print Matrix "   }\n";
print Matrix "   if(dim == 3)\n";
print Matrix "   {\n";
print Matrix "      for(i = 0; i <= Nx1; i++)\n";
print Matrix "      {\n";
print Matrix "         for(j = 0; j <= Nx2; j++)\n";
print Matrix "         {\n";
print Matrix "            for(k = 0; k <= Nx3; k++)\n";
print Matrix "            {\n";
print Matrix "               D = uu[c3(0,i,j,k)];\n";
print Matrix "               t = uu[c3(1,i,j,k)];\n";
print Matrix "               m = uu[c3(2,i,j,k)];\n";
print Matrix "               n = uu[c3(3,i,j,k)];\n";
print Matrix "               o = uu[c3(4,i,j,k)];\n";
print Matrix "\n";
print Matrix "               x1 = X1[i];\n";
print Matrix "               x2 = X2[j];\n";
print Matrix "               x3 = X3[k];\n";
print Matrix "\n";
print Matrix "               SS = $c[3]";
print Matrix "\n";
print Matrix "               theta_0 = 0.0;\n";
print Matrix "               f = 2.0;\n";
print Matrix "\n";
print Matrix "               h   = 1.0;\n";
print Matrix "               lor  = sqrt(1.0 + SS / pow(D*h,2.0));\n";
print Matrix "\n";
print Matrix "               a[c3(0,i,j,k)] = D / lor;\n";
print Matrix "               a[c3(1,i,j,k)] = 0.0;\n";
print Matrix "               a[c3(2,i,j,k)] = m / (D * h * lor);\n";
print Matrix "               a[c3(3,i,j,k)] = n / (D * h * lor);\n";
print Matrix "               a[c3(4,i,j,k)] = o / (D * h * lor);\n";
print Matrix "            }\n";
print Matrix "         }\n";
print Matrix "      }\n";
print Matrix "   }\n";
print Matrix "\n";
print Matrix "   return 0;\n";
print Matrix "}\n";
close(Matrix);

open(Matrix, ">gauge.c");
print Matrix "#include<stdio.h>\n";
print Matrix "#include<math.h>\n";
print Matrix "#include\"../Headers/vector.h\"\n";
print Matrix "#include\"../Headers/main.h\"\n";
print Matrix "\n";
print Matrix "int GAUGE(double *a, double g1, double g2, double g3)\n";
print Matrix "{\n";
print Matrix "   x1 = g1;\n";
print Matrix "   x2 = g2;\n";
print Matrix "   x3 = g3;\n";
print Matrix "\n";
print Matrix "   a[0] = $c[18]";
print Matrix "   a[1] = $c[19]";
print Matrix "   a[2] = $c[20]";
print Matrix "\n";
print Matrix "   return 0;\n";
print Matrix "}\n";
close(Matrix);
