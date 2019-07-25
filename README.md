# Poincare-Gauge-Cosmology

(* Title: PGC, Poincare Gauge Cosmology *)

(* Author: Hongchao Zhang *)

(* Summary: Symbolic computing package on Poincare Gauge Cosmology *)

(* Copyright (C) 2018-2019 Hongchao Zhang *)

(* Email: zhanghc@mail.dlut.edu.cn *)

(* Affiliation: Dalian University of Technology & Penn State University *)

(* Package Version: 1.2.1 *)

(* Brief Discussion: 
   - Define a 4D metric-affinely connected manifold.
   - Define a Chart and a time-like observer.
   - Define the components of metric and torsion tensor on the FLRW cosmology
   - Vary the given gravitational Lagrangian with respect the metric and torsion, to get the Einstein and the Cartan field equation, respectively
   - Calculate the cosmological equations on the FLRW background.
   - Export and store the results in files in the working directory.
*)

(* Users are required to cite the PGC papers: *)

(* arXiv: 1904.03545, arXiv: 1906.04340 *)

(* Acknowledgement: This package is based on the Wolfram Mathematica and the xAct series.
   For more information about the xAct, please transfer to the site: http://www.xact.es/
*)

**************************************
INSTALLATION NOTES FOR THE PGC PACKAGE
**************************************

Prerequisites
-------------

* Wolfram Mathematica.

* xAct packages: http://www.xact.es/

When uncompressed, the archive files give a number of files hanging from a directory called PGC121/. This directory must be placed at (or linked from) one of the places Mathematica prepares for external applications. You can find the actual paths in your Mathematica installation in the variables $BaseDirectory and $UserBaseDirectory.

Linux:

   - single-user installation:

        $HOME/.Mathematica/Applications/

Mac OS:

   - single-user installation:

        /Users/<user>/Library/Mathematica/Applications/

MSWindows:

   - single-user installation:

        C:\Documents and settings\<user>\Application Data\Mathematica\Applications\

   Beware that in Windows these directories might be hidden!

For more detail about how to install and load Mathematica packages, please read the "install" file of xAct. In the file "PGC121_test_0.nb" attached to PGC121 archive, you can find some examples and usages.