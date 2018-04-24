// gmatch relies heavily on a Mata class
// before using the  command for the first time, you need to compile the Mata code.
// This short file compiles the code and saves it in your PERSONAL directory (see help adopath)
// in a file named ~/lgmatch.mlib


//****************************************************************************/
*! $Id: gmatch_one_time_setup.do,v e38398bcb0ac 2018/02/14 08:16:27 kkranker $
*! Generalization of IPW and CBPS estimators
*! Compile gmatch() class definition
//
*! By Keith Kranker
// Last updated $Date: 2018/02/14 08:16:27 $
//
// Copyright (C) Mathematica Policy Research, Inc.
// This code cannot be copied, distributed or used without the express written
// permission of Mathematica Policy Research, Inc.
//*****************************************************************************/

// gmatch.ado relies heavily on a Mata class named gmatch().  Stata offers two options for using Mata code:
// (1) re-compiling Mata code each time Stata is launched versus (2) compiling once and saving the reults
// in a .mlib file.  I prefer the later. Before, using the command for the first time, you will need to compile
// the Mata code by running this program.
//
// One-time setup instructions:
//
// 1. Open Stata.  Type "adopath" into the command line.
//    It will show you were your PERSONAL directory is located.
//
// 2. Drop these files into your PERSONAL directory (in a subfolder named /g/)
//      <PERSONAL>/g/gmatchclass.mata
//      <PERSONAL>/g/gmatch.ado
//      <PERSONAL>/g/gmatchcall.ado
//      <PERSONAL>/l/lmbuild.ado
//
// 3. Update the filepath to gmatchclass.mata below.
//
// 4. Run this program (gmatch_one_time_setup.do)
//
// 5. (Possibly) cd out of the default directory Stata opens in
//    If your version opens in something like C:\Program Files (x86)\Stata15
//    And you get an error like "file lST_1c84_000001_tmp.mlib could not be opened", try cd-ing

clear mata
do "C:\ado\personal\g\gmatchclass.mata"
lmbuild lgmatch, dir(PERSONAL) replace
mata: mata describe using lgmatch
mata: mata mlib index
