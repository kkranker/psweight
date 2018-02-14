// gmatch relies heavily on a Mata class
// before using the  command for the first time, you need to compile the Mata code.
// This short file compiles the code and savis it in your PERSONAL director (see help adopath)
// in a file named ~/l/lgmatch.mlib


//****************************************************************************/
*! $Id$
*! Generalization of IPW and CBPS estimators
*! Compile gmatch() class definition
//
*! By Keith Kranker
// Last updated $Date$
//
// Copyright (C) Mathematica Policy Research, Inc.
// This code cannot be copied, distributed or used without the express written
// permission of Mathematica Policy Research, Inc.
//*****************************************************************************/

// gmatch.ado relies heavily on a Mata class named gmatch().  Stata offers two options for using Mata code:
// (1) re-compiling Mata code each time Stata is launched versus (2) compiling once and saving the reults
// in a .mlib file.  I prefer the later. Before, using the command for the first time, you wll need to compile 
// the Mata code by running this program.
//
// One-time setup instructions:
// 
// 1. Open Stata.  Type "adopath" into the command line. 
//    It will show you were your PERSONAL directory is located.
// 
// 2. Drop these files into your PERSONAL directory (in a subfolder named /g/)
//      <<PERSONAL>>/g/gmatchclass.mata
//      <<PERSONAL>>/g/gmatch.ado
//      
// 3. Update the filepath to gmatchclass.mata below.
// 
// 4. Run this program (gmatch_one_time_setup.do)

clear mata
do "C:\Users\kkranker\Documents\Stata\Ado\Devel\gmatch\gmatchclass.mata"
lmbuild lgmatch, dir(PERSONAL) replace
mata: mata describe using lgmatch
mata: mata mlib index
