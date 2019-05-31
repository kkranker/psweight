//****************************************************************************/
*! $Id$
*! Generalization of IPW and CBPS estimators
*! Compile psweight() class definition
//
*! By Keith Kranker
// Last updated $Date$
//
// Copyright (C) Mathematica Policy Research, Inc.
// This code cannot be copied, distributed or used without the express written
// permission of Mathematica Policy Research, Inc.
//*****************************************************************************/

// psweight.ado relies heavily on a Mata class named psweight().  Stata offers two options for using Mata code:
// (1) re-compiling Mata code each time Stata is launched versus (2) compiling once and saving the reults
// in a .mlib file.  I prefer the latter. Before, using the command for the first time, you wll need to compile
// the Mata code by running this program.
//
// One-time setup instructions:
//
// 1. Open Stata.  Type "adopath" into the command line.
//    It will show you were your PERSONAL directory is located.
//
// 2. Drop these files into your PERSONAL directory (in a subfolder named p/)
//      <PERSONAL>/p/psweight.mata
//      <PERSONAL>/p/psweight.ado
//      <PERSONAL>/p/psweightcall.ado
//
// 3. Update the filepath to psweight.mata below.
//
// 4. Run this program (psweight_one_time_setup.do)
//
// This short file compiles the code and savis it in your PERSONAL director (see help adopath) in a file named ~/l/lpsweight.mlib


clear mata
do "C:\Users\kkranker\Documents\Stata\Ado\Devel\gmatch\psweight.mata"
// do "`c(sysdir_personal)'/p/psweight.mata"
lmbuild lpsweight, dir(PERSONAL) replace
mata: mata describe using lpsweight
mata: mata mlib index
