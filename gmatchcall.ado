//****************************************************************************/
*! $Id$
*! Generalization of IPW and CBPS estimators
*! Stata command to call class functions on an instance of the class
*! named gmatch_ado_most_recent
//
// After calling gmatch, the data is stored in a class instance named gmatch_ado_most_recent
// You can print any of the public functions or variables to the screen with gmatchcall. For example:
//
*! By Keith Kranker
// Last updated $Date$
//
// Copyright (C) Mathematica Policy Research, Inc.
// This code cannot be copied, distributed or used without the express written
// permission of Mathematica Policy Research, Inc.
//*****************************************************************************/

program define gmatchcall, rclass
  mata: gmatch_ado_most_recent.`0'
  return add
end

