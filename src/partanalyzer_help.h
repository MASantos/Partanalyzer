/*
*    This file is part of Partanalyzer.
*
*    Partanalyzer is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    Partanalyzer is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with Partanalyzer.  If not, see <http://www.gnu.org/licenses/>.
*
*/
/** Partanalyzer
Copyright (C) Miguel A. Santos, HSC, Toronto, 2008-2011.
Licensed under GPL version 3 a later. (see http://www.gnu.org/copyleft/gpl.html )
*/

#ifndef _HELPFUNCTIONS_H
#define _HELPFUNCTIONS_H 1

#include "partanalyzer.h"
#include "partanalyzer_includes.h"

void printCopyright();

void printVersion();

void printCommandLineError(const string label="");

void printCommandLineError(char* lastSeeOption);

void systemDate();

void printHelp();

void printHelpLong();

void exitWithHelp();

#endif //END _HELPFUNCTIONS_H
