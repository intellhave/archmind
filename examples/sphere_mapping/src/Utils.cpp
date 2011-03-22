/*
  Parallel Computation of Spherical Parameterizations for Mesh Analysis
  Copyright (C) 2011 Athanasiadis Theodoros and Fudos Ioannis

  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.
*/

#include "Utils.h"
#include <fstream>
#include <iostream>

bool loadTextFile(const std::string &strFile, std::string &strText)
{
	using namespace std;

    // Open the file passed in
    ifstream fin(strFile.c_str());

    // Make sure we opened the file correctly
    if(fin == NULL)
	{
		std::cerr << "Failed to open : " << strFile << std::endl;
        return false;
	}

	std::string strLine = "";
    strText = "";

    // Go through and store each line in the text file within a "string" object
    while(getline(fin, strLine))
    {
        strText += "\n" + strLine;
    }

    // Close our file
    fin.close();

    //Success
    return true;
}