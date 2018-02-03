/*
Copyright (c) 2018, Fabian Prada
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution.

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/


#include "array-hierarchy.h"


void Array_Hierarchy::UpdateHierachy()
{
	if(active_level == 0){
		printf("Error. Level Cannot be Zero!!\n");
	}
	if(active_level >hierarchy_levels - 1){
		printf("Error. Maximum of Levels Exceded!!\n");
	}

	if(roi_absolute_scale < 0.5f*level_absolute_scale[active_level-1]) // Add a new reference level to the hierarchy
	{
		level_absolute_scale[active_level] = level_absolute_scale[active_level-1]*0.5f;
		level_absolute_coordinate[active_level] = roi_absolute_coordinate-Point<2>(1.f,1.f)*((level_absolute_scale[active_level] - roi_absolute_scale)/2.f);
		updated_level[active_level] = false;
		active_level++;
	}
	else if(roi_absolute_scale > level_absolute_scale[active_level-1]  &&  active_level > 1){
			active_level--;
	}


	level_absolute_scale[active_level] = roi_absolute_scale;
	level_absolute_coordinate[active_level]  = roi_absolute_coordinate;
	updated_level[active_level] =false;

	for(int i=active_level-1 ; i>0; i--) // Correct the levels position
	{
		Point<2> level_coord_diff = (level_absolute_coordinate[i+1] -level_absolute_coordinate[i])/(level_absolute_scale[i]);
		double level_scale_ratio = level_absolute_scale[i+1]/level_absolute_scale[i];
		double level_scale_diff =  level_absolute_scale[i] -level_absolute_scale[i+1];
		if( level_coord_diff[0] < padding_factor_threshold || level_coord_diff[1] < padding_factor_threshold || level_coord_diff[0] > 1.f - padding_factor_threshold - level_scale_ratio || level_coord_diff[1] > 1.f - padding_factor_threshold - level_scale_ratio){
			updated_level[i] =false;
			level_absolute_coordinate[i] = level_absolute_coordinate[i+1]-Point<2>(1.f,1.f)*(level_scale_diff/2.f);
		}
	}

	if (0) {
		for (int i = 0; i <= active_level; i++) {
			printf("Absolute Scale %f \n", level_absolute_scale[i]);
			printf("Absolute Position (%f,%f) \n", level_absolute_coordinate[i][0], level_absolute_coordinate[i][1]);
			printf("Level %d is updated ? %d \n", i, updated_level[i]);
		}
	}

	for(int i=1; i<=active_level; i++) // Update Level
	{
		if(!updated_level[i]){
			UpdateLevel(10,i);
			updated_level[i] = true;
		}
	}
}
