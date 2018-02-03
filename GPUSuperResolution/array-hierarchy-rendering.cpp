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

void Array_Hierarchy::drawOpenGL(Scene * scene){

	UpdateHierachy();

	for(int i=0; i<active_level; i++){
		reference_canvas->DrawBox(level_absolute_coordinate[i],level_absolute_scale[i]);
	}

	SetTextureToLevel(active_level,texture_resource);
	canvas->DrawTexture(texture_buffer);

	if(vector_field_mode == ADVECTION_VECTOR_FIELD || vector_field_mode == MEAN_GRADIENT_FIELD){
		SampleField_CUDA(vector_field_mode,HYBRID_DOMAIN,(float)roi_absolute_coordinate[0],(float)roi_absolute_coordinate[1],true,(float)roi_absolute_scale,MAXIMA_NORMALIZATION);
		GLTools::TransformFloat2ArrayToGLVectorBuffer_CUDA(vector_field_resource,float2_buffer,width,height,false,true,0.1f,20.f);
		canvas->DrawTriangles(vector_field_buffer,height*width);
	}
}

void Array_Hierarchy::drawROI(const Point<2> p_roi_absolute_coordinate, const double p_roi_absolute_scale)
{
	roi_absolute_coordinate = p_roi_absolute_coordinate;
	roi_absolute_scale = p_roi_absolute_scale;
	drawOpenGL();
}

void Array_Hierarchy::HandleKeyboardEvent(Prompt  & prompt, Scene * scene)
{
	switch (prompt.active_command){
	case 'v':
		ChangeVectorField_CallBack(prompt);
	break;
	}
}

void Array_Hierarchy::setupOpenGL()
{
	Initialize_CUDA_Advection_Object();

	cudaChannelFormatDesc channelDesc2 = cudaCreateChannelDesc(32,32,32,32,cudaChannelFormatKindFloat); 
	hierarchical_color_field_array = new cudaArray*[hierarchy_levels];
	for(int i=1; i<hierarchy_levels; i++){
		cudaMallocArray(&hierarchical_color_field_array[i],&channelDesc2,width,height);
	}

	hierarchical_color_field_array[0] = color_field_array;
	updated_level[0] = true;	

}