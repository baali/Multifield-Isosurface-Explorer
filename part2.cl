#pragma OPENCL EXTENSION cl_khr_fp64 : enable 

#define BINS 100


float4 SortTetVertices ( float4 v[], float fs[])
{
  for (int i = 3; i > 0; i--)
    {
      for (int j = 0; j < i; j++)
	{
	  if (fs[j] > fs[j + 1])
	    {
	      float swapf = fs[j];
	      fs[j] = fs[j + 1];
	      fs[j + 1] = swapf;

	      float4 swapv = v[j];
	      v[j] = v[j + 1];
	      v[j + 1] = swapv;
	    }
	}
    }
}

float4
FindGradient ( float4 *verts, float *fus)
{
  //Sort first
  float4 vn[4];
  float fsn[4];
  for (int i = 0; i < 4; i++)
    {
      vn[i] = verts[i];
      fsn[i] = fus[i];
    }
  // sorting fs and changing values of v correspondingly.
  SortTetVertices (vn, fsn);
  float4 grad = {0, 0, 0, 0};
  float4 base[4];
  int baseCount = 0;
  if (fsn[0] < fsn[3])
    // this condition should always hold true unless the gradient is 0?
    {
      // mean of extreme values?
      // float equality cant be checked?
      float func = (fsn[0] + fsn[3]) / 2;
      
      // It takes care of case when the avg coincides with scalar value at one of vertices
      if (func == fsn[1])
        func = (fsn[0] + fsn[1]) / 2;	
      if (func == fsn[2])
        func = (fsn[0] + fsn[2]) / 2;

      float t;
      // is this calculation of slope or something
      for (int i = 3; i > 0; i--)
      	{
      	  for (int j = 0; j < i; ++j)
      	    // these two loops are for all possible combinations of indexes.
      	    {	
      	      //j & i
      	      if (fsn[j] < func && fsn[i] > func)
      		{
      		  t = (func - fsn[j]) / (fsn[i] - fsn[j]);
      		  if (baseCount < 3)
      		    {
		      float x = (1 - t) * (vn[j].x) + t * (vn[i].x);
	      	      float y = (1 - t) * (vn[j].y) + t * (vn[i].y);
	      	      float z = (1 - t) * (vn[j].z) + t * (vn[i].z);
		      base[baseCount] = (float4)(x, y, z, 1);
      		      baseCount++;		      
      		    }
      		}
	      else if (fsn[j] > func && fsn[i] > func)
		{
		  // Do nothing
		}
	      else if(fsn[i] > func)
		{
		  // Case when func value is equal to one of value at vertices
		  t = (func - fsn[j]) / (fsn[i] - fsn[j]);
      		  if (baseCount < 3)
      		    {
		      float x = (1 - t) * (vn[j].x) + t * (vn[i].x);
	      	      float y = (1 - t) * (vn[j].y) + t * (vn[i].y);
	      	      float z = (1 - t) * (vn[j].z) + t * (vn[i].z);
		      base[baseCount] = (float4)(x, y, z, 1);
      		      baseCount++;		      
      		    }

		}
      	    }
      	}

      if (baseCount < 3)
	return grad;
      float4 v1v2 = base[1] - base[0];
      float4 v1v3 = base[2] - base[0];
      float4 crossp = cross(v1v2, v1v3);
      // We have to check this.
      float norm = sqrt(dot(crossp, crossp));

      if (norm > 0)
	{
	  crossp = normalize(crossp);
	  float4 l = base[0] - vn[0];
	  float dist = dot (l, crossp);
	  if (dist < 0)
	    dist = -dist;
	  //if (dist < 0.000001f)
	  //  {
	  //    return grad;
	  //  }

	  if (dist > 0)
	    {
	      float mag = (func - fsn[0]) / dist;
	      grad = crossp;
	      grad *= mag;
	    }
	}
    }
  return grad;
}

void PlotKappaForTet ( float4 v[4], float fs[4], float gs[4], float hs[4], float kappa,__global float increment, __global float minmax[], __global float *bins,__global float *binst)
{
  SortTetVertices (v, hs);

  float totVolume = fabs(dot((v[0] - v[3]), cross((v[1] - v[3]), (v[2] - v[3]))))/6;

  float4 prevVertices[6];
  int prevVertexCount = 0;

  float4 currentVertices[6];
  int currentVertexCount = 0;

  float prevF = hs[0];
  int index = (hs[0] - minmax[0]) / increment;

  float currentF = minmax[0] + increment * (float) index;
  currentF += increment;

  float4 prevAppex[2];
  float4 currentAppex[2];
  int currentCount, prevCount;
  prevCount = 0;
  float finalVol;
  int incFlag = 0;
  int mark = 0;
  // Removed = from the comparison, I guess float equivalence cant be used on GPU.
  while (currentF < prevF)
    {
      currentF += increment;
      index += 1;
    }
    prevVertices[0] = v[0];
    prevVertexCount = 1;
    prevCount = 8;
    int flag = 0; // for checking the use of prevAppex for the first time.
    for (int m = 1; m < 4; ++m)
      {
	while (currentF < hs[m])
	  {	    
	    currentVertexCount = 0;
	    for (int i = 3; i > 0; i--)
	      {
		for (int j = 0; j < i; ++j)
		  {
		    if (hs[j] <= currentF && hs[i] >= currentF)
		      {
			float t = (currentF - hs[j]) / (hs[i] - hs[j]);
			float x = (1 - t) * (v[j].x) + t * (v[i].x);
			float y = (1 - t) * (v[j].y) + t * (v[i].y);
			float z = (1 - t) * (v[j].z) + t * (v[i].z);
			currentVertices[currentVertexCount] =
			  (float4)(x, y, z, 1);
			currentVertexCount++;
		      }
		  }
	      }
	    if (currentVertexCount >= 3)
	      {
		float4 n;
		n = cross((currentVertices[1] - currentVertices[0]), (currentVertices[2] - currentVertices[0]));
		currentCount = 0;
		for ( int signCount = 0; signCount < 4; signCount++)
		  {
		    float sign = dot(n, (v[signCount] - currentVertices[0]));
		    currentCount = currentCount << 1;
		    if (sign < 0)
		      currentCount += 0;
		    else
		      currentCount += 1;
 		  }

		if (currentCount == 7 || currentCount == 8)
		  currentAppex[0] = v[0];
		else if (currentCount == 11 || currentCount == 4)
		  currentAppex[0] = v[1];
		else if (currentCount == 13 || currentCount == 2)
		  currentAppex[0] = v[2];
		else if (currentCount == 14 || currentCount == 1)
		  currentAppex[0] = v[3];
		else if(currentCount == 3 || currentCount == 12)
		  {
		    currentAppex[0] = v[0];
		    currentAppex[1] = v[1];
		  }
		else if(currentCount == 10 || currentCount == 5)
		  {
		    currentAppex[0] = v[0];
		    currentAppex[1] = v[2];
		  }
		else if(currentCount == 6 || currentCount == 9)
		  {
		    currentAppex[0] = v[0];
		    currentAppex[1] = v[3];
		  }		
	      }

	    if (currentVertexCount == 3 || prevVertexCount == 3)
	      {
		if(prevVertexCount == 1)
		  {
		    finalVol = fabs(dot((currentVertices[0] - prevVertices[0]), cross((currentVertices[1] - prevVertices[0]), (currentVertices[2] - prevVertices[0]))))/6;
		  }

		else if(currentVertexCount == 1)
		  {
		    finalVol = fabs(dot((prevVertices[0] - currentVertices[0]),cross((prevVertices[1] - currentVertices[0]), (prevVertices[2] - currentVertices[0]))))/6;
		  }

		else if(prevVertexCount == 2)
		  {		    
		    mark = 1;
		    float vol = fabs(dot((currentVertices[0] - currentAppex[0]), cross((currentVertices[1] - currentAppex[0]), (currentVertices[2] - currentAppex[0]))))/6;
		    finalVol = totVolume - vol;
		  }

		else if(currentVertexCount == 2)
		  {		    
		    mark = 1;
		    float vol = fabs(dot((prevVertices[0] - prevAppex[0]), cross((prevVertices[1] - prevAppex[0]), (prevVertices[2] - prevAppex[0]))))/6;
		    finalVol = totVolume - vol;
		  }

		else if (currentCount == prevCount && prevVertexCount == currentVertexCount)
		  {
		    float vol1 = fabs(dot((currentVertices[0] - prevVertices[1]), cross((currentVertices[1] - prevVertices[1]), (currentVertices[2] - prevVertices[1]))))/6;
		    float vol2 = fabs(dot((currentVertices[0] - prevVertices[1]), cross((prevVertices[0] - prevVertices[1]), (currentVertices[2] - prevVertices[1]))))/6;
		    float vol3 = fabs(dot((currentVertices[2] - prevVertices[0]), cross((prevVertices[0] - prevVertices[1]), (prevVertices[2] - prevVertices[1]))))/6;

		    finalVol = vol1 + vol2 + vol3;
		  }

		else if(currentVertexCount == prevVertexCount && currentCount != prevCount)
		  {
		    if (flag == 0)
		      flag = 1;

		    float vol1 = fabs(dot((currentVertices[0] - prevVertices[1]), cross((currentVertices[1] - prevVertices[1]), (currentVertices[2] - prevVertices[1]))))/6;
		    float vol2 = fabs(dot((currentVertices[0] - prevVertices[1]), cross((prevVertices[0] - prevVertices[1]), (currentVertices[2] - prevVertices[1]))))/6;
		    float vol3 = fabs(dot((currentVertices[2] - prevVertices[0]), cross((prevVertices[0] - prevVertices[1]), (prevVertices[2] - prevVertices[1]))))/6;
		    finalVol = vol1 + vol2 + vol3;
		  }

		else if(currentVertexCount == 4)
		  {
		    mark = 1;
		    currentAppex[0] = prevAppex[0];
		    int x = prevCount ^ currentCount;
		    if (x == 1)
		      currentAppex[1] = v[0];
		    else if (x == 2)
		      currentAppex[1] = v[1];
		    else if (x == 4)
		      currentAppex[1] = v[2];
		    else if (x == 8)
		      currentAppex[1] = v[3];
		    
		    float volP1 = fabs(dot((currentVertices[0] - currentAppex[0]), cross((currentVertices[1] - currentAppex[0]), (currentVertices[2] - currentAppex[0]))))/6;
		    float volP2 = fabs(dot((currentVertices[0] - currentAppex[0]), cross((currentVertices[2] - currentAppex[0]), (currentVertices[3] - currentAppex[0]))))/6;
		    float volP3 = fabs(dot((currentVertices[0] - currentAppex[1]), cross((currentAppex[0] - currentAppex[1]), (currentVertices[3] - currentAppex[1]))))/6;
		    float volT = fabs(dot((prevVertices[0] - prevAppex[0]), cross((prevVertices[1] - prevAppex[0]), (prevVertices[2] - prevAppex[0]))))/6;
		    finalVol = (volP1 + volP2 + volP3 - volT);
		  }
		
		else if(prevVertexCount == 4)
		  {
		    mark = 1;
		    prevAppex[0] = currentAppex[0];
		    int x = prevCount ^ currentCount;
		    if (x == 1)
		      prevAppex[1] = v[0];
		    else if (x == 2)
		      prevAppex[1] = v[1];
		    else if (x == 4)
		      prevAppex[1] = v[2];
		    else if (x == 8)
		      prevAppex[1] = v[3];
		    float volP1 = fabs(dot((prevVertices[0] - prevAppex[0]), cross((prevVertices[1] - prevAppex[0]), (prevVertices[2] - prevAppex[0]))))/6;
		    float volP2 = fabs(dot((prevVertices[0] - prevAppex[0]), cross((prevVertices[2] - prevAppex[0]), (prevVertices[3] - prevAppex[0]))))/6;
		    float volP3 = fabs(dot((prevVertices[0] - prevAppex[1]), cross((prevAppex[0] - prevAppex[1]), (prevVertices[3] - prevAppex[1]))))/6;
		    float volT = fabs(dot((currentVertices[0] - currentAppex[0]), cross((currentVertices[1] - currentAppex[0]), (currentVertices[2] - currentAppex[0]))))/6;
		    finalVol = (volP1 + volP2 + volP3 - volT);
		  }
		
		else if(prevVertexCount == 5)
		  {
		    if (incFlag == 1)
		      {		
			float volT = fabs(dot((currentVertices[0] - currentAppex[0]), cross((currentVertices[1] - currentAppex[0]), (currentVertices[2] - currentAppex[0]))))/6;
			float volP1 = fabs(dot((prevVertices[3] - currentAppex[0]), cross((prevVertices[2] - currentAppex[0]), (prevVertices[4] - currentAppex[0]))))/6;
			float volP2 = fabs(dot((prevVertices[0] - currentAppex[0]), cross((prevVertices[1] - currentAppex[0]), (prevVertices[3] - currentAppex[0]))))/6;
			float volP3 = fabs(dot((prevVertices[0] - currentAppex[0]), cross((prevVertices[2] - currentAppex[0]), (prevVertices[3] - currentAppex[0]))))/6;
			finalVol = volP1 + volP2 + volP3 - volT;			
			incFlag = 0;
		      }
		    else if (incFlag == 2)
		      {		
			float volT1 = fabs(dot((prevVertices[0] - prevAppex[0]), cross((prevVertices[1] - prevAppex[0]), (prevVertices[2] - prevAppex[0]))))/6;
			float volT2 = fabs(dot((currentVertices[0] - currentAppex[0]), cross((currentVertices[1] - currentAppex[0]), (currentVertices[2] - currentAppex[0]))))/6;
			finalVol = totVolume - volT1 - volT2;			
			incFlag = 0;
		      }
		  }
		  
	      } // end of if condition handling when any of vertex count is 3
	    
	    else if (currentVertexCount == 4 || prevVertexCount == 4)
	      {
		if (currentVertexCount == 1)
		  {
		    mark = 1;
		    float4 n;
		    n = (cross((prevVertices[1] - prevVertices[0]), (prevVertices[2] - prevVertices[0])));
		    int sign = dot(n, (currentVertices[0] - prevVertices[0])) >= 0 ? 1 : 0;
		    if (sign > 0)
		      {
			int count = prevCount;
			int ind = 0;
			for (int counta = 0; counta < 4; counta++)
			  {
			    if ((count & 0x0001) ^ 0x0001 == 0)
			      currentAppex[ind++] = prevVertices[counta];
			    count = count >> 1;
			  }			  
		      }
		    else
		      {
			int count = prevCount;
			int ind = 0;
			for (int counta = 0; counta < 4; counta++)
			  {
			    if ((count & 0x0001) ^ 0x0001 == 1)
			      currentAppex[ind++] = prevVertices[counta];
			    count = count >> 1;
			  }
		      }
		    float volP1 = fabs(dot((prevVertices[0] - currentAppex[0]), cross((prevVertices[1] - currentAppex[0]),(currentAppex[1] - currentAppex[0]))))/6;
		    float volP2 = fabs(dot((prevVertices[0] - currentAppex[1]), cross((prevVertices[1] - currentAppex[1]),(prevVertices[2] - currentAppex[1]))))/6;
		    float volP3 = fabs(dot((prevVertices[0] - currentAppex[1]), cross((prevVertices[1] - currentAppex[1]),(prevVertices[3] - currentAppex[1]))))/6;
		    finalVol = (volP1 + volP2 + volP3);
		  }

		else if (prevVertexCount == 1)
		  {
		    float4 n;
		    mark = 1;
		    n = cross((currentVertices[1] - currentVertices[0]), (currentVertices[2] - currentVertices[0]));
		    int sign = dot(n, (prevVertices[0] - currentVertices[0])) >= 0 ? 1 : 0;
		    if (sign > 0)
		      {
			int count = currentCount;
			int ind = 0;
			for (int counta = 0; counta < 4; counta++)
			  {
			    if ((count & 0x0001) ^ 0x0001 == 0)
			      prevAppex[ind++] = v[3 - counta];
			    count = count >> 1;
			  }			  
		      }
		    else
		      {
			int count = currentCount;
			int ind = 0;
			for (int counta = 0; counta < 4; counta++)
			  {
			    if ((count & 0x0001) ^ 0x0001 == 1)
			      prevAppex[ind++] = v[3 - counta];
			    count = count >> 1;
			  }
		      }
		    float volP1 = fabs(dot((currentVertices[0] - prevAppex[0]), cross((currentVertices[1] - prevAppex[0]), (prevAppex[1] - prevAppex[0]))))/6;
		    float volP2 = fabs(dot((currentVertices[0] - prevAppex[1]), cross((currentVertices[1] - prevAppex[1]), (currentVertices[2] - prevAppex[1]))))/6;
		    float volP3 = fabs(dot((currentVertices[0] - prevAppex[1]), cross((currentVertices[1] - prevAppex[1]), (currentVertices[3] - prevAppex[1]))))/6;
		    finalVol = (volP1 + volP2 + volP3);
		  }
		
		else if(currentVertexCount == 2)
		  {
		    float volP1 = fabs(dot((prevVertices[1] - currentVertices[0]), cross((prevVertices[3] - currentVertices[0]), (currentVertices[1] - currentVertices[0]))))/6;
		    float volP2 = fabs(dot((prevVertices[0] - currentVertices[0]), cross((prevVertices[1] - currentVertices[0]), (prevVertices[3] - currentVertices[0]))))/6;
		    float volP3 = fabs(dot((prevVertices[0] - currentVertices[0]), cross((prevVertices[2] - currentVertices[0]), (prevVertices[3] - currentVertices[0]))))/6;
		    finalVol = (volP1 + volP2 + volP3);
		  }
		
		else if (prevVertexCount == 2)
		  {
		    float volP1 = fabs(dot((currentVertices[1] - prevVertices[0]), cross((currentVertices[3] - prevVertices[0]), (prevVertices[1] - prevVertices[0]))))/6;
		    float volP2 = fabs(dot((currentVertices[0] - prevVertices[0]), cross((currentVertices[1] - prevVertices[0]), (currentVertices[3] - prevVertices[0]))))/6;
		    float volP3 = fabs(dot((currentVertices[0] - prevVertices[0]), cross((currentVertices[2] - prevVertices[0]), (currentVertices[3] - prevVertices[0]))))/6;
		    finalVol = (volP1 + volP2 + volP3);
		  }

		else if (prevVertexCount == currentVertexCount )
		  {
		    if (incFlag != 0)
		      {				
			float volT = fabs(dot((prevVertices[0] - prevAppex[0]), cross((prevVertices[1] - prevAppex[0]), (prevVertices[2] - prevAppex[0]))))/6;
			float volP1 = fabs(dot((currentVertices[3] - prevAppex[0]), cross((currentVertices[1] - prevAppex[0]), (prevVertices[3] - prevAppex[0]))))/6;
			float volP2 = fabs(dot((currentVertices[0] - prevAppex[0]), cross((currentVertices[1] - prevAppex[0]), (currentVertices[3] - prevAppex[0]))))/6;
			float volP3 = fabs(dot((currentVertices[0] - prevAppex[0]), cross((currentVertices[2] - prevAppex[0]), (currentVertices[3] - prevAppex[0]))))/6;
			finalVol = volP1 + volP2 + volP3 - volT;			
			incFlag = 0;
		      }
		    else
		      {
			float volP1 = fabs(dot((prevVertices[3] - currentAppex[0]), cross((prevVertices[1] - currentAppex[0]), (currentAppex[1] - currentAppex[0]))))/6;
			float volP2 = fabs(dot((prevVertices[0] - currentAppex[0]), cross((prevVertices[1] - currentAppex[0]), (prevVertices[3] - currentAppex[0]))))/6;
			float volP3 = fabs(dot((prevVertices[0] - currentAppex[0]), cross((prevVertices[2] - currentAppex[0]), (prevVertices[3] - currentAppex[0]))))/6;
			float vol1 = volP1 + volP2 + volP3;
			volP1 = fabs(dot((currentVertices[3] - currentAppex[0]), cross((currentVertices[1] - currentAppex[0]), (currentAppex[1] - currentAppex[0]))))/6;
			volP2 = fabs(dot((currentVertices[0] - currentAppex[0]), cross((currentVertices[1] - currentAppex[0]), (currentVertices[3] - currentAppex[0]))))/6;
			volP3 = fabs(dot((currentVertices[0] - currentAppex[0]), cross((currentVertices[2] - currentAppex[0]), (currentVertices[3] - currentAppex[0]))))/6;
			float vol2 = volP1 + volP2 + volP3;
			
			finalVol = (vol2 > vol1) ? (vol2 - vol1) : (vol1 - vol2);
		      }
		  }		
	      }

	    bins[index] += finalVol * kappa;
	    binst[index] += finalVol;
	    prevF = currentF;
	    currentF += increment;
	    index++;
	    for (int i = 0; i < currentVertexCount; ++i)
	      {
		prevVertices[i] = currentVertices[i];
	      }
	    prevVertexCount = currentVertexCount;
	    prevCount = currentCount;
	    prevAppex[0] = currentAppex[0];
	    prevAppex[1] = currentAppex[1];

	  }
	if (prevVertexCount >= 3) incFlag += 1;
	prevVertices[prevVertexCount] = v[m];
	prevVertexCount++;
	prevCount |= (0x1 << (4-m));
      }
    if (prevVertexCount == 3)
      {
	mark = 1;
	finalVol = fabs(dot((prevVertices[0] - v[3]), cross((prevVertices[1] - v[3]), (prevVertices[2] - v[3]))))/6;
      }

    else if (prevVertexCount == 4)
      {
	if (flag == 0)
	  {
	    finalVol = fabs(dot((prevVertices[0] - v[3]), cross((prevVertices[1] - v[3]), (prevVertices[2] - v[3]))))/6;
	  }
	else
	  {
	    // float volP1 = fabs(dot((prevVertices[0] - currentAppex[0]), cross((prevVertices[1] - currentAppex[0]), (currentAppex[1] - currentAppex[0])))/6;
	    float volP2 = fabs(dot((prevVertices[0] - currentAppex[0]), cross((prevVertices[1] - currentAppex[0]), (prevVertices[2] - currentAppex[0]))))/6;
	    float volP3 = fabs(dot((prevVertices[0] - currentAppex[0]), cross((prevVertices[1] - currentAppex[0]), (prevVertices[3] - currentAppex[0]))))/6;	
	    finalVol = (volP2 + volP3);
	  } 

      }// done case for vertex count == 4
    else if(prevVertexCount == 6)
      {
	if (currentVertexCount == 3)
	  {
	    float vol = fabs(dot((prevVertices[0] - v[0]), cross((prevVertices[1] - v[0]), (prevVertices[2] - v[0]))))/6;
	    finalVol = totVolume - vol;
	  }
	else
	  {
	    float volP1 = fabs(dot((prevVertices[0] - prevVertices[4]), cross((prevVertices[1] - prevVertices[4]), (prevVertices[5] - prevVertices[4]))))/6;
	    float volP2 = fabs(dot((prevVertices[0] - prevVertices[4]), cross((prevVertices[1] - prevVertices[4]), (prevVertices[3] - prevVertices[4]))))/6;
	    float volP3 = fabs(dot((prevVertices[0] - prevVertices[4]), cross((prevVertices[2] - prevVertices[4]), (prevVertices[3] - prevVertices[4]))))/6;
	    finalVol = (volP1 + volP2 + volP3);
	  }
      }
      // have to do initialisation yaar.
      bins[index] += finalVol * kappa;
      binst[index] += finalVol;
}

__kernel void part2(__global float4 (*vg)[8], __global float (*fsg)[8], __global float (*gsg)[8], __global float (*hsg)[8],__global int *kappaFlag, __global float range[],__global float *inc, __global float (*bins)[110], __global float (*binst)[110])
{
    //get our index in the array
    unsigned int i = get_global_id(0);
    
    float4 v[24];
    float fs[24], gs[24], hs[24];
    int indices[] = {0, 1, 5, 6, 0, 1, 2, 6, 0, 2, 3, 6, 0, 4, 5, 6, 0, 4, 7, 6, 0, 3, 7, 6};
    for (int j = 0; j < 24; j++)
       {
	 v[j] = vg[i][indices[j]];
	 fs[j] = fsg[i][indices[j]];
	 gs[j] = gsg[i][indices[j]];
	 hs[j] = hsg[i][indices[j]];
       }
    for (int p = 0; p < 6; p++)
       {         
         float kappa = 1;
    	 if (*kappaFlag == 1)
       	    {
	      float4 gradF = FindGradient (&v[p*4], &fs[p*4]);
              float4 gradG = FindGradient (&v[p*4], &gs[p*4]);
              float4 crossp = cross (gradF, gradG);
              kappa = sqrt(dot(crossp, crossp));
       	    }    
	 PlotKappaForTet (&v[p*4], &fs[p*4], &gs[p*4], &hs[p*4], kappa, *inc, range, bins[i], binst[i]);
       }
}

__kernel void initial(__global float (*bins)[110], __global float (*binst)[110])
{
    //get our index in the array
    unsigned int i = get_global_id(0);
    for(int j = 0; j < 110; j++)
      {
        bins[i][j] = 0; 
        binst[i][j] = 0; 
      }
}

__kernel void summ(__global float (*binR)[110], __global float (*binstR)[110], __global float bins[110], __global float binst[110], __global int *p)
{
    //get our index in the array
    unsigned int i = get_global_id(0);    
    // bins[i] = binst[i] = 0;
    for(int j = 0; j < *p; j++)
      {
        bins[i] += binR[j][i]; 
        binst[i] += binstR[j][i]; 
      }
}
