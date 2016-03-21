#version 450
uniform samplerCube cubemap;
uniform samplerCube minmaxdepthmap;
uniform float dist;
uniform float exposure;

in data
{
    vec2 sensorPos;
};

layout(location = 0) out vec4 col;
#define FLT_MAX 1e10
const float lambda = 0.550f;

//init.h
const float lens_outer_pupil_radius = 165.000000; // scene facing radius in mm
const float lens_inner_pupil_radius = 20.000000; // sensor facing radius in mm
const float lens_length = 236.087021; // overall lens length in mm
const float lens_focal_length = 22.787001; // approximate lens focal length in mm (BFL)
const float lens_aperture_pos = 134.100006; // distance aperture -> outer pupil in mm
const float lens_aperture_housing_radius = 15.000000; // lens housing radius at the aperture
const float lens_outer_pupil_curvature_radius = 204.582993; // radius of curvature of the outer pupil
const float lens_field_of_view = -0.845934; // cosine of the approximate field of view assuming a 35mm image

void lens_sphereToCs(vec2 inpos, vec2 indir, out vec3 outpos, out vec3 outdir, float sphereCenter, float sphereRad)
{
  vec3 normal = vec3(inpos/sphereRad, sqrt(max(0.0f, sphereRad*sphereRad-dot(inpos, inpos)))/abs(sphereRad));
  vec3 tempDir = vec3(indir, sqrt(max(0.0f, 1.0f - dot(indir, indir))));

  vec3 ex = normalize(vec3(normal.z, 0, -normal.x));
  vec3 ey = cross(normal, ex);
  
  outdir = tempDir.x * ex + tempDir.y * ey + tempDir.z * normal;
  outpos = vec3(inpos, normal.z * sphereRad + sphereCenter);
}

float lens_ipow(const float a, const int exp)
{
    float ret = 1;
    for(int i = 0; i < exp; i++)
        ret *= a;
    return ret;
}

float eval(vec4 sensor, out vec4 outer)
{
    const float x = sensor.x+dist*sensor.p;
    const float y = sensor.y+dist*sensor.q;
    const float dx = sensor.p;
    const float dy = sensor.q;
//pt_evaluate.h
const float out_x =  + 11.5082 *dx + -7.83728 *x + -533.923 *dx*lens_ipow(dy, 2) + -541.852 *lens_ipow(dx, 3) + 0.000868343 *lens_ipow(x, 2)*lens_ipow(y, 2)*dx;
const float out_y =  + 11.8302 *dy + -7.83372 *y + -544.861 *lens_ipow(dy, 3) + -536.503 *lens_ipow(dx, 2)*dy + 0.000876488 *lens_ipow(x, 2)*lens_ipow(y, 2)*dy;
const float out_dx =  + -0.0644007 *dx + -0.0647349 *x + -3.46031 *dx*lens_ipow(dy, 2) + -2.8413 *lens_ipow(dx, 3) + 3.83374e-05 *lens_ipow(x, 3);
const float out_dy =  + -0.176798 *dy + -0.064167 *y + -1.89333 *lens_ipow(dy, 3) + 4.10222e-05 *lens_ipow(y, 3) + 8.35355e-05 *lens_ipow(x, 2)*y;
const float out_transmittance =  + 0.353663  + -8.85901e-05 *lens_ipow(x, 2) + 0.0898054 *lens_ipow(lambda, 3) + -3.63356e-07 *lens_ipow(y, 4) + -8.90844e-07 *lens_ipow(x, 2)*lens_ipow(y, 2);
    outer = vec4(out_x, out_y, out_dx, out_dy);
    return out_transmittance;
}

void sample_ap(inout vec4 sensor, inout vec4 aperture)
{
float x = sensor.x, y = sensor.y, dx = sensor.z, dy = sensor.w;
float out_x = aperture.x, out_y = aperture.y, out_dx = aperture.z, out_dy = aperture.w;
//pt_sample_aperture.h
float pred_x;
float pred_y;
float pred_dx;
float pred_dy;
float sqr_err = FLT_MAX;
for(int k=0;k<5&&sqr_err > 1e-4f;k++)
{
  const float begin_x = x + dist * dx;
  const float begin_y = y + dist * dy;
  const float begin_dx = dx;
  const float begin_dy = dy;
  const float begin_lambda = lambda;
  pred_x =  + 38.3526 *begin_dx + 0.114175 *begin_x + -14.6495 *lens_ipow(begin_dx, 3) + 0.527942 *begin_y*begin_dx*begin_dy + -0.286615 *begin_x*lens_ipow(begin_dy, 2);
  pred_y =  + 38.357 *begin_dy + 0.114619 *begin_y + -14.602 *lens_ipow(begin_dy, 3) + -0.291044 *begin_y*lens_ipow(begin_dx, 2) + 0.5251 *begin_x*begin_dx*begin_dy;
  pred_dx =  + -0.810404 *begin_dx + -0.0302534 *begin_x + -1.17548 *begin_dx*lens_ipow(begin_dy, 2) + 0.0338555 *begin_x*lens_ipow(begin_dx, 2) + -5.26566e-06 *lens_ipow(begin_x, 3);
  pred_dy =  + -0.809669 *begin_dy + -0.0302208 *begin_y + -1.17342 *lens_ipow(begin_dx, 2)*begin_dy + 0.0342015 *begin_y*lens_ipow(begin_dy, 2) + -5.46446e-06 *lens_ipow(begin_y, 3);
  float dx1_domega0[2][2];
  dx1_domega0[0][0] =  + 38.3526  + -43.9485 *lens_ipow(begin_dx, 2) + 0.527942 *begin_y*begin_dy+0.0f;
  dx1_domega0[0][1] =  + 0.527942 *begin_y*begin_dx + -0.573229 *begin_x*begin_dy+0.0f;
  dx1_domega0[1][0] =  + -0.582087 *begin_y*begin_dx + 0.5251 *begin_x*begin_dy+0.0f;
  dx1_domega0[1][1] =  + 38.357  + -43.8059 *lens_ipow(begin_dy, 2) + 0.5251 *begin_x*begin_dx+0.0f;
  float invJ[2][2];
  const float invdet = 1.0f/(dx1_domega0[0][0]*dx1_domega0[1][1] - dx1_domega0[0][1]*dx1_domega0[1][0]);
  invJ[0][0] =  dx1_domega0[1][1]*invdet;
  invJ[1][1] =  dx1_domega0[0][0]*invdet;
  invJ[0][1] = -dx1_domega0[0][1]*invdet;
  invJ[1][0] = -dx1_domega0[1][0]*invdet;
  const float dx1[2] = {out_x - pred_x, out_y - pred_y};
  for(int i=0;i<2;i++)
  {
    dx += invJ[0][i]*dx1[i];
    dy += invJ[1][i]*dx1[i];
  }
  sqr_err = dx1[0]*dx1[0] + dx1[1]*dx1[1];
}
out_dx = pred_dx;
out_dy = pred_dy;

sensor = vec4(x, y, dx, dy);
aperture = vec4(out_x, out_y, out_dx, out_dy);
}

#define NUM_SAMPLES_X 10
#define NUM_SAMPLES_Y 10

#define SWAP(type, a, b){type swp = a; a = b; b = swp;}

vec4 traceRays(vec3 positions[NUM_SAMPLES_X*NUM_SAMPLES_Y+1], vec3 directions[NUM_SAMPLES_X*NUM_SAMPLES_Y+1], float transmittance[NUM_SAMPLES_X*NUM_SAMPLES_Y+1])
{
  #if 0
  vec4 ret = vec4(0,0,0,1);
  //state of the ray (t_min, t_max, done)
  vec3 raystate[NUM_SAMPLES_X*NUM_SAMPLES_Y+1];
  int faceidx[NUM_SAMPLES_X*NUM_SAMPLES_Y+1];
  vec2 texture_dimensions = textureSize(minmaxdepthmap, 0);

  //initialize raystate
  for(int i = 0; i < NUM_SAMPLES_X*NUM_SAMPLES_Y + 1; i++) raystate[i] = vec3(2*lens_outer_pupil_curvature_radius, FLT_MAX, 0);
  
  //for each face of the cubemap:
  //for(int face = 0; face < 6; face++)
  {
  
    int numLevels = textureQueryLevels(minmaxdepthmap);
    //iterate through mipmap levels containing minimum and maximum depth
    for(int j = 0; j < numLevels; j++)
    {
      vec3 dir = vec3(0);
      for(int i = 0; i < NUM_SAMPLES_X*NUM_SAMPLES_Y + 1; i++) dir += directions[i];
      int face = abs(dir.x)>abs(dir.y)?(abs(dir.x)>abs(dir.z)?0:2):(abs(dir.y)>abs(dir.z)?1:2);
      /*
      for(int i = 0; i < NUM_SAMPLES_X*NUM_SAMPLES_Y + 1; i++)
        faceidx[i] = abs(directions[i].x)>abs(directions[i].y)?(abs(directions[i].x)>abs(directions[i].z)?0:2):(abs(directions[i].y)>abs(directions[i].z)?1:2);
      //renormalize directions
      */
      for(int i = 0; i < NUM_SAMPLES_X*NUM_SAMPLES_Y + 1; i++) directions[i] /= abs(directions[i][face]);
      
      //calculate footprint of ray-bundle
      vec3 tcMax = vec3(-FLT_MAX);
      vec3 tcMin = vec3(FLT_MAX);
      for(int i = 0; i < NUM_SAMPLES_X*NUM_SAMPLES_Y + 1; i++)
      {
        if(raystate[i].z > 0)
          continue;
        //intersections of ray with planes at tmin and tmax
        vec3 rayMin = vec3(positions[i] + directions[i]*(raystate[i].x-positions[i][face]));
        vec3 rayMax = vec3(positions[i] + directions[i]*(raystate[i].y-positions[i][face]));/*
        rayMin /= abs(rayMin[face]);
        rayMax /= abs(rayMin[face]);*/
        tcMin = min(tcMin, rayMin);
        tcMin = min(tcMin, rayMax);
        tcMax = max(tcMax, rayMin);
        tcMax = max(tcMax, rayMax);
      }
      
      //get minimum and maximum depth of cubemap at that point (and mipmap level)
      float lod = log2(length(tcMax-tcMin));
      vec2 minmaxDepth = textureLod(minmaxdepthmap, 0.5f * (tcMin + tcMax), lod).rg;
      vec4 color = textureLod(cubemap, 0.5f * (tcMin + tcMax), lod);
      //return textureLod(minmaxdepthmap, positions[NUM_SAMPLES_X*NUM_SAMPLES_Y], dist);
      //vec4 color = texture(cubemap, 0.5f * (tcMin + tcMax));
      //vec2 minmaxDepth = textureLod(minmaxdepthmap, positions[0] + 1000 * directions[0], 0).rg;
      //calculate intersections of the rays with planes at these depth-values
      //and clamp the rays to this range
      for(int i = 0; i < NUM_SAMPLES_X*NUM_SAMPLES_Y + 1; i++)
      {
        vec2 tNew = max(vec2(0), minmaxDepth - vec2(positions[i][face]));
        raystate[i].x = max(raystate[i].x, tNew.x);
        raystate[i].y = min(raystate[i].y, tNew.y);
        //if the ray hit a scene point, accumulate the color
        if(raystate[i].y <= raystate[i].x && raystate[i].z < 1)
        //if(j == 0)
        {
          ret.rgb += color.rgb * transmittance[i];
          ret.a += 1;
          raystate[i].z = 1;
        }
      }
    }
  }
  return ret;
  #else
  vec4 ret = vec4(0,0,0,1);
  vec3 centerRayPos;
  for(int i = NUM_SAMPLES_X*NUM_SAMPLES_Y; i >= 0; i--)
  //int i = NUM_SAMPLES_X*NUM_SAMPLES_Y;
  {
    vec3 pos = positions[i];
    vec3 dir = directions[i];
    
    float near = 2*lens_outer_pupil_curvature_radius;
    float far = 20000;
    float t[12];
    t[0] = (-pos.x-near)/dir.x;
    t[1] = (-pos.x+near)/dir.x;
    t[2] = (-pos.y-near)/dir.y;
    t[3] = (-pos.y+near)/dir.y;
    t[4] = (-pos.z-near)/dir.z;
    t[5] = (-pos.z+near)/dir.z;
    
    t[6] = (-pos.x-far)/dir.x;
    t[7] = (-pos.x+far)/dir.x;
    t[8] = (-pos.y-far)/dir.y;
    t[9] = (-pos.y+far)/dir.y;
    t[10] = (-pos.z-far)/dir.z;
    t[11] = (-pos.z+far)/dir.z;
    
    int mint = 0;
        for(int i = 1; i < 6; i++)
        if((t[mint] < 0 || t[mint] > t[i]) && t[i] > 0)
            mint = i;

    int maxt = mint;
    for(int i = 0; i < 6; i++)
        if((t[maxt+6] < 0 || t[maxt+6] > t[i+6]) && t[i+6] > 0)
            maxt = i;
    
    int faceIdx = mint>>1;
    int faceIdx2 = maxt>>1;
    int sig = sign(mint&1)*2-1;
    int sig2 = sign(maxt&1)*2-1;
    
    vec3 r0 = (pos + t[mint] * dir)/(near*2)+0.5;
    r0[faceIdx] = 0;
    vec3 rmax = (pos + t[mint+6] * dir)/(far*2)+0.5;
    rmax[faceIdx] = sig;
    vec3 rd = rmax - r0;
    rd = normalize(rd);
    
    vec3 lookupOff = vec3(-1);
    lookupOff[faceIdx] = sig;
    
    vec3 lookupFac = vec3(2);
    lookupFac[faceIdx] = 0;
    
    //--------------------------------------------------------
    
    vec3 r02 = (pos + t[maxt] * dir)/(near*2)+0.5;
    r02[faceIdx2] = 0;
    vec3 rmax2 = (pos + t[maxt+6] * dir)/(far*2)+0.5;
    rmax2[faceIdx2] = sig2;
    vec3 rd2 = rmax2 - r02;
    rd2 = normalize(rd2);
    
    vec3 lookupOff2 = vec3(-1);
    lookupOff2[faceIdx2] = sig2;
    
    vec3 lookupFac2 = vec3(2);
    lookupFac2[faceIdx2] = 0;
    
    const int maxLod = 9;
    int lod = maxLod;
    vec4 color = vec4(1);
    // factor for size of one texel (1 = whole image in 1 texel, 
    // 1/2 = 2 texels for whole image, ...) texelSize = 1/(2^(10-lod))
    int numIterations = 0;
    
    for(int i = 0; i < 20; i++,numIterations++)
    {
        if(r0[(faceIdx+1)%3] > 1 || r0[(faceIdx+2)%3] > 1 || r0[(faceIdx+1)%3] < 0 || r0[(faceIdx+2)%3] < 0)
        {
            if(faceIdx == faceIdx2)
                break;
            float prevDepth = r0[faceIdx];
            SWAP(int, faceIdx, faceIdx2);
            SWAP(int, sig, sig2);
            SWAP(vec3, lookupOff, lookupOff2);
            SWAP(vec3, lookupFac, lookupFac2);
            SWAP(vec3, rd, rd2);
            SWAP(vec3, r0, r02);
            
            float t = (prevDepth - abs(r0[faceIdx])) / abs(rd[faceIdx]);
            r0 += t * rd;
            
            lod = maxLod;
        }
        
        float texelSize = 1.0f/(1<<(9-lod));
        
        vec2 minmaxDepth = (textureLod(minmaxdepthmap, r0 * lookupFac + lookupOff, lod).rg)/(far);
        //TODO dependency on dimensions of ray bundle
        float lod = 0;
        color.rgb = textureLod(cubemap, r0 * lookupFac + lookupOff, lod).rgb;
        if(i == NUM_SAMPLES_X*NUM_SAMPLES_Y)
          centerRayPos = r0/abs(r0[faceIdx]);
        else
        {
          lod = length(centerRayPos - r0/abs(r0[faceIdx]))*5;
          //color.rgb = vec3(lod);
        }
        
        color.rgb = textureLod(cubemap, r0 * lookupFac + lookupOff, lod).rgb;
        float dist = -1;
        
        // if the current position is before the min-plane, choose max step size s.t. we reach the min-plane,
        // if current ray position is between min and max depth, look at higher resolution depth (i.e. lod = lod-1)
        // otherwise (current depth greater than the maximum) go to next texel but choose lower resolution for next iteration
        // in any case the maximum step size is also limited by the texel size.
        if(abs(r0[faceIdx]) < minmaxDepth.r)
            dist = (minmaxDepth.r - abs(r0[faceIdx])) / abs(rd[faceIdx]);
        else if(/*abs(r0[faceIdx]) >= minmaxDepth.r*1000 && */abs(r0[faceIdx]) < minmaxDepth.g)
        {
            //color.r = i/20.0f;
            //color.gb = vec2(0);
            if(minmaxDepth.g - minmaxDepth.r < 1e-2)
                break;
            lod = max(lod-1, 0);numIterations--;
            continue;
        }
        else
        {
            lod = min(lod+1, maxLod);numIterations--;
            continue;
        }
        
        vec2 inTexelPos = mod(vec2(r0[(faceIdx+1)%3], r0[(faceIdx+2)%3]), texelSize);
        vec2 inTexelDir = vec2(rd[(faceIdx+1)%3], rd[(faceIdx+2)%3]);
        
        vec2 texelBorder = vec2(sign(inTexelDir.x)*texelSize, sign(inTexelDir.y)*texelSize);
        texelBorder = clamp(texelBorder, 0, texelSize);
        vec2 texelIntersect = (texelBorder-inTexelPos) / (inTexelDir);
        
        r0 += dist*(1+1e-7) * rd;
    }
    color.rgb *= transmittance[i];
    //color.rgb = vec3(r0);
    if(i < NUM_SAMPLES_X*NUM_SAMPLES_Y)
      ret += color;
  }
  return ret;
  #endif
}

//copied from image_head.cl
int clip_aperture(const float x, const float y, const float radius)
{
    // early out
    //if(x*x + y*y > radius*radius) return 0;
    const int num_blades = 8;
    float xx = radius;
    float yy = 0.0f;
    for(int b=1;b<num_blades+1;b++)
    {
        float tmpx, tmpy;
        tmpy = sin(2.0f*3.14159/num_blades * b);
        tmpx = cos(2.0f*3.14159/num_blades * b);
        tmpx *= radius;
        tmpy *= radius;
        const float normalx = xx + tmpx;
        const float normaly = yy + tmpy;
        float dot0 = (normalx)*(x-xx) + (normaly)*(y-yy);
        if(dot0 > 0.0f) return 0;
        xx = tmpx;
        yy = tmpy;
    }
    return 1;
}

void main()
{    
    vec4 center = vec4(sensorPos, -sensorPos/(lens_length-lens_aperture_pos+dist));
    vec4 aperture = vec4(0);
    sample_ap(center, aperture);
    
    vec4 outer = vec4(0);
    float t = eval(center, outer);
    
    vec3 p, d;
    lens_sphereToCs(outer.xy, outer.zw, p, d, 0, lens_outer_pupil_curvature_radius);
    if(length(outer.xy) > lens_outer_pupil_radius)
        t = 0;
    vec4 color;
    vec3 ray_positions[NUM_SAMPLES_X*NUM_SAMPLES_Y+1];
    vec3 ray_directions[NUM_SAMPLES_X*NUM_SAMPLES_Y+1];
    float transmittance[NUM_SAMPLES_X*NUM_SAMPLES_Y+1];
    transmittance[NUM_SAMPLES_X*NUM_SAMPLES_Y] = t;
    ray_positions[NUM_SAMPLES_X*NUM_SAMPLES_Y] = p;
    ray_directions[NUM_SAMPLES_X*NUM_SAMPLES_Y] = d;
    
    for(int i = 0; i < NUM_SAMPLES_X; i++)
        for(int j = 0; j < NUM_SAMPLES_Y; j++)
        {
            vec4 ray = center;
            //sample aperture (as disk)
            vec4 aperture = vec4(sin(6.28f*j/NUM_SAMPLES_Y), cos(6.28f*j/NUM_SAMPLES_Y), 0, 0)*sqrt(lens_aperture_housing_radius*i/NUM_SAMPLES_X);
            sample_ap(ray, aperture);
            vec4 outer;
            float t = eval(ray, outer);
            lens_sphereToCs(outer.xy, outer.zw, p, d, 0, lens_outer_pupil_curvature_radius);
            if(clip_aperture(aperture.x, aperture.y, lens_aperture_housing_radius) == 0 || length(outer.xy) > lens_outer_pupil_radius)
                t = 0;
            transmittance[i*NUM_SAMPLES_Y+j] = t;
            ray_positions[i*NUM_SAMPLES_Y+j] = p;
            ray_directions[i*NUM_SAMPLES_Y+j] = d;
        }
    
    color = traceRays(ray_positions, ray_directions, transmittance);
    col = color/color.a/exposure;
}
