#version 450
uniform samplerCube cubemap;
uniform samplerCube minmaxdepthmap;
uniform float dist;
uniform float exposure;
uniform mat4 invprojectionmatrix;

float lambda = .55f;

in data
{
    vec2 sensorPos;
};

layout(location = 0) out vec4 col;

#define NUM_SAMPLES_X 0
#define NUM_SAMPLES_Y 0
#define SWAP(type, a, b){type swp = a; a = b; b = swp;}

// a few static defines about general lens geometry
#pragma include "init.h"

float lens_ipow(const float x, const int exp)
{
  return pow(x, exp);
}

void lens_sphereToCs(vec2 inpos, vec2 indir, out vec3 outpos, out vec3 outdir, float sphereCenter, float sphereRad)
{
  vec3 normal = vec3(inpos/sphereRad, sqrt(max(0, sphereRad*sphereRad-length(inpos))/abs(sphereRad)));
  vec3 tempDir = vec3(indir, sqrt(max(0, 1-length(indir))));

  vec3 ex = vec3(normal.z, 0, -normal.x);
  ex = normalize(ex);
  vec3 ey = cross(normal, ex);
  
  outdir = tempDir.x * ex + tempDir.y * ey + tempDir.z * normal;
  outpos = vec3(inpos, normal.z * sphereRad + sphereCenter);
}

void lens_csToSphere(vec3 inpos, vec3 indir, out vec2 outpos, out vec2 outdir, float sphereCenter, float sphereRad)
{
  vec3 normal = vec3(inpos.xy, abs((inpos.z-sphereCenter)/sphereRad));
  vec3 tempDir = normalize(indir);

  vec3 ex = vec3(normal.z, 0, -normal.x);
  ex = normalize(ex);
  vec3 ey = cross(normal, ex);
  outpos = inpos.xy;
  outdir = vec2(dot(tempDir, ex), dot(tempDir, ey));
}

// evaluates from sensor (in) to outer pupil (out).
// input arrays are 5d [x,y,dx,dy,lambda] where dx and dy are the direction in
// two-plane parametrization (that is the third component of the direction would be 1.0).
// units are millimeters for lengths and micrometers for the wavelength (so visible light is about 0.4--0.7)
// returns the transmittance computed from the polynomial.
float lens_evaluate(vec4 sensor, out vec4 outer_pupil, float lambda, float dist)
{
  float x = sensor[0], y = sensor[1], dx = sensor[2], dy = sensor[3];
#pragma include "pt_evaluate.h"
  outer_pupil[0] = out_x; outer_pupil[1] = out_y; outer_pupil[2] = out_dx; outer_pupil[3] = out_dy;
  return max(0.0f, out_transmittance);
}

// solves for the two directions [dx,dy], keeps the two positions [x,y] and the
// wavelength, such that the path through the lens system will be valid, i.e.
// lens_evaluate_aperture(in, out) will yield the same out given the solved for in.
// in: point on sensor. out: point on aperture.
float lens_pt_sample_aperture(inout vec4 sensor, inout vec4 aperture, float lambda, float dist)
{
  float out_x = aperture[0], out_y = aperture[1], out_dx = aperture[2], out_dy = aperture[3], out_transmittance = 1.0f;
  float x = sensor[0], y = sensor[1], dx = sensor[2], dy = sensor[3];
#pragma include "pt_sample_aperture.h"
  // directions may have changed, copy all to be sure.
  aperture[0] = out_x; aperture[1] = out_y; aperture[2] = out_dx; aperture[3] = out_dy;
  sensor[0] = x; sensor[1] = y; sensor[2] = dx; sensor[3] = dy;
  return max(0.0f, out_transmittance);
}

vec4 traceRays(vec3 positions[NUM_SAMPLES_X*NUM_SAMPLES_Y+1], vec3 directions[NUM_SAMPLES_X*NUM_SAMPLES_Y+1])
{
    vec3 pos = positions[NUM_SAMPLES_X*NUM_SAMPLES_Y];
    vec3 dir = directions[NUM_SAMPLES_X*NUM_SAMPLES_Y];
    
    float near = 2*lens_outer_pupil_curvature_radius;
    float far = 20000;
    float t[12];
    t[0] = (-pos.x-near)/dir.x;
    t[1] = (-pos.x+near)/dir.x;
    t[2] = (-pos.y-near)/dir.y;
    t[3] = (-pos.y+near)/dir.y;
    t[4] = (-pos.z-near)/dir.z;
    t[5] = (-pos.z+near)/dir.z;
    
    //return vec4(color);
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
    
    mat4 proj = mat4(vec4(1,0,0,0), vec4(0,1,0,0), vec4(0,0,far/(far-near),1), vec4(0,0,-far*near/(far-near),0));

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
        color.rgb = textureLod(cubemap, r0 * lookupFac + lookupOff, 0).rgb;
        color.rgb = textureLod(minmaxdepthmap, r0 * lookupFac + lookupOff, 0).rrr;
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
    color.rgba = vec4(0);
    float depth = r0[faceIdx]*far+(1-r0[faceIdx])*near;
    for(int i = 0; i < NUM_SAMPLES_X*NUM_SAMPLES_Y; i++)
    {
        //positions[i] += depth * directions[i];
        if(any(greaterThan(abs(positions[i]), vec3(0))) || any(greaterThan(abs(directions[i]), vec3(0))))
        {
            color.rgba += vec4(texture(cubemap, positions[i] + depth * directions[i]).rgb, 1);
        }
    }
    
    return color;
}

void main()
{    
    vec4 sensor_center = vec4(sensorPos, -sensorPos/(lens_length-lens_aperture_pos+dist));
    vec4 aperture_center, out_center;
    lens_pt_sample_aperture(sensor_center, aperture_center, lambda, dist);
    float transmittance = lens_evaluate(sensor_center, out_center, lambda, dist);
    vec3 pos, dir;
    lens_sphereToCs(out_center.xy, out_center.zw, pos, dir, lens_length-lens_outer_pupil_curvature_radius, lens_outer_pupil_curvature_radius);
    
    vec4 color;
    vec3 ray_positions[NUM_SAMPLES_X*NUM_SAMPLES_Y+1];
    vec3 ray_directions[NUM_SAMPLES_X*NUM_SAMPLES_Y+1];
    ray_positions[NUM_SAMPLES_X*NUM_SAMPLES_Y] = pos;
    ray_directions[NUM_SAMPLES_X*NUM_SAMPLES_Y] = dir;
    
    for(int i = 0; i < NUM_SAMPLES_X; i++)
        for(int j = 0; j < NUM_SAMPLES_Y; j++)
        {
            vec4 sensor_pos = sensor_center;
            vec4 aperture_pos = vec4(sin(1.0f*j/NUM_SAMPLES_Y), cos(1.0f*j/NUM_SAMPLES_Y), 0, 0)*sqrt(lens_aperture_housing_radius*i/NUM_SAMPLES_X);
            vec4 out_pos;
            lens_pt_sample_aperture(sensor_pos, aperture_pos, lambda, dist);
            
            float transmittance = lens_evaluate(sensor_pos, out_pos, lambda, dist);
            lens_sphereToCs(out_pos.xy, out_pos.zw, ray_positions[i*NUM_SAMPLES_Y+j], ray_directions[i*NUM_SAMPLES_Y+j], lens_length-lens_outer_pupil_curvature_radius, lens_outer_pupil_curvature_radius);
        }
    color = traceRays(ray_positions, ray_directions);
    col = color/color.a/exposure;
}
