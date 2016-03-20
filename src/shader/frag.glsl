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
const float lens_outer_pupil_radius = 26.000000; // scene facing radius in mm
const float lens_inner_pupil_radius = 11.000000; // sensor facing radius in mm
const float lens_length = 106.056999; // overall lens length in mm
const float lens_focal_length = 41.099998; // approximate lens focal length in mm (BFL)
const float lens_aperture_pos = 47.270996; // distance aperture -> outer pupil in mm
const float lens_aperture_housing_radius = 9.000000; // lens housing radius at the aperture
const float lens_outer_pupil_curvature_radius = 85.000000; // radius of curvature of the outer pupil
const float lens_field_of_view = -0.741060; // cosine of the approximate field of view assuming a 35mm image

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
const float out_x =  + 6.37118 *dx + -1.51108 *x + -1.44991 *y*dx*dy + -0.0011872 *x*lens_ipow(y, 2) + 0.0829203 *lens_ipow(x, 2)*dx;
const float out_y =  + 6.33625 *dy + -1.51159 *y + 0.0832582 *lens_ipow(y, 2)*dy + -1.4688 *x*dx*dy + -0.00118749 *lens_ipow(x, 2)*y;
const float out_dx =  + -0.174918 *dx + -0.0828335 *x + 0.00080286 *lens_ipow(y, 2)*dx + 4.52295e-05 *x*lens_ipow(y, 2) + 5.10244e-05 *lens_ipow(x, 3);
const float out_dy =  + -0.1742 *dy + -0.082897 *y + 5.17692e-05 *lens_ipow(y, 3) + 0.000611217 *lens_ipow(x, 2)*dy + 7.6413e-05 *lens_ipow(x, 2)*y;
const float out_transmittance =  + 0.359772  + -0.000299792 *lens_ipow(y, 2) + -0.000424313 *lens_ipow(x, 2) + 0.0958115 *lens_ipow(lambda, 3) + -5.96e-10 *lens_ipow(x, 2)*lens_ipow(y, 6);
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
  pred_x =  + 32.4104 *begin_dx + 0.451791 *begin_x + 0.447088 *begin_x*lens_ipow(begin_dx, 2) + -0.000244526 *begin_x*lens_ipow(begin_y, 2) + -0.000302011 *lens_ipow(begin_x, 3);
  pred_y =  + 32.4102 *begin_dy + 0.451611 *begin_y + 0.449621 *begin_y*lens_ipow(begin_dy, 2) + -0.00030145 *lens_ipow(begin_y, 3) + -0.000240359 *lens_ipow(begin_x, 2)*begin_y;
  pred_dx =  + -0.741653 *begin_dx + -0.0407318 *begin_x + -0.870564 *lens_ipow(begin_dx, 3) + -1.72083e-05 *begin_x*lens_ipow(begin_y, 2) + -1.7307e-05 *lens_ipow(begin_x, 3);
  pred_dy =  + -0.741947 *begin_dy + -0.0407371 *begin_y + -0.859427 *lens_ipow(begin_dy, 3) + -1.72345e-05 *lens_ipow(begin_y, 3) + -1.72529e-05 *lens_ipow(begin_x, 2)*begin_y;
  float dx1_domega0[2][2];
  dx1_domega0[0][0] =  + 32.4104  + 0.894175 *begin_x*begin_dx+0.0f;
  dx1_domega0[0][1] = +0.0f;
  dx1_domega0[1][0] = +0.0f;
  dx1_domega0[1][1] =  + 32.4102  + 0.899242 *begin_y*begin_dy+0.0f;
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

#define NUM_SAMPLES_X 2
#define NUM_SAMPLES_Y 3

#define SWAP(type, a, b){type swp = a; a = b; b = swp;}

vec4 traceRays(vec3 positions[NUM_SAMPLES_X*NUM_SAMPLES_Y+1], vec3 directions[NUM_SAMPLES_X*NUM_SAMPLES_Y+1], float transmittance[NUM_SAMPLES_X*NUM_SAMPLES_Y+1])
{
  vec4 ret = vec4(0);
  //for(int i = 0; i < NUM_SAMPLES_X*NUM_SAMPLES_Y + 1; i++)
  int i = NUM_SAMPLES_X*NUM_SAMPLES_Y;
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
        //color.rgb = textureLod(minmaxdepthmap, r0 * lookupFac + lookupOff, 0).rrr;
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
    ret += color;
  }
  return ret;
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
