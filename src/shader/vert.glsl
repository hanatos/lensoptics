#version 450
#if 0
#define WIDTH (1080)
#define HEIGHT (720)
#endif
uniform float minDepth;

layout(location = 0) in vec2 vert;
layout(location = 1) in vec2 tc;

#if 0
layout(location = 0) out data
{
    //vec2 tcout;
    vec2 apCoords;
};

void main()
{
    vec2 size = vec2(2);//vec2(36/float(WIDTH), 24/float(HEIGHT));
    //tcout = tc*size - size/2;
    apCoords = tc * size - size / 2;
    //vec2 dim = vec2(WIDTH/minDepth, HEIGHT/minDepth);
    vec2 dim = vec2(WIDTH, HEIGHT);
    gl_Position = vec4(
        (vert.x+1)/dim.x-1 + 2*(gl_InstanceID%int(dim.x))/dim.x,
        (vert.y+1)/dim.y-1 + 2*(gl_InstanceID/dim.x)/dim.y, 
        0.0f, 1.0f);
}
#else
layout(location = 0) out data
{
    vec2 sensorPos;
};

void main()
{
    sensorPos = (tc - 0.5f) * vec2(36, 24);
    gl_Position = vec4(vert, 0, 1);
}
#endif
