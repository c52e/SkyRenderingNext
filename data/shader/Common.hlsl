#ifndef _COMMON_H
#define _COMMON_H

#define PI 3.1415926535897932384626433832795
#define INV_PI (1.0 / PI)

uint WangHash(uint seed)
{
    seed = (seed ^ 61) ^ (seed >> 16);
    seed *= 9;
    seed = seed ^ (seed >> 4);
    seed *= 0x27d4eb2d;
    seed = seed ^ (seed >> 15);
    return seed;
}

// https://www.reedbeta.com/blog/hash-functions-for-gpu-rendering/
uint PCGHash(uint seed)
{
    uint state = seed * 747796405u + 2891336453u;
    uint word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
    return (word >> 22u) ^ word;
}

struct RandomGenerater
{
    uint seed;
    
    float Random()
    {
        float res = asfloat(seed >> 9 | 0x3f800000u) - 1.0f;
        seed = PCGHash(seed);
        return res;
    }
};

RandomGenerater NewRandomGenerater(uint seed)
{
    RandomGenerater rnd;
    rnd.seed = seed;
    return rnd;
}

float Linear2sRGB(float C_lin)
{
    float C_srgb;
    if (C_lin <= 0.0031308)
        C_srgb = C_lin * 12.92;
    else
        C_srgb = 1.055 * pow(C_lin, 1.0 / 2.4) - 0.055;
    return C_srgb;
}

float3 Linear2sRGB(float3 C_lin)
{
    return float3(Linear2sRGB(C_lin.x), Linear2sRGB(C_lin.y), Linear2sRGB(C_lin.z));
}

float Max3(float3 v)
{
    return max(v.x, max(v.y, v.z));
}

void CreateOrthonormalBasis(float3 N, out float3 t0, out float3 t1)
{
    // https://graphics.pixar.com/library/OrthonormalB/paper.pdf
    float s = (N.z >= 0.0 ? 1.0 : -1.0);
    float a = -1.0 / (s + N.z);
    float b = N.x * N.y * a;
    t0 = float3(1.0 + s * N.x * N.x * a, s * b, -s * N.x);
    t1 = float3(b, s + N.y * N.y * a, -N.y);
}

float3 UniformSphereSample(inout RandomGenerater rnd)
{
    float phi = 2.0 * PI * rnd.Random();
    float cosTheta = 1.0 - 2.0 * rnd.Random();
    float sinTheta = sqrt(clamp(1.0 - cosTheta * cosTheta, 0, 1));
    return float3(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta);
}

RayDesc NewRayDesc(float3 Origin, float3 Direction, float TMin = 0.0001, float TMax = 10000.0)
{
    RayDesc ray;
    ray.Origin = Origin;
    ray.Direction = Direction;
    ray.TMin = TMin;
    ray.TMax = TMax;
    return ray;
}

float3 GetRayPosition(RayDesc ray, float t)
{
    return ray.Origin + ray.Direction * t;
}

#define Ray RayDesc
#define NewRay NewRayDesc

#endif

