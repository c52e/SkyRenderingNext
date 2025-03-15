#include "ShaderResources.hlsl"
#include "Common.hlsl"

DECLARE_ROOT(NoiseRoot)

struct NoiseCreateInfo
{
    uint seed;
    uint baseFrequency;
    float remapMin;
    float remapMax;
};

static const NoiseCreateInfo noise2DCreateInfo =
{
    0, 10, 1.0f, 0.062f
};

static const NoiseCreateInfo noise3DCreateInfo =
{
    0, 4, 0.039f, 1.0f
};

static const float3 kPerlinGradients[16] =
{
    { 1, 1, 0 },
    { -1, 1, 0 },
    { 1, -1, 0 },
    { -1, -1, 0 },
    { 1, 0, 1 },
    { -1, 0, 1 },
    { 1, 0, -1 },
    { -1, 0, -1 },
    { 0, 1, 1 },
    { 0, -1, 1 },
    { 0, 1, -1 },
    { 0, -1, -1 },
    { 1, 1, 0 },
    { -1, 1, 0 },
    { 0, -1, 1 },
    { 0, -1, -1 }
};

float3 GetPerlinGradients(uint i, uint j, uint k, uint seed)
{
    return kPerlinGradients[WangHash(seed + WangHash(i + WangHash(j + WangHash(k)))) & 0xf];
}

// https://www.semanticscholar.org/paper/Improving-noise-Perlin/a6fd5071b73f542c79bd08d409c5f73de38dac5d
// Tiled Perlin Noise 3d
float PerlinNoise(float3 p, uint freq, uint seed)
{
    p *= float(freq);
    uint3 ijk0 = uint3(int3(floor(p))) % freq; // Casting a negative float to uint is undefined
    uint3 ijk1 = uint3(int3(ceil(p))) % freq;
    uint i0 = ijk0.x;
    uint j0 = ijk0.y;
    uint k0 = ijk0.z;
    uint i1 = ijk1.x;
    uint j1 = ijk1.y;
    uint k1 = ijk1.z;

    float3 t = p - floor(p);
    float3 uvw = t * t * t * (t * (t * 6.0f - 15.0f) + 10.0f);
    float u = uvw.x;
    float v = uvw.y;
    float w = uvw.z;

    float x0 = t.x;
    float y0 = t.y;
    float z0 = t.z;
    float x1 = t.x - 1.0f;
    float y1 = t.y - 1.0f;
    float z1 = t.z - 1.0f;

    return lerp(lerp(lerp(dot(GetPerlinGradients(i0, j0, k0, seed), float3(x0, y0, z0)),
                          dot(GetPerlinGradients(i1, j0, k0, seed), float3(x1, y0, z0)), u),
                     lerp(dot(GetPerlinGradients(i0, j1, k0, seed), float3(x0, y1, z0)),
                          dot(GetPerlinGradients(i1, j1, k0, seed), float3(x1, y1, z0)), u), v),
                lerp(lerp(dot(GetPerlinGradients(i0, j0, k1, seed), float3(x0, y0, z1)),
                          dot(GetPerlinGradients(i1, j0, k1, seed), float3(x1, y0, z1)), u),
                     lerp(dot(GetPerlinGradients(i0, j1, k1, seed), float3(x0, y1, z1)),
                          dot(GetPerlinGradients(i1, j1, k1, seed), float3(x1, y1, z1)), u), v), w);
}

// Tiled Worley Noise 3d
float WorleyNoise(float3 p, uint freq, uint seed)
{
    p *= float(freq);
    uint3 ijk = uint3(floor(p));
    p += float(freq);
    float min_dist = 1e10;
    for (uint di = freq - 1; di <= freq + 1; ++di)
    {
        for (uint dj = freq - 1; dj <= freq + 1; ++dj)
        {
            for (uint dk = freq - 1; dk <= freq + 1; ++dk)
            {
                uint3 grid = ijk + uint3(di, dj, dk);
                uint3 gridseed = grid % freq;
                uint rnd0 = WangHash(seed + WangHash(gridseed.x + WangHash(gridseed.y + WangHash(gridseed.z))));
                uint rnd1 = WangHash(seed + WangHash(gridseed.x + WangHash(gridseed.y + WangHash(gridseed.z + 1))));
                uint rnd2 = WangHash(seed + WangHash(gridseed.x + WangHash(gridseed.y + WangHash(gridseed.z + 2))));
                float3 g = float3(grid) + float3(rnd0, rnd1, rnd2) / 4294967296.0f;
                float dist = distance(p, g);
                min_dist = min(min_dist, dist);
            }
        }
    }
    return min_dist;
}

float RemapTo01(float x, float x0, float x1)
{
    return saturate((x - x0) / (x1 - x0));
}

float PerlinFBM(float3 p, NoiseCreateInfo createInfo)
{
    float res = 0.0;
    uint f = createInfo.baseFrequency;
    float a = 0.5;
    float sum_a = 0.0;
    for (uint c = 0; c < 8; ++c)
    {
        float noise = PerlinNoise(p, f, createInfo.seed) * 0.5 + 0.5;
        res += RemapTo01(noise, createInfo.remapMin, createInfo.remapMax) * a;
        sum_a += a;
        f *= 2;
        a *= 0.5;
    }
    res /= sum_a;
    return res;
}

float WorleyFBM(float3 p, NoiseCreateInfo createInfo)
{
    float res = 0.0;
    uint f = createInfo.baseFrequency;
    float a = 0.5;
    float sum_a = 0.0;
    for (uint c = 0; c < 8; ++c)
    {
        float noise = WorleyNoise(p, f, createInfo.seed);
        res += RemapTo01(noise, createInfo.remapMin, createInfo.remapMax) * a;
        sum_a += a;
        f *= 2;
        a *= 0.5;
    }
    res /= sum_a;
    return res;
}

[numthreads(8, 8, 1)]
void CSNoise2D(int3 ThreadID : SV_DispatchThreadID)
{
    DECLARE_RESOURCE(RWTexture2D<float>, noise2DRW)
    
    int2 dim;
    noise2DRW.GetDimensions(dim.x, dim.y);
    
    float2 coord = (float2(ThreadID.xy) + 0.5) / float2(dim);
    float density = PerlinFBM(float3(coord, 0.0), noise2DCreateInfo);

    noise2DRW[ThreadID.xy] = density;
}


[numthreads(4, 4, 4)]
void CSNoise3D(uint3 ThreadID : SV_DispatchThreadID)
{
    DECLARE_RESOURCE(RWTexture3D<float>, noise3DRW)
    
    int3 dim;
    noise3DRW.GetDimensions(dim.x, dim.y, dim.z);
    
    float3 coord = (float3(ThreadID) + 0.5) / float3(dim);
    float worley = WorleyFBM(coord, noise3DCreateInfo);
    noise3DRW[ThreadID] = worley;
}
