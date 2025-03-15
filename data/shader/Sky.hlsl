#ifndef _SKY_H
#define _SKY_H

#include "AtmosphereCommon.hlsl"

#define DECALRE_CLOUD_FLOAT_PARAMETER(name) static const float name = cGlobal.atmosphere.cloud.name;

static const float cloudBottomRadius = cGlobal.atmosphere.cloud.bottomRadius;
static const float cloudTopRadius = cGlobal.atmosphere.cloud.topRadius;
static const float cloudExtinction = cGlobal.atmosphere.cloud.extinction;
DECALRE_CLOUD_FLOAT_PARAMETER(freq2D0)
DECALRE_CLOUD_FLOAT_PARAMETER(freq2D1)
DECALRE_CLOUD_FLOAT_PARAMETER(freq3D0)
DECALRE_CLOUD_FLOAT_PARAMETER(freq3D1)
static const float cloudCoverage = cGlobal.atmosphere.cloud.coverage;
DECALRE_CLOUD_FLOAT_PARAMETER(cloudBaseHeight)
DECALRE_CLOUD_FLOAT_PARAMETER(lowFrequencyNoise)
DECALRE_CLOUD_FLOAT_PARAMETER(highFrequencyNoise)
DECALRE_CLOUD_FLOAT_PARAMETER(taperPosition)
DECALRE_CLOUD_FLOAT_PARAMETER(taperFloor)
static const int maxBounce = cGlobal.atmosphere.cloud.maxBounce;

#undef DECALRE_CLOUD_PARAMETER

float SampleCloudExtinction(float3 P)
{
    DECLARE_RESOURCE_INDIRECT(Texture2D<float>, noise2D, cGlobal.atmosphere.cloud.noise2DTex);
    DECLARE_RESOURCE_INDIRECT(Texture3D<float>, noise3D, cGlobal.atmosphere.cloud.noise3DTex);
    
    float r = length(P - earth_center);
    float altitude01 = saturate((r - cloudBottomRadius) / (cloudTopRadius - cloudBottomRadius));
    
    float2 uv0 = P.xz * freq2D0;
    float noise2d0 = noise2D.Sample(bilinearRepeatSampler, uv0);
    float2 uv1 = P.xz * freq2D1;
    float noise2d1 = noise2D.Sample(bilinearRepeatSampler, uv1);
    
    float baseDensity = saturate(noise2d0 * noise2d1 - cloudCoverage);
    
    if (baseDensity == 0)
    {
        return 0;
    }
    
    float3 uvw0 = P * freq3D0;
    float cellLowFreq = noise3D.Sample(bilinearRepeatSampler, uvw0);
    float3 uvw1 = P * freq3D1;
    float cellHighFreq = noise3D.Sample(bilinearRepeatSampler, uvw1);
    
    float heightMask = saturate((noise2d0 * 0.45 + 1.0));
    heightMask *= (altitude01 > cloudBaseHeight ? (altitude01 - cloudBaseHeight) / (1 - cloudBaseHeight) : 1.0 - altitude01 / cloudBaseHeight);

    float highFreq = (cellHighFreq - 0.5) * (cellLowFreq * highFrequencyNoise);
    float lowFreq0 = lerp(cellLowFreq, 0.6, heightMask) * lowFrequencyNoise;
    float lowFreq1 = saturate(saturate(heightMask - taperPosition * noise2d1) + taperFloor);
    
    float lowFreq = baseDensity - lowFreq0 * lowFreq1;
    
    float detail = (highFreq + lowFreq) ;
    
    float base0 = heightMask;
    float base1 = saturate(1.0 - heightMask * heightMask);
    float base = base0 * base1 * noise2d0;
    
    float res = saturate((base * detail) / 0.001);
    
    return res * cloudExtinction;
}

float4 RayShellIntersect(float r, float mu, float bottomRadius, float topRadius)
{
    float4 secs = 0;
    float detBottom = r * r * (mu * mu - 1.0) + bottomRadius * bottomRadius;
    float detTop = r * r * (mu * mu - 1.0) + topRadius * topRadius;
    float sqrtDetBottom = sqrt(max(detBottom, 0));
    float sqrtDetTop = sqrt(max(detTop, 0));
    if (r < bottomRadius)
    {
        secs.x = -r * mu + sqrtDetBottom;
        secs.y = -r * mu + sqrtDetTop;
    }
    else if (r < topRadius)
    {
        if (detBottom >= 0.0 && mu < 0.0)
        {
            secs.y = -r * mu - sqrtDetBottom;
            secs.z = -r * mu + sqrtDetBottom;
            secs.w = -r * mu + sqrtDetTop;
        }
        else
        {
            secs.y = -r * mu + sqrtDetTop;
        }
    }
    else
    {
        if (detBottom >= 0.0 && mu < 0.0)
        {
            secs.x = -r * mu - sqrtDetTop;
            secs.y = -r * mu - sqrtDetBottom;
            secs.z = -r * mu + sqrtDetBottom;
            secs.w = -r * mu + sqrtDetTop;
        }
        else if (detTop >= 0.0 && mu < 0.0)
        {
            secs.x = -r * mu - sqrtDetTop;
            secs.y = -r * mu + sqrtDetTop;
        }
    }
    return secs;
}

void TraceCloudSegment(Ray ray, float t0, float t1, float u, inout float opticalLength)
{
    float t = max(t1 - t0, 0);
    if (t < 0.01)
    {
        return;
    }
    // https://research.nvidia.com/publication/2021-06_unbiased-ray-marching-transmittance-estimator
    // biased ray marching
    float c = cloudExtinction * t;
    float M = ceil(pow((0.015 + c) * (0.65 + c) * (60.3 + c), 1.0 / 3.0));
    float dt = t / M;
    ray.Origin += (t0 + u * dt) * ray.Direction;
    for (float i = 0; i < M; ++i)
    {
        opticalLength += SampleCloudExtinction(ray.Origin) * dt;
        ray.Origin += dt * ray.Direction;
        if (opticalLength > -log(1e-10))
        {
            break;
        }
    }
}

float TraceCloudTransmittance(Ray ray, inout RandomGenerater rnd)
{
    float r = length(ray.Origin - earth_center);
    float mu = dot(ray.Direction, normalize(ray.Origin - earth_center));
    float4 secs = RayShellIntersect(r, mu, cloudBottomRadius, cloudTopRadius);
    float opticalLength = 0;
    float u = rnd.Random();
    TraceCloudSegment(ray, secs.x, secs.y, u, opticalLength);
    TraceCloudSegment(ray, secs.z, secs.w, u, opticalLength);
    return exp(-opticalLength);
}

bool FromSpaceIntersectTopAtmosphereBoundary(float r, float mu, out float near_distance)
{
    float discriminant = r * r * (mu * mu - 1.0) + top_radius * top_radius;
    if (mu < 0.0 && discriminant >= 0.0)
    {
        near_distance = ClampDistance(-r * mu - SafeSqrt(discriminant));
        return true;
    }
    near_distance = 0;
    return false;
}

float3 GetSunIrradiance(float3 position, float3 sun_direction, inout RandomGenerater rnd)
{
    DECLARE_RESOURCE_INDIRECT(Texture2D<float4>, transmittanceTex, cGlobal.atmosphere.transmittanceTex);
    
    float r = length(position - earth_center);
    float3 object_up_direction = normalize(position - earth_center);
    float mu_s = dot(sun_direction, object_up_direction);
    
    float3 sun_visibility;
    if (r > top_radius)
    {
        float near_distance;
        if (FromSpaceIntersectTopAtmosphereBoundary(r, mu_s, near_distance))
        {
            position += near_distance * sun_direction;
            r = length(position - earth_center);
            object_up_direction = normalize(position - earth_center);
            mu_s = dot(sun_direction, object_up_direction);
            sun_visibility = GetSunVisibility(transmittanceTex, r, mu_s);
        }
        else
        {
            sun_visibility = 1.0f;
        }
    }
    else
    {
        sun_visibility = GetSunVisibility(transmittanceTex, r, mu_s);
    }
    
    if (Max3(sun_visibility) > 0)
    {
        sun_visibility *= TraceCloudTransmittance(NewRay(position, sun_direction), rnd);
    }

    return solar_illuminance * sun_visibility;
}

float3 GetSkyRadiance(float3 P, float3 V, float3 L, out float3 transmittance)
{
    DECLARE_RESOURCE_INDIRECT(Texture2D<float4>, transmittanceTex, cGlobal.atmosphere.transmittanceTex);
    DECLARE_RESOURCE_INDIRECT(Texture3D<float4>, scatteringTex, cGlobal.atmosphere.scatteringTex);
    DECLARE_RESOURCE_INDIRECT(Texture3D<float4>, singleMieScatteringTex, cGlobal.atmosphere.singleMieScatteringTex);
    return GetSkyRadiance(transmittanceTex, scatteringTex, singleMieScatteringTex,
                P - earth_center, V, L, transmittance);
}

float3 GetSunRadiance(float LdotV)
{
    if (LdotV < cos(sun_angular_radius))
    {
        return 0;
    }
    
    // https://media.contentapi.ea.com/content/dam/eacom/frostbite/files/s2016-pbs-frostbite-sky-clouds-new.pdf
    // B Sun limb darkening astro-physical models
    float3 u = float3(1.0, 1.0, 1.0);
    float3 a = float3(0.397, 0.503, 0.652);

    float sinLV = sqrt(1.0 - LdotV * LdotV);
        
    float centerToEdge = clamp(sinLV / sin(sun_angular_radius), 0.0, 1.0);
    float mu = sqrt(1.0 - centerToEdge * centerToEdge);

    float3 factor = 1.0 - u * (1.0 - pow((float3) mu, a));

    return solar_illuminance / (PI * sun_angular_radius * sun_angular_radius) * factor;
}

struct PathTraceSkyContext
{
    RandomGenerater rnd;
    Ray ray;
    float3 mainLightDir;
    int layer;
    float sigmaMajor;
    float3 res;
    float3 throughput;
};

int GetLayerFromPosition(float3 P)
{
    float r = length(P - earth_center);
    return r < cloudBottomRadius ? 0 : r < cloudTopRadius ? 1 : 2;
}

float GetLayerMajor(int layer)
{
    float3 sigmaR = rayleigh_scattering;
    float3 sigmaMie = (mie_scattering + mie_absorption);
    float3 sigmaOzone = ozone_absorption;
    if (layer > 0)
    {
        float bottomAltitude = (layer == 1 ? cloudBottomRadius : cloudTopRadius) - bottom_radius;
        sigmaR *= GetRayleighDensity(bottomAltitude);
        sigmaMie *= GetMieDensity(bottomAltitude);
        if (layer == 1)
        {
            sigmaMie += cloudExtinction;
        }
    }
    float3 sigma = sigmaR + sigmaMie + sigmaOzone;
    return Max3(sigma);
}

// sign means intersect bottom or top
float ShellIntersect(Ray ray, float3 center, float rBottom, float rTop)
{
    float r = length(ray.Origin - center);
    float mu = dot(ray.Direction, normalize(ray.Origin - center));
    
    float detBottom = r * r * (mu * mu - 1.0) + rBottom * rBottom;
    float detTop = r * r * (mu * mu - 1.0) + rTop * rTop;
    
    bool intersectBottom = detBottom >= 0 && mu < 0;
    //assert(detTop >= 0)
    float dist = -r * mu + (intersectBottom ? -sqrt(detBottom) : sqrt(max(detTop, 0)));
    dist = max(dist, 0);
    if (intersectBottom)
    {
        dist = -dist;
    }
    return dist;
}

#define EVENT_TYPE_MASK 0x3
#define EVENT_TYPE_ABSORPTION 0x0
#define EVENT_TYPE_THROUGH 0x1
#define EVENT_TYPE_SCATTERING 0x2

#define EVENT_SCATTERING_MASK (0x1 << 2)
#define EVENT_SCATTERING_RAYLEIGH (0x0 << 2)
#define EVENT_SCATTERING_MIE (0x1 << 2)

#define EVENT_THROUGH_MASK (0x1 << 2)
#define EVENT_THROUGH_BOTTOM (0x0 << 2)
#define EVENT_THROUGH_TOP (0x1 << 2)

#define CHECK_EVENT(value, type, flag) ((value & EVENT_##type##_MASK) == EVENT_##type##_##flag)

// https://jannovak.info/publications/SDTracking/SDTracking.pdf
int SpectralTracking(inout PathTraceSkyContext ctx, out float t)
{
    float bottomRadius, topRadius;
    if (ctx.layer == 0)
    {
        bottomRadius = bottom_radius;
        topRadius = cloudBottomRadius;
    }
    else if (ctx.layer == 1)
    {
        bottomRadius = cloudBottomRadius;
        topRadius = cloudTopRadius;
    }
    else
    {
        bottomRadius = cloudTopRadius;
        topRadius = top_radius;
    }
    float tMax = ShellIntersect(ctx.ray, earth_center, bottomRadius, topRadius);
    t = 0;
    while (1)
    {
        t += -log(1 - ctx.rnd.Random()) / ctx.sigmaMajor;
        if (t > abs(tMax))
        {
            t = abs(tMax);
            return EVENT_TYPE_THROUGH | (tMax > 0 ? EVENT_THROUGH_TOP : EVENT_THROUGH_BOTTOM);
        }
        float3 P = GetRayPosition(ctx.ray, t);
        float altitude = length(P - earth_center) - bottom_radius;
        
        float3 sigmaSR = GetRayleighDensity(altitude) * rayleigh_scattering;
        float3 sigmaSM = GetMieDensity(altitude) * mie_scattering;
        if (ctx.layer == 1)
        {
            sigmaSM += SampleCloudExtinction(P);
        }
        float3 sigmaA = GetMieDensity(altitude) * mie_absorption
                        + GetOzoneDensity(altitude) * ozone_absorption;
        float3 sigmaNull = max(ctx.sigmaMajor - sigmaSR - sigmaSM - sigmaA, 0);
        
        // average is not needed
        float Pa = dot(sigmaA, ctx.throughput);
        float PsR = dot(sigmaSR, ctx.throughput);
        float PsM = dot(sigmaSM, ctx.throughput);
        float Pn = dot(sigmaNull, ctx.throughput);
        float c = Pa + PsR + PsM + Pn;
        Pa /= c;
        PsR /= c;
        PsM /= c;
        Pn /= c;
        
        float eta = ctx.rnd.Random();

        if (eta < Pa)
        {
            ctx.throughput = 0;
            return EVENT_TYPE_ABSORPTION;
        }
        else if (eta < Pa + PsR)
        {
            ctx.throughput *= sigmaSR / max(ctx.sigmaMajor * PsR, 1e-9);
            return EVENT_TYPE_SCATTERING | EVENT_SCATTERING_RAYLEIGH;
        }
        else if (eta < Pa + PsR + PsM)
        {
            ctx.throughput *= sigmaSM / max(ctx.sigmaMajor * PsM, 1e-9);
            return EVENT_TYPE_SCATTERING | EVENT_SCATTERING_MIE;
        }
        else
        {
            ctx.throughput *= sigmaNull / max(ctx.sigmaMajor * Pn, 1e-9);
        }
    }
    return 0;
}

float HenyeyGreensteinInvertcdf(float eta, float g)
{
    float one_plus_g2 = 1.0 + g * g;
    float one_minus_g2 = 1.0 - g * g;
    float one_over_2g = 0.5 / g;
    float t = (one_minus_g2) / (1.0 - g + 2.0 * g * eta);
    return one_over_2g * (one_plus_g2 - t * t);
}

float3 GenerateRayleighSample(inout RandomGenerater rnd, float3 dir, out float weight)
{
    float3 outDir = UniformSphereSample(rnd);
    float3 value = RayleighPhaseFunction(dot(outDir, dir));
    float3 pdf = 1.0 / (4.0 * PI);
    weight = value / pdf;
    return outDir;
}

float3 GenerateHGSample(inout RandomGenerater rnd, float3 dir, out float3 weight)
{
    float g = mie_phase_g;
    float cosTheta = HenyeyGreensteinInvertcdf(rnd.Random(), g);
    float sinTheta = sqrt(saturate(1.0 - cosTheta * cosTheta));
    float3 t0, t1;
    CreateOrthonormalBasis(dir, t0, t1);
    float phi = 2.0 * PI * rnd.Random();
    weight = 1.0f;
    return float3(sinTheta * sin(phi) * t0 + sinTheta * cos(phi) * t1 + cosTheta * dir);
}

float3 GenerateLambertSample(inout RandomGenerater rnd, float3 N, float3 albedo, out float3 weight)
{
    float sinTheta = sqrt(rnd.Random());
    float cosTheta = sqrt(saturate(1.0 - sinTheta * sinTheta));
    float3 t0, t1;
    CreateOrthonormalBasis(N, t0, t1);
    float phi = 2.0 * PI * rnd.Random();
//    value = albedo * INV_PI * cosTheta;
//    pdf = INV_PI * cosTheta;
    weight = albedo;
    return float3(sinTheta * sin(phi) * t0 + sinTheta * cos(phi) * t1 + cosTheta * N);
}


float3 PathTraceSky(float3 P, float3 L, float3 V, out float3 outTransmittance, inout RandomGenerater rnd)
{
    PathTraceSkyContext ctx;
    ctx.rnd = rnd;
    ctx.ray = NewRayDesc(P, V);
    ctx.mainLightDir = L;
    ctx.layer = GetLayerFromPosition(P);
    ctx.sigmaMajor = GetLayerMajor(ctx.layer);
    
    ctx.res = 0.0f;
    ctx.throughput = 1.0f;
    
    bool hasScattered = false;
    for (int i = 0; i < maxBounce; ++i)
    {
        float t;
        int event = SpectralTracking(ctx, t);
        if (CHECK_EVENT(event, TYPE, THROUGH))
        {
            if (ctx.layer == 2 && CHECK_EVENT(event, THROUGH, TOP))
            {
                break;
            }
            ctx.ray.Origin += ctx.ray.Direction * t;
            if (ctx.layer == 0 && CHECK_EVENT(event, THROUGH, BOTTOM))
            {
                hasScattered = true;
                float3 groundNormal = normalize(ctx.ray.Origin - earth_center);
                float3 lightBsdf = INV_PI * ground_albedo;
                float NdotL = dot(groundNormal, ctx.mainLightDir);
                ctx.res += ctx.throughput * GetSunIrradiance(ctx.ray.Origin, ctx.mainLightDir, ctx.rnd) * (lightBsdf * NdotL);

                float3 weight;
                ctx.ray.Direction = GenerateLambertSample(rnd, groundNormal, ground_albedo, weight);
                ctx.ray.Origin += ctx.ray.Direction * 1.0f;
                ctx.throughput *= weight;
            }
            else
            {
                ctx.layer += CHECK_EVENT(event, THROUGH, TOP) ? 1 : -1;
            }
            ctx.sigmaMajor = GetLayerMajor(ctx.layer);
        }
        else if (CHECK_EVENT(event, TYPE, SCATTERING))
        {
            hasScattered = true;
            ctx.ray.Origin += ctx.ray.Direction * t;
            float nu = dot(ctx.mainLightDir, ctx.ray.Direction);
            float lightBsdf = CHECK_EVENT(event, SCATTERING, RAYLEIGH) ?
                    RayleighPhaseFunction(nu) : HenyeyGreenstein(nu, mie_phase_g);
            ctx.res += ctx.throughput * GetSunIrradiance(ctx.ray.Origin, ctx.mainLightDir, ctx.rnd) * lightBsdf;
            
            float weight;
            if (CHECK_EVENT(event, SCATTERING, RAYLEIGH))
            {
                ctx.ray.Direction = GenerateRayleighSample(ctx.rnd, ctx.ray.Direction, weight);
            }
            else
            {
                ctx.ray.Direction = GenerateHGSample(ctx.rnd, ctx.ray.Direction, weight);
            }
            ctx.throughput *= weight;
        }
        else // if (CHECK_EVENT(event, TYPE, ABSORPTION))
        {
            break;
        }
        
        float p = saturate(Max3(ctx.throughput));
        if (rnd.Random() >= p)
            break;
        ctx.throughput /= p;
    }
    
    outTransmittance = hasScattered ? 0 : ctx.throughput;
    
    rnd = ctx.rnd;
    return ctx.res;
}



#define SKY_TYPE 2

float HeightFogSkyTrans(float3 V)
{
    float fogOpticalDepth = V.y > 0 ? /*H - H * exp2(-inf/H)*/min(1.0 / V.y, 1e2) : /*inf*/1e2;
    float density = 0.3 /* /H */;
    return exp2(-density * fogOpticalDepth);
}

float3 SkyRadiance(float3 P, float3 L, float3 V, out float3 transmittance, inout RandomGenerater rnd)
{
#if SKY_TYPE == 2
    return PathTraceSky(P, L, V, transmittance, rnd);
#elif SKY_TYPE == 1
    return GetSkyRadiance(P, V, L, transmittance);
#else
    float fogTrans = HeightFogSkyTrans(V);
    float3 fakeRayleigh = (1.0 - fogTrans) * float3(0.2, 0.3, 0.8);
    float3 fakeMie = fogTrans * pow(saturate(dot(L, V)) * 0.5 + 0.5, 30) * lerp(float3(1.0, 0.4, 0.2), 1.0, saturate(L.y));
    transmittance = fogTrans;
    return fakeMie + fakeRayleigh * saturate(L.y + 0.2);
#endif
}

float3 SunIrradiance(float3 P, float3 L, inout RandomGenerater rnd)
{
#if SKY_TYPE > 0
    return GetSunIrradiance(P, L, rnd);
#else
    return solar_illuminance * HeightFogSkyTrans(L);
#endif
}

#endif
