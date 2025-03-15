#include "CommonRayTracing.hlsl"
#include "Common.hlsl"
#include "Sky.hlsl"

DECLARE_ROOT(PathTracerRoot)

float Pow5(float x)
{
    float x2 = x * x;
    return x2 * x2 * x;
}

float3 F_Schlick(float HdotV, float3 F0)
{
    return F0 + (1.0 - F0) * Pow5(1.0 - HdotV);
}

float D_GGX(float a, float NdotH)
{
    float a2 = a * a;
    a2 = max(a2, 1e-6f);
    float d = (NdotH * a2 - NdotH) * NdotH + 1.0;
    return a2 / max(PI * d * d, 1e-5);
}

float SmithJointApprox_G1_Over4NdotV(float a, float NdotV)
{
    float Vis_SmithV = (NdotV * (1 - a) + a) + NdotV;
    return 0.5 / max(Vis_SmithV, 1e-9);
}

float Vis_SmithJointApprox(float a, float NdotV, float NdotL)
{
    float Vis_SmithV = NdotL * (NdotV * (1 - a) + a);
    float Vis_SmithL = NdotV * (NdotL * (1 - a) + a);
    return 0.5 / max(Vis_SmithV + Vis_SmithL, 1e-9);
}

float3 DirectLightRadiance(RTASRootConstant rtasRoot, float3 pos, float3 L, float3 V, HitProperty hitProperty, inout RandomGenerater rnd)
{
    RayDesc ray = NewRayDesc(pos, L, 0.0001, 10000.0);
    
    float3 N = hitProperty.normalWS;
    float3 H = normalize(L + V);
    
    float NdotL = saturate(dot(N, L));
    float NdotV = saturate(dot(N, V));
    float NdotH = saturate(dot(N, H));
    float HdotV = saturate(dot(H, V));
    
    float3 diffuse = INV_PI * hitProperty.diffuse.rgb;
    float3 F = F_Schlick(HdotV, hitProperty.F0);
    float ggx = D_GGX(hitProperty.roughness, NdotH);
    float vis = Vis_SmithJointApprox(hitProperty.roughness, NdotV, NdotL);
    float3 specular = min(ggx * vis * F, 10.0);
    
    float3 radiance = SunIrradiance(pos, L, rnd) * saturate(dot(L, hitProperty.normalWS)) * (diffuse + specular);
    if (Max3(radiance) > 0)
    {
        float3 t1, t2;
        CreateOrthonormalBasis(ray.Direction, t1, t2);
        
        const float tanRadius = 0.01f;
        float angle = 2.0f * PI * rnd.Random();
        float r = tanRadius * sqrt(rnd.Random());
        ray.Direction = normalize(ray.Direction + r * mul(float2(cos(angle), sin(angle)), float2x3(t1, t2)));
        radiance *= TraceShadowRay(rtasRoot, ray);
    }
    return radiance;
}

float4 CosineSample(float2 E, float3 N)
{
	// https://cseweb.ucsd.edu/~tzli/cse272/wi2023/lectures/28_misc.pdf
    float theta = 2.0 * PI * E.x;
    E.y = 2.0 * E.y - 1.0;
    float3 spherePoint = float3(sqrt(1.0 - E.y * E.y) * float2(cos(theta), sin(theta)), E.y);
    float3 res = normalize(N + spherePoint);
    float pdf = INV_PI * dot(res, N);
    return float4(res, pdf);
}

// https://hal.archives-ouvertes.fr/hal-01509746/document
float3 GgxVndf(float3 wo, float roughness, float u1, float u2)
{
    const float3 YAxis_ = float3(0.0, 1.0, 0.0);
    const float3 XAxis_ = float3(1.0, 0.0, 0.0);
    
    // -- Stretch the view vector so we are sampling as though
    // -- roughness==1
    float3 v = normalize(float3(wo.x * roughness,
                                wo.y,
                                wo.z * roughness));

    // -- Build an orthonormal basis with v, t1, and t2
    float3 t1 = (v.y < 0.999f) ? normalize(cross(v, YAxis_)) : XAxis_;
    float3 t2 = cross(t1, v);

    // -- Choose a point on a disk with each half of the disk weighted
    // -- proportionally to its projection onto direction v
    float a = 1.0f / (1.0f + v.y);
    float r = sqrt(u1);
    float phi = (u2 < a) ? (u2 / a) * PI
                         : PI + (u2 - a) / (1.0f - a) * PI;
    float p1 = r * cos(phi);
    float p2 = r * sin(phi) * ((u2 < a) ? 1.0f : v.y);

    // -- Calculate the normal in this stretched tangent space
    float3 n = p1 * t1 + p2 * t2
             + sqrt(max(0.0f, 1.0f - p1 * p1 - p2 * p2)) * v;

    // -- unstretch and normalize the normal
    return normalize(float3(roughness * n.x,
                            max(0.0f, n.y),
                            roughness * n.z));
}

float3 GgxVndf(float3 wo, float3 N, float roughness, float u1, float u2)
{
    float3 t0, t1;
    CreateOrthonormalBasis(N, t0, t1);
    float3 woTS = mul(float3x3(t0, N, t1), wo);
    float3 wiTS = GgxVndf(woTS, roughness, u1, u2);
    return mul(wiTS, float3x3(t0, N, t1));
}

// mixture sampling
// https://dl.acm.org/doi/pdf/10.1145/3592435
// https://pbr-book.org/4ed/Monte_Carlo_Integration/Improving_Efficiency#MultipleImportanceSampling

// ggx vndf sampling
// https://schuttejoe.github.io/post/ggximportancesamplingpart2/
float3 GenerateNextDir(inout RandomGenerater rnd, float3 V, HitProperty hitProperty, out float3 weight)
{
#if 0
    float3 nextDir = CosineSample(float2(rnd.Random(), rnd.Random()), hitProperty.normalWS).xyz;
    weight = hitProperty.albedo.rgb;
    return nextDir;
#else
    float3 N = hitProperty.normalWS;
    float NdotV = saturate(dot(N, V));
    float specularWeight = Max3(F_Schlick(NdotV, hitProperty.F0));
    float diffuseWeight = Max3(hitProperty.diffuse);
    float P_specular = specularWeight / max(1e-5, specularWeight + diffuseWeight);
    float3 L;
    float c;
    if (rnd.Random() < P_specular)
    {
        float3 H = GgxVndf(V, N, hitProperty.roughness, rnd.Random(), rnd.Random());
        L = reflect(-V, H);
    }
    else
    {
        L = CosineSample(float2(rnd.Random(), rnd.Random()), N).xyz;
    }
    
    float NdotL = saturate(dot(N, L));
    if (NdotL <= 0)
    {
        weight = 0;
        return L;
    }
    
    float3 H = normalize(L + V);
    float NdotH = saturate(dot(N, H));
    float HdotV = saturate(dot(H, V));
            
    float diffusePdf = INV_PI * NdotL;
    float specularPdf = SmithJointApprox_G1_Over4NdotV(
        hitProperty.roughness, NdotV)* D_GGX(hitProperty.roughness, NdotH);
    float pdfMix = lerp(diffusePdf, specularPdf, P_specular);
    
    
    float3 diffuse = INV_PI * hitProperty.diffuse.rgb;
    
    float3 F = F_Schlick(HdotV, hitProperty.F0);
    float ggx = D_GGX(hitProperty.roughness, NdotH);
    float vis = Vis_SmithJointApprox(hitProperty.roughness, NdotV, NdotL);
    float3 specular = ggx * vis * F;

    weight = NdotL * (diffuse + specular) / (pdfMix + 1e-9);
    return L;
#endif
}

void ApplyLens(inout RayDesc ray, float3 frontDir, float lensDistance, float apertureRadius, float focusingDistance, inout RandomGenerater rnd)
{
    float distanceScale = rcp(dot(frontDir, ray.Direction));
    float3 focusingPoint = focusingDistance * distanceScale * ray.Direction;
    float3 lensCenterPoint = lensDistance * distanceScale * ray.Direction;
    
    float3 t1, t2;
    CreateOrthonormalBasis(frontDir, t1, t2);
    float angle = 2.0f * PI * rnd.Random();
    float r = apertureRadius * sqrt(rnd.Random());
    float3 lensPoint = lensCenterPoint + r * mul(float2(cos(angle), sin(angle)), float2x3(t1, t2));
    ray.Origin += lensPoint;
    ray.Direction = normalize(focusingPoint - lensPoint);
}

[numthreads(8,8,1)]
void CS(uint3 ThreadID: SV_DispatchThreadID)
{
    DECLARE_RESOURCE(RWTexture2D<float4>, texRW)
    DECLARE_RESOURCE(ConstantBuffer<PathTracerConstant>, cb)
    
    const float3 camPos = cGlobal.camPos;
    const float3 mainLightDir = cGlobal.mainLightDir;
    const uint2 PixelID = ThreadID.xy + cb.dispatchGroupOffset * 8;
    
    RandomGenerater rnd = NewRandomGenerater(PCGHash(PixelID.x ^ PCGHash(PixelID.y ^ PCGHash(cb.accFrameCount))));
#if REALTIME_PATH_TRACING
    float2 jitter = 0.5;
#else
    float2 jitter = float2(rnd.Random(), rnd.Random());
#endif
    float2 uv = (PixelID.xy + jitter) / float2(cGlobal.targetSize);
    float3 viewDir = normalize(mul(cGlobal.invVP, float4(uv * 2.0 - 1.0, 0.0, 1.0)).xyz);
    
    float3 Lo = 0;
    float3 throughput = 1;
    
    RayDesc ray = NewRayDesc(cGlobal.camPos, viewDir, 0.0001, 10000.0);
    
#if REALTIME_PATH_TRACING
    float depth = 0;
    float3 diffuseAlbedo = 0;
    float3 specularAlbedo = 0;
    float4 normalRoughness = 0;
#endif
    if (cGlobal.EnableDOF)
    {
        ApplyLens(ray, cGlobal.camFront, cGlobal.lensDistance, cGlobal.apertureRadius, cGlobal.focusingDistance, rnd);
    }
    
    for (int i = 0; i < 8; ++i)
    {
        HitProperty hitProperty;
        if (!TraceRay(cb.rtasRoot, ray, hitProperty))
        {
            float3 transmittance;
            float3 skyRadiance = SkyRadiance(ray.Origin, mainLightDir, ray.Direction, transmittance, rnd);
            if (i == 0)
            {
                skyRadiance += transmittance * GetSunRadiance(dot(mainLightDir, ray.Direction));
            }
            Lo += skyRadiance * throughput;
            break;
        }
        if (cb.debugView > 0)
        {
            if (cb.debugView == 1)
                Lo = hitProperty.albedo;
            else if (cb.debugView == 2)
                Lo = float4(hitProperty.normalWS * 0.5 + 0.5, 1);
            else if (cb.debugView == 3)
                Lo = float4(hitProperty.emissive, 1);
            else if (cb.debugView == 4)
                Lo = hitProperty.metallic;
            else if (cb.debugView == 5)
                Lo = hitProperty.perceptualRoughness;
            else if (cb.debugView == 6)
                Lo = hitProperty.ao;
            break;
        }
        
        ray.Origin += hitProperty.hitT * ray.Direction;
        Lo += DirectLightRadiance(cb.rtasRoot, ray.Origin, mainLightDir, -ray.Direction, hitProperty, rnd) * throughput;
        Lo += hitProperty.emissive * throughput;
        
#if REALTIME_PATH_TRACING
        if (i == 0)
        {
            float4 coord = mul(cGlobal.VP, float4(ray.Origin, 1.0));
            depth = coord.z / coord.w;
            diffuseAlbedo = hitProperty.diffuse;
            specularAlbedo = hitProperty.F0;
            normalRoughness = float4(hitProperty.normalWS, hitProperty.roughness);
        }
#endif
            
        float3 weight;
        ray.Direction = GenerateNextDir(rnd, -ray.Direction, hitProperty, weight);
        throughput *= weight;
        
        float p = saturate(Max3(throughput));
        if (rnd.Random() >= p)
            break;
        throughput /= p;
    }
    
    float4 outColor = float4(Lo, 1);
    
#if REALTIME_PATH_TRACING
    DECLARE_RESOURCE(RWTexture2D<float>, depthRW);
    DECLARE_RESOURCE(RWTexture2D<float2>, mvRW);
    DECLARE_RESOURCE(RWTexture2D<float3>, diffuseAlbedoRW);
    DECLARE_RESOURCE(RWTexture2D<float3>, specularAlbedoRW);
    DECLARE_RESOURCE(RWTexture2D<float4>, normalRoughnessRW);
        
    float4 preCoord = mul(cGlobal.PreVP, mul(cGlobal.invVP, float4(uv * 2.0 - 1.0, depth, 1)));
    float2 preUV = preCoord.xy / preCoord.w * 0.5 + 0.5;
    float2 velocity = preUV - uv - cGlobal.jitterMotion * 0.5;
        
    depthRW[PixelID.xy] = depth;
    mvRW[PixelID.xy] = velocity;
    diffuseAlbedoRW[PixelID.xy] = diffuseAlbedo;
    specularAlbedoRW[PixelID.xy] = specularAlbedo;
    normalRoughnessRW[PixelID.xy] = normalRoughness;
    
    texRW[PixelID.xy] = outColor;
#else
    if (cb.accFrameCount <= 1)
    {
        texRW[PixelID.xy] = outColor;
    }
    else
    {
        texRW[PixelID.xy] = lerp(texRW[PixelID.xy], outColor, 1.0 / cb.accFrameCount);
    }
#endif
}

