#ifndef _SHADER_RESOURCES_H
#define _SHADER_RESOURCES_H

#ifdef __cplusplus

namespace ShaderResource
{
using uint = uint32_t;
using uint2 = Math::uint2;
using uint3 = Math::uint3;
using uint4 = Math::uint4;
using int2 = Math::int2;
using int3 = Math::int3;
using int4 = Math::int4;
using float2 = Math::float2;
using float3 = Math::float3;
using float4 = Math::float4;
using float4x4 = Math::float4x4;

#else

SamplerState anisotropicSampler : register(s0);
SamplerState bilinearSampler : register(s1);
SamplerState bilinearRepeatSampler : register(s2);

#define DECLARE_ROOT(name) ConstantBuffer<name> cRoot : register(b0);

#define DECLARE_RESOURCE(type, name) type name = ResourceDescriptorHeap[cRoot.name];
#define DECLARE_RESOURCE_INDIRECT(type, name, index) type name = ResourceDescriptorHeap[NonUniformResourceIndex(index)];

#endif

#define SCATTERING_TEXTURE_R_SIZE 32
#define SCATTERING_TEXTURE_MU_SIZE 128
#define SCATTERING_TEXTURE_MU_S_SIZE 32
#define SCATTERING_TEXTURE_NU_SIZE 8

#define SCATTERING_TEXTURE_WIDTH (SCATTERING_TEXTURE_NU_SIZE * SCATTERING_TEXTURE_MU_S_SIZE)
#define SCATTERING_TEXTURE_HEIGHT SCATTERING_TEXTURE_MU_SIZE
#define SCATTERING_TEXTURE_DEPTH SCATTERING_TEXTURE_R_SIZE

struct CloudConstant
{
    float bottomRadius;
    float topRadius;
    float extinction;
    float freq2D0;

    float freq2D1;
    float freq3D0;
    float freq3D1;
    float coverage;
    
    float cloudBaseHeight;
    float lowFrequencyNoise;
    float highFrequencyNoise;
    float taperPosition;
    
    float taperFloor;
    int maxBounce;
    uint noise2DTex;
    uint noise3DTex;
};

struct AtmosphereConstant
{
    float3 solarIlluminance;
    float sunAngularRadius;

    float3 rayleighScattering;
    float invRayleighExponentialDistribution;

    float3 mieScattering;
    float invMieExponentialDistribution;

    float3 mieAbsorption;
    float ozoneCenterAltitude;

    float3 ozoneAbsorption;
    float invOzoneWidth;

    float3 groundAlbedo;
    float miePhaseG;

    float bottomRadius;
    float topRadius;
    float transmittanceSteps;
    float multiscatteringSteps;

    uint transmittanceTex;
    uint multiscatteringTex;
    uint scatteringTex;
    uint singleMieScatteringTex;

    CloudConstant cloud;
};

struct GlobalConstant
{
    float4x4 VP;
    float4x4 invVP;
    float4x4 PreVP;
    float3 camPos;
    float lensDistance;
    float3 mainLightDir;
    float apertureRadius;
    float3 camFront;
    float focusingDistance;
    uint2 targetSize;
    uint2 screenSize;
    float2 jitterMotion;
    uint EnableDOF;
    float pad0;
    AtmosphereConstant atmosphere;
};

#ifndef __cplusplus
ConstantBuffer<GlobalConstant> cGlobal : register(b1);
#endif

struct Material
{
    float4 baseColorFactor;
    float3 emissiveFactor;
    float metallicFactor;
    float roughnessFactor;
    float aoFactor;
    uint baseColorTex;
    uint normalTex;
    uint metallicRoughnessTex;
    uint emissiveTex;
    uint aoTex;
    float alphaCutoff;
};

struct Geometry
{
    uint modelIndex;
    int materialIndex;
    uint vertexOffset;
    uint indexOffset;
};

struct Vertex
{
    float3 position;
    float3 normal;
    float2 uv;
    float4 tangent;
};

struct Model
{
    uint indexBuffer;
    uint vertexBuffer;
    uint pad0;
    uint pad1;
};

struct RTASRootConstant
{
    uint tlas;
    uint modelBuffer;
    uint materialBuffer;
    uint instanceOffsetBuffer;
    uint geometryBuffer;
    uint pad0;
    uint pad1;
    uint pad2;
};

struct BlitRoot
{
    uint tex;
    float exposure;
};

struct PathTracerRoot
{
    uint cb;
    uint texRW;
    uint depthRW;
    uint mvRW;
    uint diffuseAlbedoRW;
    uint specularAlbedoRW;
    uint normalRoughnessRW;
};

struct PathTracerConstant
{
    RTASRootConstant rtasRoot;
    uint accFrameCount;
    int debugView;
    uint2 dispatchGroupOffset;
};

struct FocusingRoot
{
    RTASRootConstant rtasRoot;
    int2 focusingPos;
    uint focusingInfoBufferRW;
};

struct AtmosphereRoot
{
    uint transmittanceRW;
    uint multiscatteringRW;
    uint scatteringRW;
    uint singleMieScatteringRW;
};

struct NoiseRoot
{
    uint noise2DRW;
    uint noise3DRW;
};

#ifdef __cplusplus
} // namespace ShaderResource
#endif

#endif
