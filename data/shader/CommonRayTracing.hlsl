#ifndef _COMMON_RAY_TRACING_H
#define _COMMON_RAY_TRACING_H

#include "ShaderResources.hlsl"

struct HitProperty
{
    float4 albedo;
    float3 normalWS;
    float3 emissive;
    float metallic;
    float perceptualRoughness;
    float ao;
    float alphaCutoff;
    float hitT;
    float roughness;
    float3 diffuse;
    float3 F0;
};

void InterpolateVertex(StructuredBuffer<Vertex> vertexBuffer, uint3 index3, float2 barycentrics, out float2 uv, out float3 normal, out float4 tangent)
{
    Vertex vertex3[3];
    vertex3[0] = vertexBuffer[index3[0]];
    vertex3[1] = vertexBuffer[index3[1]];
    vertex3[2] = vertexBuffer[index3[2]];
    
    uv = vertex3[0].uv + barycentrics.x * (vertex3[1].uv - vertex3[0].uv) + barycentrics.y * (vertex3[2].uv - vertex3[0].uv);
    normal = vertex3[0].normal + barycentrics.x * (vertex3[1].normal - vertex3[0].normal) + barycentrics.y * (vertex3[2].normal - vertex3[0].normal);
    tangent = vertex3[0].tangent + barycentrics.x * (vertex3[1].tangent - vertex3[0].tangent) + barycentrics.y * (vertex3[2].tangent - vertex3[0].tangent);
}

HitProperty CalcHitProperty(RTASRootConstant rtasRoot, uint instanceIndex, uint geometryIndex, uint primitiveIndex, float2 barycentrics, float3x4 worldToObject3x4 = 0)
{
    HitProperty res;
    
    DECLARE_RESOURCE_INDIRECT(StructuredBuffer<Model>, modelBuffer, rtasRoot.modelBuffer);
    DECLARE_RESOURCE_INDIRECT(StructuredBuffer<Material>, materialBuffer, rtasRoot.materialBuffer);
    DECLARE_RESOURCE_INDIRECT(StructuredBuffer<uint>, instanceOffsetBuffer, rtasRoot.instanceOffsetBuffer);
    DECLARE_RESOURCE_INDIRECT(StructuredBuffer<Geometry>, geometryBuffer, rtasRoot.geometryBuffer);
    
    uint instanceOffset = instanceOffsetBuffer[instanceIndex];
    Geometry geometry = geometryBuffer[instanceOffset + geometryIndex];
    Model model = modelBuffer[geometry.modelIndex];
    
    DECLARE_RESOURCE_INDIRECT(StructuredBuffer<uint>, indexBuffer, model.indexBuffer);
    DECLARE_RESOURCE_INDIRECT(StructuredBuffer<Vertex>, vertexBuffer, model.vertexBuffer);
    
    uint3 index3 = geometry.vertexOffset;
    uint baseIndexOffset = geometry.indexOffset + primitiveIndex * 3;
    index3[0] += indexBuffer[baseIndexOffset + 0];
    index3[1] += indexBuffer[baseIndexOffset + 1];
    index3[2] += indexBuffer[baseIndexOffset + 2];
    
    float2 uv;
    float3 normal;
    float4 tangent;
    InterpolateVertex(vertexBuffer, index3, barycentrics, uv, normal, tangent);
    float3 bitangent = normalize(cross(tangent.xyz, normal)) * tangent.w;
    float3x3 TBN = transpose(float3x3(tangent.xyz, bitangent, normal));
    
    Material material = materialBuffer[geometry.materialIndex];
    
    DECLARE_RESOURCE_INDIRECT(Texture2D<float4>, albedoTex, material.baseColorTex);
    DECLARE_RESOURCE_INDIRECT(Texture2D<float3>, normalTex, material.normalTex);
    DECLARE_RESOURCE_INDIRECT(Texture2D<float3>, metallicRoughnessTex, material.metallicRoughnessTex);
    DECLARE_RESOURCE_INDIRECT(Texture2D<float3>, emissiveTex, material.emissiveTex);
    DECLARE_RESOURCE_INDIRECT(Texture2D<float>, aoTex, material.aoTex);
    
    res.albedo = albedoTex.Sample(anisotropicSampler, uv) * float4(material.baseColorFactor.rgb, 1.0);
    float3 normalTS = normalTex.Sample(anisotropicSampler, uv) * 2.0 - 1.0;
    float3 normalOS = mul(TBN, normalTS);
    res.normalWS = normalize(mul(normalOS, (float3x3)worldToObject3x4));
    //res.normalWS = normalize(mul(invModelT, normal));
    //res.normalWS = normalTex.Sample(anisotropicSampler, uv) * 2.0 - 1.0;
    res.emissive = emissiveTex.Sample(anisotropicSampler, uv) * material.emissiveFactor;
    float2 metallicRoughness = metallicRoughnessTex.Sample(anisotropicSampler, uv).bg;
    res.metallic = metallicRoughness.r * material.metallicFactor;
    res.perceptualRoughness = metallicRoughness.g * material.roughnessFactor;
    res.ao = aoTex.Sample(anisotropicSampler, uv) * material.aoFactor;
    res.alphaCutoff = material.alphaCutoff;
    res.hitT = 0;
    
    res.F0 = lerp(0.04, res.albedo.rgb, res.metallic);
    res.diffuse = res.albedo.rgb - res.albedo.rgb * res.metallic;
    res.roughness = res.perceptualRoughness * res.perceptualRoughness;
    
    return res;
}

void ProceedWithAlphaTest(RayQuery< RAY_FLAG_NONE> rq, RTASRootConstant rtasRoot)
{
    while (rq.Proceed())
    {
        if (rq.CandidateType() == CANDIDATE_NON_OPAQUE_TRIANGLE)
        {
            HitProperty property = CalcHitProperty(rtasRoot, rq.CandidateInstanceIndex(), rq.CandidateGeometryIndex(), rq.CandidatePrimitiveIndex(), rq.CandidateTriangleBarycentrics());
            if (property.albedo.a > property.alphaCutoff)
            {
                rq.CommitNonOpaqueTriangleHit();
            }
        }
    }
}

float TraceShadowRay(RTASRootConstant rtasRoot, RayDesc ray)
{
    RayQuery<RAY_FLAG_NONE> rq;
    DECLARE_RESOURCE_INDIRECT(RaytracingAccelerationStructure, tlas, rtasRoot.tlas);
    rq.TraceRayInline(tlas, RAY_FLAG_ACCEPT_FIRST_HIT_AND_END_SEARCH, 0xFF, ray);
    ProceedWithAlphaTest(rq, rtasRoot);
    return rq.CommittedStatus() == COMMITTED_TRIANGLE_HIT ? 0.0 : 1.0;
}

bool TraceRay(RTASRootConstant rtasRoot, RayDesc ray, out HitProperty property)
{
    RayQuery < RAY_FLAG_NONE > rq;
    DECLARE_RESOURCE_INDIRECT(RaytracingAccelerationStructure, tlas, rtasRoot.tlas);
    rq.TraceRayInline(tlas, RAY_FLAG_NONE, 0xFF, ray);
    ProceedWithAlphaTest(rq, rtasRoot);
    if (rq.CommittedStatus() == COMMITTED_TRIANGLE_HIT)
    {
        property = CalcHitProperty(rtasRoot, rq.CommittedInstanceIndex(), rq.CommittedGeometryIndex(), rq.CommittedPrimitiveIndex(), rq.CommittedTriangleBarycentrics()
            , rq.CommittedWorldToObject3x4());
        property.hitT = rq.CommittedRayT();
        if (!rq.CommittedTriangleFrontFace())
        {
            property.normalWS = -property.normalWS;
        }
        return true;
    }
    return false;
}

#endif

