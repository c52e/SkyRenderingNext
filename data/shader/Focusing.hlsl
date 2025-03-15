#include "CommonRayTracing.hlsl"

DECLARE_ROOT(FocusingRoot)

[numthreads(1,1,1)]
void CS(uint3 ThreadID: SV_DispatchThreadID)
{
    DECLARE_RESOURCE(RWStructuredBuffer<float>, focusingInfoBufferRW)
    
    const float3 camPos = cGlobal.camPos;
    float2 uv = (cRoot.focusingPos + 0.5f) / float2(cGlobal.screenSize);
    uv.y = 1.0 - uv.y;
    float3 viewDir = normalize(mul(cGlobal.invVP, float4(uv * 2.0 - 1.0, 0.0, 1.0)).xyz);
    
    RayDesc ray;
    ray.Origin = cGlobal.camPos;
    ray.Direction = viewDir;
    ray.TMin = 0.1;
    ray.TMax = 1000;
    
    HitProperty hitProperty;
    hitProperty.hitT = 10000.0;
    TraceRay(cRoot.rtasRoot, ray, hitProperty);
    
    float focusingDistance = hitProperty.hitT * dot(viewDir, cGlobal.camFront);
    
    focusingInfoBufferRW[0] = focusingDistance;
}

