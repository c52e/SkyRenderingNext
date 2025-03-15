#include "ShaderResources.hlsl"
#include "Common.hlsl"

DECLARE_ROOT(BlitRoot)

struct VertexOut
{
	float4 pos : SV_POSITION;
    float2 uv : COLOR;
};

VertexOut VS(uint vertex: SV_VertexID)
{
	VertexOut vout;

    float2 uvs[3] =
    {
        { 0, 2 },
        { 2, 0 },
        { 0, 0 }
    };
    float4 positions[3] =
    {
        { -1, 3, 0, 1 },
        { 3, -1, 0, 1 },
        { -1, -1, 0, 1 }
    };
	
	vout.pos = positions[vertex];
    vout.uv = uvs[vertex];
    
    float targetScreenXRatio = float(cGlobal.targetSize.x) * float(cGlobal.screenSize.y) 
                            / (float(cGlobal.screenSize.x) * float(cGlobal.targetSize.y));
    if (targetScreenXRatio > 1)
    {
        vout.uv.y = vout.uv.y * targetScreenXRatio - 0.5 * targetScreenXRatio + 0.5;
    }
    else
    {
        vout.uv.x = vout.uv.x / targetScreenXRatio - 0.5 / targetScreenXRatio + 0.5;
    }
    
    return vout;
}

float3 ToneMapping(float3 luminance, float exposure)
{
	// https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
    const float A = 2.51f;
    const float B = 0.03f;
    const float C = 2.43f;
    const float D = 0.59f;
    const float E = 0.14f;

    luminance *= exposure;
    return (luminance * (A * luminance + B)) / (luminance * (C * luminance + D) + E);
}

float4 PS(VertexOut pin) : SV_Target
{
    if (any(saturate(pin.uv) != pin.uv))
    {
        return 0;
    }
    DECLARE_RESOURCE(Texture2D<float4>, tex);
    
    float4 color = tex.Sample(bilinearSampler, pin.uv);
    color.rgb = ToneMapping(color.rgb, cRoot.exposure);
    
    color.rgb = Linear2sRGB(color.rgb);
    
    RandomGenerater rnd = NewRandomGenerater(PCGHash(int(pin.pos.x) ^ PCGHash(pin.pos.y)));
    color.rgb += rnd.Random() / 255.0;
    
    return color;
}
