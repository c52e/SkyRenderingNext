#include "Atmosphere.h"

#include "Graphics/Shader.h"

static const Math::uint2 kTransmittanceTextureDim{256, 64};
static const Math::uint2 kMultiscatteringTextureDim{32, 32};

static const Math::uint2 kNoise2DSize{512, 512};
static const Math::uint3 kNoise3DSize{64, 64, 64};

Atmosphere::Atmosphere()
{
    m_TransimittanceTexture =
        RenderTexture(kTransmittanceTextureDim, DXGI_FORMAT_R16G16B16A16_FLOAT, true, L"Transmittance LUT");
    m_MultiscatteringTexture =
        RenderTexture(kMultiscatteringTextureDim, DXGI_FORMAT_R16G16B16A16_FLOAT, true, L"MultiScattering LUT");
    m_ScatteringTexture = RenderTexture({SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT, SCATTERING_TEXTURE_DEPTH},
                                        DXGI_FORMAT_R16G16B16A16_FLOAT, true, L"Scattering LUT");
    m_SingleMieScatteringTexture =
        RenderTexture({SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT, SCATTERING_TEXTURE_DEPTH},
                      DXGI_FORMAT_R16G16B16A16_FLOAT, true, L"Single Mie Scattering LUT");

    CompileComputePSO(m_TransmittancePSO, SHADER_PATH("AtmosphereLUT.hlsl"), L"CSTransmittance");
    CompileComputePSO(m_MultiscatteringPSO, SHADER_PATH("AtmosphereLUT.hlsl"), L"CSMultiscattering");
    CompileComputePSO(m_ScatteringPSO, SHADER_PATH("AtmosphereLUT.hlsl"), L"CSScattering");

    m_Noise2D = RenderTexture(kNoise2DSize, DXGI_FORMAT_R8_UNORM, true, L"Noise 2D");
    m_Noise3D = RenderTexture(kNoise3DSize, DXGI_FORMAT_R8_UNORM, true, L"Noise 3D");

    CompileComputePSO(m_Noise2DPSO, SHADER_PATH("NoiseGen.hlsl"), L"CSNoise2D");
    CompileComputePSO(m_Noise3DPSO, SHADER_PATH("NoiseGen.hlsl"), L"CSNoise3D");
}

ShaderResource::AtmosphereConstant Atmosphere::UpdateParameters(const AtmosphereParameters &parameters)
{
    ShaderResource::AtmosphereConstant res{};

    res.solarIlluminance = parameters.solarIlluminance;
    res.sunAngularRadius = Math::radians(parameters.sunAngularRadius);

    res.rayleighScattering = parameters.rayleighScattering * parameters.rayleighScatteringScale;
    res.mieScattering = parameters.mieScattering * parameters.mieScatteringScale;
    res.mieAbsorption = parameters.mieAbsorption * parameters.mieAbsorptionScale;
    res.ozoneAbsorption = parameters.ozoneAbsorption * parameters.ozoneAbsorptionScale;

    res.invRayleighExponentialDistribution = 1.0f / parameters.rayleighExponentialDistribution;
    res.invMieExponentialDistribution = 1.0f / parameters.mieExponentialDistribution;
    res.ozoneCenterAltitude = parameters.ozoneCenterAltitude;
    res.invOzoneWidth = 1.0f / parameters.ozoneWidth;
    res.groundAlbedo = parameters.groundAlbedo;
    res.bottomRadius = parameters.bottomRadius;
    res.topRadius = parameters.bottomRadius + parameters.thickness;
    res.miePhaseG = parameters.miePhaseG;

    res.transmittanceSteps = parameters.transmittanceSteps;
    res.multiscatteringSteps = parameters.multiscatteringSteps;

    res.transmittanceTex = m_TransimittanceTexture.GetSrvIndex();
    res.multiscatteringTex = m_MultiscatteringTexture.GetSrvIndex();
    res.scatteringTex = m_ScatteringTexture.GetSrvIndex();
    res.singleMieScatteringTex = m_SingleMieScatteringTexture.GetSrvIndex();

    res.cloud.bottomRadius = parameters.bottomRadius + parameters.cloud.bottomAltitude;
    res.cloud.topRadius = parameters.bottomRadius + parameters.cloud.bottomAltitude + parameters.cloud.thickness;
    res.cloud.extinction = parameters.cloud.extinction;
    res.cloud.freq2D0 = 1.0 / parameters.cloud.baseWidth2D;
    res.cloud.freq2D1 = 1.0 / (parameters.cloud.baseWidth2D * sqrt(2));
    res.cloud.freq3D0 = 1.0 / parameters.cloud.baseWidth3D;
    res.cloud.freq3D1 = 1.0 / (parameters.cloud.baseWidth3D / sqrt(47));
    res.cloud.coverage = parameters.cloud.coverage;
    res.cloud.cloudBaseHeight = parameters.cloud.cloudBaseHeight;
    res.cloud.lowFrequencyNoise = parameters.cloud.lowFrequencyNoise;
    res.cloud.highFrequencyNoise = parameters.cloud.highFrequencyNoise;
    res.cloud.taperPosition = parameters.cloud.taperPosition;
    res.cloud.taperFloor = parameters.cloud.taperFloor;
    res.cloud.maxBounce = parameters.cloud.maxBounce;

    res.cloud.noise2DTex = m_Noise2D.GetSrvIndex();
    res.cloud.noise3DTex = m_Noise3D.GetSrvIndex();

    return res;
}

void Atmosphere::UpdateLuts(GraphicsCommandList *cmd)
{
    ShaderResource::AtmosphereRoot root {
        .transmittanceRW = m_TransimittanceTexture.GetUavIndex(),
        .multiscatteringRW = m_MultiscatteringTexture.GetUavIndex(),
        .scatteringRW = m_ScatteringTexture.GetUavIndex(),
        .singleMieScatteringRW = m_SingleMieScatteringTexture.GetUavIndex(),
    };

    cmd->SetComputeRootConstant(root);

    m_TransimittanceTexture.BarrierTransitionToUav(cmd);
    cmd->SetPipelineState(m_TransmittancePSO.Get());
    cmd->Dispatch(kTransmittanceTextureDim, {8, 8});
    m_TransimittanceTexture.BarrierTransitionToSrv(cmd, false);

    m_MultiscatteringTexture.BarrierTransitionToUav(cmd);
    cmd->SetPipelineState(m_MultiscatteringPSO.Get());
    cmd->Dispatch(kMultiscatteringTextureDim);
    m_MultiscatteringTexture.BarrierTransitionToSrv(cmd, false);

    m_ScatteringTexture.BarrierTransitionToUav(cmd, false);
    m_SingleMieScatteringTexture.BarrierTransitionToUav(cmd);
    cmd->SetPipelineState(m_ScatteringPSO.Get());
    cmd->Dispatch({SCATTERING_TEXTURE_WIDTH, SCATTERING_TEXTURE_HEIGHT, SCATTERING_TEXTURE_DEPTH}, {8, 8, 1});
    m_ScatteringTexture.BarrierTransitionToSrv(cmd, false);
    m_SingleMieScatteringTexture.BarrierTransitionToSrv(cmd);

    
    ShaderResource::NoiseRoot noiseRoot{
        .noise2DRW = m_Noise2D.GetUavIndex(),
        .noise3DRW = m_Noise3D.GetUavIndex(),
    };

    cmd->SetComputeRootConstant(noiseRoot);

    m_Noise2D.BarrierTransitionToUav(cmd);
    cmd->SetPipelineState(m_Noise2DPSO.Get());
    cmd->Dispatch(kNoise2DSize, {8, 8});
    m_Noise2D.BarrierTransitionToSrv(cmd, false);

    m_Noise3D.BarrierTransitionToUav(cmd);
    cmd->SetPipelineState(m_Noise3DPSO.Get());
    cmd->Dispatch(kNoise3DSize, {4, 4, 4});
    m_Noise3D.BarrierTransitionToSrv(cmd);
}
