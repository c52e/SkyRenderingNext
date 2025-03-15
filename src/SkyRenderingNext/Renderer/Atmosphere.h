#pragma once

#include "Core/Common.h"
#include "Core/Serialization.h"
#include "Graphics/GraphicsDevice.h"
#include "Graphics/RenderTexture.h"
#include "Graphics/Texture2D.h"
#include "Graphics/ShaderResources.h"

struct CloudParameters
{
    float bottomAltitude = 1300.0f;
    float thickness = 1500.0f;
    float extinction = 0.04f;

    float baseWidth2D = 35e3f;
    float baseWidth3D = 2e3f;
    float coverage = 0.33f;
    float cloudBaseHeight = 0.1f;
    float lowFrequencyNoise = 0.3f;
    float highFrequencyNoise = 0.08f;
    float taperPosition = 0.1f;
    float taperFloor = 0.05f;
    int maxBounce = 128;

    auto operator<=>(const CloudParameters &) const = default;
};

struct AtmosphereParameters
{
    Math::float3 solarIlluminance{1.f, 1.f, 1.f};
    float sunAngularRadius = 0.5334f * 0.5f;

    float bottomRadius = 6360e3f;
    float thickness = 6420e3f - 6360e3f;
    Math::float3 groundAlbedo{0.4f, 0.4f, 0.4f};

    float rayleighExponentialDistribution = 8e3f;
    float rayleighScatteringScale = 0.0331e-3f;
    Math::float3 rayleighScattering{0.175287f, 0.409607f, 1.0f};

    float mieExponentialDistribution = 1.2e3f;
    float miePhaseG = 0.8f;
    float mieScatteringScale = 0.003996e-3f;
    Math::float3 mieScattering{1.0f, 1.0f, 1.0f};
    float mieAbsorptionScale = 0.000444e-3f;
    Math::float3 mieAbsorption{1.0f, 1.0f, 1.0f};

    float ozoneCenterAltitude = 25e3f;
    float ozoneWidth = 15e3f;
    float ozoneAbsorptionScale = 0.001881e-3f;
    Math::float3 ozoneAbsorption{0.345561f, 1.0f, 0.045189f};

    float transmittanceSteps = 400.0f;
    float multiscatteringSteps = 100.0f;

    CloudParameters cloud{};

    auto operator<=>(const AtmosphereParameters &) const = default;
};

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(CloudParameters
    , bottomAltitude
    , bottomAltitude
    , thickness
    , extinction
    , baseWidth2D
    , baseWidth3D
    , coverage
    , cloudBaseHeight
    , lowFrequencyNoise
    , highFrequencyNoise
    , taperPosition
    , taperFloor
    , maxBounce)
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(AtmosphereParameters
    , solarIlluminance
    , sunAngularRadius
    , bottomRadius
    , thickness
    , groundAlbedo
    , rayleighExponentialDistribution
    , rayleighScatteringScale
    , rayleighScattering
    , mieExponentialDistribution
    , miePhaseG
    , mieScatteringScale
    , mieScattering
    , mieAbsorptionScale
    , mieAbsorption
    , ozoneCenterAltitude
    , ozoneWidth
    , ozoneAbsorptionScale
    , ozoneAbsorption
    , transmittanceSteps
    , multiscatteringSteps
    , cloud)

class Atmosphere
{
  public:
    Atmosphere();

    ShaderResource::AtmosphereConstant UpdateParameters(const AtmosphereParameters &parameters);

    void UpdateLuts(GraphicsCommandList *cmd);

  private:
    ResPtr<ID3D12PipelineState> m_TransmittancePSO{};
    ResPtr<ID3D12PipelineState> m_MultiscatteringPSO{};
    ResPtr<ID3D12PipelineState> m_ScatteringPSO{};

    RenderTexture m_TransimittanceTexture{};
    RenderTexture m_MultiscatteringTexture{};
    RenderTexture m_ScatteringTexture{};
    RenderTexture m_SingleMieScatteringTexture{};

    ResPtr<ID3D12PipelineState> m_Noise2DPSO{};
    ResPtr<ID3D12PipelineState> m_Noise3DPSO{};

    RenderTexture m_Noise2D{};
    RenderTexture m_Noise3D{};
};