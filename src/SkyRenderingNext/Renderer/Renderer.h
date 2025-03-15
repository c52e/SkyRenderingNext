#pragma once

#include <deque>
#include <chrono>

#include "Core/Common.h"
#include "Graphics/Tlas.h"
#include "Graphics/GraphicsDevice.h"
#include "Graphics/RenderBuffer.h"
#include "Graphics/RenderTexture.h"
#include "Graphics/Texture2D.h"
#include "Atmosphere.h"
#include "World.h"

class FocusingRenderer
{
  public:
    FocusingRenderer();

    void Execute(GraphicsCommandList *cmd, const Tlas &tlas);

    void RenderGui();

  private:
    ResPtr<ID3D12PipelineState> m_FocusingPSO{};

    StructuredBuffer m_FocusingInfoBuffer{};
    ReadbackBuffer m_FocusingInfoReadbackBuffer{};
    Math::int2 m_FocusingPos{};
    int m_FocusingState = -1;
};

class FpsStablizer
{
  public:
    void Reset();

    bool FrameUpdate(Math::uint2 totalGroupCount, float targetFPS, Math::uint2 &outGroupPerTile,
                     Math::uint2 &outGroupOffset);

  private:
    static constexpr float kMaxTimeDiff = 0.2f;

    Math::uint2 m_TileIndex{0, 0};
    Math::uint2 m_TileCount{1, 1};
    std::optional<std::chrono::high_resolution_clock::time_point> m_LastFrameTime{};
};

struct LoadResourceInfo
{
    std::string path;
    Math::float4x4 matrix;
};

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(LoadResourceInfo, path, matrix)

class Renderer
{
  public:
    Renderer();

    void LoadScene(const char *path);

    void Render(GraphicsCommandList *cmd);

    void RenderGui();

    void LoadResource(std::wstring_view path)
    {
        LoadResource(std::filesystem::path(path).string(), Math::float4x4(1.0f));
    }

    void LoadResource(const std::string& path, const Math::float4x4 &matrix)
    {
        m_LoadQueue.push_back({path, matrix});
    }

  private:
    void StoreScene(const char *path);

    void LoadShader();

    void SetGlobalConstant(GraphicsCommandList *cmd);

    void ResetHistory()
    {
        m_AccFrameCount = 0;
        m_FpsStablizer.Reset();
    }

    void SaveScreenShot(GraphicsCommandList *cmd);

    ResPtr<ID3D12PipelineState> m_BlitPSO{};
    ResPtr<ID3D12PipelineState> m_PathTracerPSO{};
    ResPtr<ID3D12PipelineState> m_RealtimePathTracerPSO{};

    RenderTexture m_RenderTexture{};
    RenderTexture m_RawRenderTexture{};
    RenderTexture m_DepthTexture{};
    RenderTexture m_MotionVectorTexture{};
    RenderTexture m_DiffuseAlbedoTexture{};
    RenderTexture m_SpecularAlbedoTexture{};
    RenderTexture m_NormalRoughnessTexture{};

    Texture2D m_Tex;
    bool m_TlasDirty = false;
    Tlas m_Tlas;
    std::deque<LoadResourceInfo> m_LoadQueue;

    struct FrameParams
    {
        Math::float4x4 VP{};
        float mainLightTheta{};
        float mainLightPhi{};
        float focusingDistance{};
        float fNumber{};
        bool EnableDOF{};
        AtmosphereParameters atmosphereParameters{};
        Math::uint2 targetSize{};
    } m_PreFrameParams;
    uint32_t m_AccFrameCount = 0;
    FpsStablizer m_FpsStablizer;
    int m_TargetFPSIndex = 2;
    inline static const char *kTargetFPSNames[]{"1 (Min FPS)", "10", "30", "60", "120"};
    inline static const float kTargetFPS[]{1.0f, 10.0f, 30.0f, 60.0f, 120.0f};

    bool m_EnableDOF = false;
    int m_DebugView = 0;
    float m_MainLightTheta{};
    float m_MainLightPhi{};
    float m_Exposure = 1.0f;

    FocusingRenderer m_FocusingRenderer;

    AtmosphereParameters m_AtmosphereParameters;
    std::optional<AtmosphereParameters> m_AtmosphereParametersPre;
    Atmosphere m_Atmosphere;

    char m_ScenePath[256];
    std::string m_WaitingScreenShotName;

    struct WaitingSaveImage
    {
        std::string filename;
        ReadbackBuffer buffer;
        D3D12_PLACED_SUBRESOURCE_FOOTPRINT footprint;
    };
    std::deque<WaitingSaveImage> m_WaitingSaveImages;

    Math::uint2 GetTargetSize() const;

    int m_Resolution = 2;
    inline static const char *kResolutionNames[]{"Screen Size", "256x256", "1920x1080", "3840x2160", "7040x4688"};

    NLOHMANN_DEFINE_TYPE_INTRUSIVE(Renderer
        , m_DebugView
        , m_MainLightTheta
        , m_MainLightPhi
        , m_Exposure
        , m_EnableDOF
        , m_AtmosphereParameters)
};
