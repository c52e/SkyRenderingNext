#include "Renderer.h"

#include <fstream>

#include "stb_image_write.h"

#include "Core/Gui.h"
#include "Graphics/Shader.h"
#include "Graphics/ShaderResources.h"
#include "RenderManager.h"

Renderer::Renderer()
    : m_RenderTexture(DXGI_FORMAT_R32G32B32A32_FLOAT, true, L"Render Texture"),
      m_RawRenderTexture(DXGI_FORMAT_R32G32B32A32_FLOAT, true, L"Raw Render Texture"),
      m_DepthTexture(DXGI_FORMAT_R32_FLOAT, true, L"Depth Texture"),
      m_MotionVectorTexture(DXGI_FORMAT_R32G32_FLOAT, true, L"Motion Vector Texture"),
      m_DiffuseAlbedoTexture(DXGI_FORMAT_R11G11B10_FLOAT, true, L"Diffuse Albedo Texture"),
      m_SpecularAlbedoTexture(DXGI_FORMAT_R11G11B10_FLOAT, true, L"Specular Albedo Texture"),
      m_NormalRoughnessTexture(DXGI_FORMAT_R16G16B16A16_FLOAT, true, L"Normal Roughness Texture")
{
    LoadShader();
    LoadScene("scene.json");
}

void Renderer::LoadScene(const char *path)
{
    std::ifstream fin(path);
    if (fin)
    {
        try
        {
            auto inJson = ser::json::parse(fin);
            inJson.at("scene").get_to(*this);
            inJson.at("camera").get_to(GetRenderManager().MainCamera());
            if (inJson.find("objects") != inJson.end())
            {
                std::vector<LoadResourceInfo> loadInfos;
                inJson.at("objects").get_to(loadInfos);
                m_LoadQueue.append_range(loadInfos);
            }
            GetWorld().clear();
            GetModelCache().clear();
        }
        catch (std::exception &e)
        {
            LOG_ERROR("Load scene failed: {}", e.what());
        }
    }
    else
    {
        StoreScene(path);
        LOG_WARN("\"{}\" not found. Default scene created", path);
    }
    strcpy_s(m_ScenePath, path);
}

void Renderer::StoreScene(const char* path)
{
    std::ofstream fout(path);
    if (!fout)
    {
        LOG_ERROR("Save scene failed");
        return;
    }
    ser::json outJson;
    outJson["scene"] = *this;
    outJson["camera"] = GetRenderManager().MainCamera();
    std::vector<LoadResourceInfo> loadInfos;
    for (const auto &[_, path, transform] : GetWorld().view<ResourcePath, Transform>().each())
    {
        loadInfos.emplace_back(path.value.string(), transform.m);
    }
    outJson["objects"] = loadInfos;
    fout << outJson.dump(2) << std::endl;
}

void Renderer::LoadShader()
{
    auto &device = GetGraphicsDevice();
    auto d3d12Device = device.GetD3D12Device();

    auto vs = CompileShader(SHADER_PATH("Blit.hlsl"), L"VS", L"vs_6_6");
    auto ps = CompileShader(SHADER_PATH("Blit.hlsl"), L"PS", L"ps_6_6");

    if (vs && ps)
    {
        D3D12_GRAPHICS_PIPELINE_STATE_DESC desc{
            .pRootSignature = device.GetBindlessRootSignature(),
            .VS = {reinterpret_cast<BYTE *>(vs->GetBufferPointer()), vs->GetBufferSize()},
            .PS = {reinterpret_cast<BYTE *>(ps->GetBufferPointer()), ps->GetBufferSize()},
            .BlendState{.AlphaToCoverageEnable = false,
                        .IndependentBlendEnable = false,
                        .RenderTarget{{.BlendEnable = false, .RenderTargetWriteMask = D3D12_COLOR_WRITE_ENABLE_ALL}}},
            .SampleMask = UINT_MAX,
            .RasterizerState = {.FillMode = D3D12_FILL_MODE_SOLID,
                                .CullMode = D3D12_CULL_MODE_BACK,
                                .FrontCounterClockwise = FALSE,
                                .DepthBias = D3D12_DEFAULT_DEPTH_BIAS,
                                .DepthBiasClamp = D3D12_DEFAULT_DEPTH_BIAS_CLAMP,
                                .SlopeScaledDepthBias = D3D12_DEFAULT_SLOPE_SCALED_DEPTH_BIAS,
                                .DepthClipEnable = TRUE,
                                .MultisampleEnable = FALSE,
                                .AntialiasedLineEnable = FALSE,
                                .ForcedSampleCount = 0,
                                .ConservativeRaster = D3D12_CONSERVATIVE_RASTERIZATION_MODE_OFF},
            .PrimitiveTopologyType = D3D12_PRIMITIVE_TOPOLOGY_TYPE_TRIANGLE,
            .NumRenderTargets = 1,
            .RTVFormats{kBackBufferFormat},
            .SampleDesc{.Count = 1, .Quality = 0},
        };

        d3d12Device->CreateGraphicsPipelineState(&desc, IID_PPV_ARGS(&m_BlitPSO));
    }

    CompileComputePSO(m_PathTracerPSO, SHADER_PATH("PathTracer.hlsl"), L"CS", {{L"REALTIME_PATH_TRACING", L"0"}});
    CompileComputePSO(m_RealtimePathTracerPSO, SHADER_PATH("PathTracer.hlsl"), L"CS",
                      {{L"REALTIME_PATH_TRACING", L"1"}});

    ASSERT(m_BlitPSO && m_PathTracerPSO && m_RealtimePathTracerPSO);
}

void Renderer::Render(GraphicsCommandList *cmd)
{
    auto &device = GetGraphicsDevice();
    auto targetSize = GetTargetSize();
    auto &renderManager = GetRenderManager();
    auto &cam = renderManager.MainCamera();

    cam.m_Aspect = static_cast<float>(targetSize.x) / static_cast<float>(targetSize.y);

    auto setRenderTextureScreenSizeIfNeed = [targetSize](RenderTexture &renderTexture) {
        const auto &desc = renderTexture.GetDesc();
        if (targetSize.x != desc.Width || targetSize.y != desc.Height)
        {
            renderTexture = RenderTexture(targetSize, renderTexture);
        }
    };
    setRenderTextureScreenSizeIfNeed(m_RenderTexture);
    setRenderTextureScreenSizeIfNeed(m_RawRenderTexture);
    setRenderTextureScreenSizeIfNeed(m_DepthTexture);
    setRenderTextureScreenSizeIfNeed(m_MotionVectorTexture);
    setRenderTextureScreenSizeIfNeed(m_DiffuseAlbedoTexture);
    setRenderTextureScreenSizeIfNeed(m_SpecularAlbedoTexture);
    setRenderTextureScreenSizeIfNeed(m_NormalRoughnessTexture);

    FrameParams frameParams{
        .VP = cam.ViewProjection(),
        .mainLightTheta = m_MainLightTheta,
        .mainLightPhi = m_MainLightPhi,
        .focusingDistance = cam.focusingDistance,
        .fNumber = cam.fNumber,
        .EnableDOF = m_EnableDOF,
        .atmosphereParameters = m_AtmosphereParameters,
        .targetSize = targetSize,
    };
    if (memcmp(&frameParams, &m_PreFrameParams, sizeof(frameParams)) != 0)
    {
        ResetHistory();
    }
    m_PreFrameParams = frameParams;

    while (!m_LoadQueue.empty())
    {
        auto loadInfo = std::move(m_LoadQueue.front());
        auto loadPath = std::filesystem::path(loadInfo.path);
        m_LoadQueue.pop_front();
        auto ext = loadPath.extension().string();
        auto isext = [&ext](const char *target) { return lstrcmpiA(ext.c_str(), target) == 0; };
        if (isext(".jpg") || isext(".png"))
        {
            m_Tex.CreateFromPath(cmd, loadPath);
        }
        else if (isext(".gltf") || isext(".glb"))
        {
            try
            {
                AddModel(cmd, loadPath, loadInfo.matrix);
                m_TlasDirty = true;
                ResetHistory();
            }
            catch (std::runtime_error &e)
            {
                LOG_ERROR("Loading {} failed\n\t{}", loadPath.string(), e.what());
            }
        }
    }
     
    Math::uint2 dispatchGroupSize{};
    Math::uint2 dispatchGroupOffset{};
    {
        auto totalGroupCount = (targetSize + 7u) / 8u;
        bool newFrame = m_FpsStablizer.FrameUpdate(totalGroupCount, kTargetFPS[m_TargetFPSIndex], dispatchGroupSize,
                                                   dispatchGroupOffset);
        m_AccFrameCount += newFrame;
    }

    if (m_AccFrameCount >= 100 && m_AccFrameCount % 100 == 0)
    {
        LOG_INFO("SPP: {}", m_AccFrameCount);
    }

    SetGlobalConstant(cmd);
    if (!m_AtmosphereParametersPre || *m_AtmosphereParametersPre != m_AtmosphereParameters)
    {
        m_Atmosphere.UpdateLuts(cmd);
        m_AtmosphereParametersPre = m_AtmosphereParameters;
    }

    m_FocusingRenderer.Execute(cmd, m_Tlas);

    m_RenderTexture.BarrierTransitionToUav(cmd);

    cmd->SetPipelineState(m_PathTracerPSO.Get());

    if (m_TlasDirty)
    {
        std::vector<GltfModel *> models;
        std::vector<Math::float4x4> matrices;
        models.reserve(GetWorld().view<Transform>().size());
        matrices.reserve(GetWorld().view<Transform>().size());
        for (const auto &[_, model, transform] : GetWorld().view<ModelResourceHandle, Transform>().each())
        {
            models.push_back(model.get());
            matrices.push_back(transform.m);
        }
        m_Tlas = Tlas(cmd, std::span(models), std::span(matrices));
        m_TlasDirty = false;
    }

    ShaderResource::PathTracerConstant pathTracerConstant{
        .rtasRoot = m_Tlas.GetRootConstant(),
        .accFrameCount = m_AccFrameCount,
        .debugView = m_DebugView,
        .dispatchGroupOffset = dispatchGroupOffset,
    };
    ShaderResource::PathTracerRoot pathTracerRoot{
        .cb = device.AllocTransientBuffer(pathTracerConstant).index(),
        .texRW = m_RenderTexture.GetUavIndex(),
        .depthRW = m_DepthTexture.GetUavIndex(),
        .mvRW = m_MotionVectorTexture.GetUavIndex(),
        .diffuseAlbedoRW = m_DiffuseAlbedoTexture.GetUavIndex(),
        .specularAlbedoRW = m_SpecularAlbedoTexture.GetUavIndex(),
        .normalRoughnessRW = m_NormalRoughnessTexture.GetUavIndex(),
    };
    cmd->SetComputeRootConstant(pathTracerRoot);
    cmd->Dispatch(dispatchGroupSize);

    auto jitter = cam.jitter * 0.5f * Math::float2(targetSize);

    SaveScreenShot(cmd);
    m_RenderTexture.BarrierTransitionToSrv(cmd);

    auto screenSize = GetRenderManager().GetScreenSize();
    D3D12_VIEWPORT viewport{
        .TopLeftX = 0.0f,
        .TopLeftY = 0.0f,
        .Width = static_cast<float>(screenSize[0]),
        .Height = static_cast<float>(screenSize[1]),
        .MinDepth = 0.0f,
        .MaxDepth = 1.0f,
    };
    cmd->RSSetViewports(1, &viewport);
    D3D12_RECT scissor{0, 0, static_cast<LONG>(screenSize[0]), static_cast<LONG>(screenSize[1])};
    cmd->RSSetScissorRects(1, &scissor);
    cmd->SetPipelineState(m_BlitPSO.Get());

    ShaderResource::BlitRoot blitRoot{
        .tex = m_RenderTexture.GetSrvIndex(),
        .exposure = m_Exposure,
    };
    cmd->SetGraphicsRootConstant(blitRoot);
    cmd->IASetPrimitiveTopology(D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
    cmd->DrawInstanced(3, 1, 0, 0);
}

void Renderer::RenderGui()
{
    ImGui::Begin("Render Settings");

    auto &io = ImGui::GetIO();
    ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / io.Framerate, io.Framerate);

    auto &cam = GetRenderManager().MainCamera();

    {
        static std::vector<std::string> scenePaths;
        static std::vector<const char *> scenePathsCstr;
        static int index = 0;
        if (ImGui::Button("Load scene.."))
        {
            ImGui::OpenPopup("Load scene");
            namespace fs = std::filesystem;
            scenePaths.clear();
            scenePathsCstr.clear();
            for (const auto folder : { fs::current_path(), fs::current_path() / "../data/scene" })
            {
                if (!fs::exists(folder))
                {
                    fs::create_directory(folder);
                }
                for (auto const &entry : fs::recursive_directory_iterator(folder))
                {
                    if (fs::is_regular_file(entry) && entry.path().extension() == ".json")
                    {
                        scenePaths.emplace_back(fs::relative(entry.path(), fs::current_path()).string());
                        if (scenePaths.back() == m_ScenePath)
                        {
                            index = static_cast<int>(scenePaths.size()) - 1;
                        }
                    }
                }
            }
            for (const auto& path : scenePaths)
            {
                scenePathsCstr.emplace_back(path.c_str());
            }
        }
        if (ImGui::BeginPopupModal("Load scene", nullptr, ImGuiWindowFlags_AlwaysAutoResize))
        {
            if (!scenePathsCstr.empty())
            {
                if (index >= static_cast<int>(scenePathsCstr.size()))
                {
                    index = static_cast<int>(scenePathsCstr.size()) - 1;
                }
                ImGui::ListBox("##Scenes", &index, scenePathsCstr.data(),
                               static_cast<int>(scenePathsCstr.size()), 24);
                if (ImGui::Button("OK", ImVec2(100, 0)))
                {
                    LoadScene(scenePathsCstr[index]);
                    ImGui::CloseCurrentPopup();
                }
            }
            ImGui::SameLine();
            if (ImGui::Button("Cancel", ImVec2(100, 0)))
            {
                ImGui::CloseCurrentPopup();
            }
            ImGui::EndPopup();
        }
    }
    ImGui::SameLine();

    if (ImGui::Button("Save scene.."))
    {
        ImGui::OpenPopup("Save scene");
    }
    if (ImGui::BeginPopupModal("Save scene", nullptr, ImGuiWindowFlags_AlwaysAutoResize))
    {
        ImGui::InputText("path", m_ScenePath, std::size(m_ScenePath));
        if (ImGui::Button("OK", ImVec2(120, 0)))
        {
            StoreScene(m_ScenePath);
            ImGui::CloseCurrentPopup();
        }
        ImGui::SameLine();
        if (ImGui::Button("Cancel", ImVec2(120, 0)))
        {
            ImGui::CloseCurrentPopup();
        }
        ImGui::EndPopup();
    }
    ImGui::SameLine();

    if (ImGui::Button("ScreenShot"))
    {
        time_t rawtime{};
        std::time(&rawtime);
        std::tm timeinfo{};
        localtime_s(&timeinfo, &rawtime);
        char buffer[84];
        std::strftime(buffer, std::size(buffer), "%Y-%m-%d %H-%M-%S", &timeinfo);
        m_WaitingScreenShotName = buffer;
    }

    ImGui::RatioTable("Resolution", kResolutionNames, &m_Resolution);
    ImGui::RatioTable("Target FPS", kTargetFPSNames, &m_TargetFPSIndex);

    ImGui::SliderFloat("Focal Length", &cam.focalLength, 4.0f, 200.0f);
    ImGui::Checkbox("DOF", &m_EnableDOF);
    if (m_EnableDOF)
    {
        ImGui::SliderFloat("Focusing Distance", &cam.focusingDistance, cam.focalLength * 0.001f, 30.0f);
        ImGui::SliderFloat("F Number", &cam.fNumber, 0.1f, 22.0f);
    }
    ImGui::SliderFloat("Exposure", &m_Exposure, 0.0f, 10.0f);
    ImGui::SliderFloat("Main Light Dir Theta", &m_MainLightTheta, 0.0f, 1.0f * Math::pi<float>());
    ImGui::SliderFloat("Main Light Dir Phi", &m_MainLightPhi, 0.0f, 2.0f * Math::pi<float>());
    ImGui::SliderInt("DebugView", &m_DebugView, 0, 6);

    ImGui::SliderFloat("Camera Speed", &cam.movingSpeed, 0.0f, 3000.0f, "%.1f", ImGuiSliderFlags_Logarithmic);

    ImGui::Begin("Atmosphere");
    {
        ImGui::ColorEdit3("Solar Illuminance", Math::value_ptr(m_AtmosphereParameters.solarIlluminance));
        ImGui::SliderFloat("Sun Angular Radius", &m_AtmosphereParameters.sunAngularRadius, 0, 30);
        ImGui::Separator();

        ImGui::ColorEdit3("Ground Albedo", Math::value_ptr(m_AtmosphereParameters.groundAlbedo));
        ImGui::SliderFloat("Bottom Radius", &m_AtmosphereParameters.bottomRadius, 0, 6360e3f, "%.0f",
                           ImGuiSliderFlags_Logarithmic);
        ImGui::SliderFloat("Thickness", &m_AtmosphereParameters.thickness, 0, 1000e3f, "%.0f",
                           ImGuiSliderFlags_Logarithmic);
        ImGui::Separator();

        ImGui::SliderFloat("Rayleigh Exponential Distribution", &m_AtmosphereParameters.rayleighExponentialDistribution,
                           0, 1000e3f, "%.0f", ImGuiSliderFlags_Logarithmic);
        ImGui::SliderFloat("Rayleigh Scattering Scale", &m_AtmosphereParameters.rayleighScatteringScale, 0, 1.0e-3f,
                           "%.7f", ImGuiSliderFlags_Logarithmic);
        ImGui::ColorEdit3("Rayleigh Scattering", Math::value_ptr(m_AtmosphereParameters.rayleighScattering));
        ImGui::Separator();

        ImGui::SliderFloat("Mie Exponential Distribution", &m_AtmosphereParameters.mieExponentialDistribution, 0,
                           1000e3f, "%.0f", ImGuiSliderFlags_Logarithmic);
        ImGui::SliderFloat("Mie Phase G", &m_AtmosphereParameters.miePhaseG, -1.0, 1.0);
        ImGui::SliderFloat("Mie Scattering Scale", &m_AtmosphereParameters.mieScatteringScale, 0, 1.0e-3f, "%.7f",
                ImGuiSliderFlags_Logarithmic);
        ImGui::ColorEdit3("Mie Scattering", Math::value_ptr(m_AtmosphereParameters.mieScattering));
        ImGui::SliderFloat("Mie Absorption Scale", &m_AtmosphereParameters.mieAbsorptionScale, 0, 1.0e-3f, "%.7f",
                ImGuiSliderFlags_Logarithmic);
        ImGui::ColorEdit3("Mie Absorption", Math::value_ptr(m_AtmosphereParameters.mieAbsorption));
        ImGui::Separator();

        ImGui::SliderFloat("Ozone Center Altitude", &m_AtmosphereParameters.ozoneCenterAltitude, 0, 1000e3f, "%.0f",
                           ImGuiSliderFlags_Logarithmic);
        ImGui::SliderFloat("Ozone Width", &m_AtmosphereParameters.ozoneWidth, 0, 1000e3f, "%.0f",
                           ImGuiSliderFlags_Logarithmic);
        ImGui::SliderFloat("Ozone Absorption Scale", &m_AtmosphereParameters.ozoneAbsorptionScale, 0, 1.0e-3f, "%.7f",
                           ImGuiSliderFlags_Logarithmic);
        ImGui::ColorEdit3("Ozone Absorption",Math::value_ptr( m_AtmosphereParameters.ozoneAbsorption));
    }
    ImGui::End();

    ImGui::Begin("Cloud");
    {
        ImGui::SliderFloat("Cloud Bottom", &m_AtmosphereParameters.cloud.bottomAltitude, 0.0f, 8000.0f);
        ImGui::SliderFloat("Cloud Thickness", &m_AtmosphereParameters.cloud.thickness, 0.0f, 8000.0f);
        ImGui::SliderFloat("Cloud Extinction", &m_AtmosphereParameters.cloud.extinction, 0.0f, 1.0f, "%.3f",
                           ImGuiSliderFlags_Logarithmic);

        ImGui::Separator();

        ImGui::SliderFloat("baseWidth2D", &m_AtmosphereParameters.cloud.baseWidth2D, 5e3f, 100e3f);
        ImGui::SliderFloat("baseWidth3D", &m_AtmosphereParameters.cloud.baseWidth3D, 0.5e3f, 10e3f);
        ImGui::SliderFloat("coverage", &m_AtmosphereParameters.cloud.coverage, 0.0f, 1.0f);
        ImGui::SliderFloat("cloudBaseHeight", &m_AtmosphereParameters.cloud.cloudBaseHeight, 0.0f, 1.0f);
        ImGui::SliderFloat("lowFrequencyNoise", &m_AtmosphereParameters.cloud.lowFrequencyNoise, 0.0f, 1.0f);
        ImGui::SliderFloat("highFrequencyNoise", &m_AtmosphereParameters.cloud.highFrequencyNoise, 0.0f, 1.0f);
        ImGui::SliderFloat("taperPosition", &m_AtmosphereParameters.cloud.taperPosition, 0.0f, 1.0f);
        ImGui::SliderFloat("taperFloor", &m_AtmosphereParameters.cloud.taperFloor, 0.0f, 1.0f);
        ImGui::SliderInt("max bounce", &m_AtmosphereParameters.cloud.maxBounce, 1, 128);
    }
    ImGui::End();

    ImGui::End();

    if (ImGui::IsKeyDown(ImGuiKey_LeftCtrl) && ImGui::IsKeyPressed(ImGuiKey_R, false))
    {
        LoadShader();
        ResetHistory();
    }

    m_FocusingRenderer.RenderGui();
    m_TlasDirty |= EditWorld();
}

void Renderer::SetGlobalConstant(GraphicsCommandList *cmd)
{
    auto &renderManager = GetRenderManager();
    auto &cam = renderManager.MainCamera();

    auto jitterValue = [this](uint32_t frameIndex) {
        return Math::Halton23(frameIndex, 8) * 2.0f / Math::float2(GetTargetSize());
    };

    Math::float2 jitterMotion{};
    {
        cam.jitter = {};
    }

    auto VP = cam.ViewProjection();
    auto invVP = Math::inverse(VP);

    float cosTheta = cos(m_MainLightTheta);
    float sinTheta = sin(m_MainLightTheta);
    float cosPhi = cos(m_MainLightPhi);
    float sinPhi = sin(m_MainLightPhi);

    ShaderResource::GlobalConstant globalConstant{
        .VP = VP,
        .invVP = invVP,
        .PreVP = cam.PreViewProjection(),
        .camPos = cam.m_Position,
        .lensDistance = cam.LensDistance(),
        .mainLightDir = Math::float3(sinPhi * sinTheta, cosTheta, cosPhi * sinTheta),
        .apertureRadius = cam.ApertureRadius(),
        .camFront = cam.GetFront(),
        .focusingDistance = cam.focusingDistance,
        .targetSize = GetTargetSize(),
        .screenSize = renderManager.GetScreenSize(),
        .jitterMotion = jitterMotion,
        .EnableDOF = m_EnableDOF,
        .atmosphere = m_Atmosphere.UpdateParameters(m_AtmosphereParameters),
    };

    UINT64 addr = GetGraphicsDevice().AllocTransientBufferAddress(globalConstant);
    cmd->SetGraphicsRootConstantBufferView(1, addr);
    cmd->SetComputeRootConstantBufferView(1, addr);
}

void Renderer::SaveScreenShot(GraphicsCommandList *cmd)
{
    while (!m_WaitingSaveImages.empty() && m_WaitingSaveImages.front().buffer.Available())
    {
        const auto &[name, buffer, footprint] = m_WaitingSaveImages.front();

        namespace fs = std::filesystem;
        fs::path folder = "gallery/";
        if (!fs::exists(folder))
        {
            fs::create_directory(folder);
        }
        auto width = footprint.Footprint.Width;
        auto height = footprint.Footprint.Height;
        ASSERT(footprint.Footprint.Format == DXGI_FORMAT_R32G32B32A32_FLOAT);
        std::vector<float> tightBuffer;
        tightBuffer.reserve(static_cast<UINT64>(width) * height * 3);
        auto rowData = reinterpret_cast<std::byte *>(buffer.GetPtr());
        for (UINT y = 0; y < height; ++y)
        {
            auto data = reinterpret_cast<float *>(rowData);
            for (UINT x = 0; x < width; ++x)
            {
                tightBuffer.push_back(*data++);
                tightBuffer.push_back(*data++);
                tightBuffer.push_back(*data++);
                data++; // ignore alpha
            }
            rowData += footprint.Footprint.RowPitch;
        }
        stbi_flip_vertically_on_write(true);
        stbi_write_hdr((folder / (name + ".hdr")).string().c_str(), width, height, 3, tightBuffer.data());
        m_WaitingSaveImages.pop_front();
    }
    if (!m_WaitingScreenShotName.empty())
    {
        D3D12_PLACED_SUBRESOURCE_FOOTPRINT footprint{};
        auto buffer = ReadbackRenderTexture(cmd, m_RenderTexture, footprint);
        m_WaitingSaveImages.push_back({std::move(m_WaitingScreenShotName), std::move(buffer), footprint});
    }
}

Math::uint2 Renderer::GetTargetSize() const
{
    Math::uint2 sizes[]{GetRenderManager().GetScreenSize(), {256, 256}, {1920, 1080}, {3840, 2160}, {7040, 4688}};
    return sizes[m_Resolution];
}

FocusingRenderer::FocusingRenderer()
    : m_FocusingInfoBuffer(1, sizeof(float), true), m_FocusingInfoReadbackBuffer(sizeof(float))
{
    CompileComputePSO(m_FocusingPSO, SHADER_PATH("Focusing.hlsl"), L"CS");
    ASSERT(m_FocusingPSO);
}

void FocusingRenderer::Execute(GraphicsCommandList *cmd, const Tlas &tlas)
{
    if (m_FocusingState >= 0)
    {
        if (m_FocusingState == 0)
        {
            ShaderResource::FocusingRoot focusingRoot{
                .rtasRoot = tlas.GetRootConstant(),
                .focusingPos = m_FocusingPos,
                .focusingInfoBufferRW = m_FocusingInfoBuffer.GetUavIndex(),
            };

            m_FocusingInfoBuffer.BarrierTransitionToUav(cmd);
            cmd->SetPipelineState(m_FocusingPSO.Get());
            cmd->SetComputeRootConstant(focusingRoot);
            cmd->Dispatch(1, 1, 1);

            m_FocusingInfoReadbackBuffer.ReadbackFrom(cmd, m_FocusingInfoBuffer);
            m_FocusingState = 1;
        }
        else if (m_FocusingInfoReadbackBuffer.Available())
        {
            GetRenderManager().MainCamera().focusingDistance =
                *reinterpret_cast<float *>(m_FocusingInfoReadbackBuffer.GetPtr());
            m_FocusingState = -1;
        }
    }
}

void FocusingRenderer::RenderGui()
{
    if (ImGui::IsMouseClicked(ImGuiMouseButton_Left) && ImGui::IsKeyDown(ImGuiKey_LeftCtrl))
    {
        POINT p{};
        if (GetCursorPos(&p))
        {
            ScreenToClient(GetRenderManager().GetHwnd(), &p);
            m_FocusingState = 0;
            m_FocusingPos = {p.x, p.y};
            LOG_INFO("Focusing ({} {})", m_FocusingPos.x, m_FocusingPos.y);
        }
    }
}

void FpsStablizer::Reset()
{
    m_TileIndex = {0, 0};
    m_LastFrameTime = {};
}

bool FpsStablizer::FrameUpdate(Math::uint2 totalGroupCount, float targetFPS, Math::uint2 &outGroupPerTile,
                               Math::uint2 &outGroupOffset)
{
    bool newFrame = (m_TileIndex.x == 0 && m_TileIndex.y == 0);
    if (newFrame)
    {
        auto now = std::chrono::high_resolution_clock::now();
        if (m_LastFrameTime.has_value())
        {
            auto duration = now - *m_LastFrameTime;
            float seconds = static_cast<float>(duration.count()) * 1e-9f;
            auto oldTileCount = m_TileCount;
            float secondsPerTile = seconds / (m_TileCount.x * m_TileCount.y);
            while (secondsPerTile > (1.0f + kMaxTimeDiff) / targetFPS)
            {
                auto groupCountPerTile = (totalGroupCount + m_TileCount - 1u) / m_TileCount;
                int index = groupCountPerTile.x <= groupCountPerTile.y;
                if (m_TileCount[index] + 1 <= totalGroupCount[index])
                {
                    m_TileCount[index] += 1;
                }
                else if (m_TileCount[1 - index] + 1 <= totalGroupCount[1 - index])
                {
                    m_TileCount[1 - index] += 1;
                }
                else
                {
                    break;
                }
                secondsPerTile = seconds / (m_TileCount.x * m_TileCount.y);
            }
            while (secondsPerTile < 1.0f / (1.0f + kMaxTimeDiff) / targetFPS)
            {
                auto groupCountPerTile = (totalGroupCount + m_TileCount - 1u) / m_TileCount;
                int index = groupCountPerTile.x <= groupCountPerTile.y;
                if (m_TileCount[index] > 1)
                {
                    m_TileCount[index] -= 1;
                }
                else if (m_TileCount[1 - index] > 1)
                {
                    m_TileCount[1 - index] -= 1;
                }
                else
                {
                    break;
                }
                secondsPerTile = seconds / (m_TileCount.x * m_TileCount.y);
            }
            if (m_TileCount.x != oldTileCount.x || m_TileCount.y != oldTileCount.y)
            {
                LOG_INFO("Change tile count to ({}, {})", m_TileCount.x, m_TileCount.y);
            }
        }
        m_LastFrameTime = now;
    }

    outGroupPerTile = (totalGroupCount + m_TileCount - 1u) / m_TileCount;
    outGroupOffset = outGroupPerTile * ((m_TileIndex + m_TileCount / 2u) % m_TileCount); // Render center first

    m_TileIndex.x = (m_TileIndex.x + 1) % m_TileCount.x;
    m_TileIndex.y = (m_TileIndex.y + (m_TileIndex.x == 0)) % m_TileCount.y;

    return newFrame;
}
