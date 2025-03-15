#include "RenderManager.h"

#include "imgui.h"
#include "imgui_impl_dx12.h"
#include "imgui_impl_win32.h"
#include <DirectXColors.h>

#include "Core/Log.h"
#include "Graphics/GraphicsDevice.h"
#include "World.h"

RenderManager &GetRenderManager()
{
    static RenderManager renderManager;
    return renderManager;
}

LRESULT WINAPI WndProc(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam)
{
    extern IMGUI_IMPL_API LRESULT ImGui_ImplWin32_WndProcHandler(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam);
    if (ImGui_ImplWin32_WndProcHandler(hWnd, msg, wParam, lParam))
        return true;

    switch (msg)
    {
    case WM_CREATE:
        // Allow the window to accept file drops
        DragAcceptFiles(hWnd, TRUE);
        return 0;
    case WM_SIZE:
        if (wParam != SIZE_MINIMIZED)
        {
            GetRenderManager().WMResize((UINT)LOWORD(lParam), (UINT)HIWORD(lParam));
        }
        return 0;
    case WM_PAINT:
        GetRenderManager().LoopBody();
        return 0;
    case WM_DROPFILES: {
        HDROP hDrop = (HDROP)wParam;
        TCHAR szFileName[MAX_PATH];
        UINT fileCount = DragQueryFile(hDrop, 0xFFFFFFFF, NULL, 0);
        DragQueryFile(hDrop, 0, szFileName, MAX_PATH);
        GetRenderManager().DropEvent(szFileName);
        // Release memory
        DragFinish(hDrop);
        return 0;
    }
    case WM_DESTROY:
        ::PostQuitMessage(0);
        return 0;
    }
    return ::DefWindowProcW(hWnd, msg, wParam, lParam);
}

RenderManager::RenderManager()
{

    // Init window
    {
        WNDCLASS wc = {
            .style = CS_CLASSDC,
            .lpfnWndProc = WndProc,
            .cbClsExtra = 0,
            .cbWndExtra = 0,
            .hInstance = GetModuleHandle(nullptr),
            .hIcon = nullptr,
            .hCursor = nullptr,
            .hbrBackground = nullptr,
            .lpszMenuName = nullptr,
            .lpszClassName = L"SkyRenderingNext",
        };
        RegisterClassW(&wc);
        m_Hwnd = CreateWindowW(wc.lpszClassName, L"SkyRenderingNext", WS_OVERLAPPEDWINDOW, 100, 100, m_ScreenSize[0],
                               m_ScreenSize[1], nullptr, nullptr, wc.hInstance, nullptr);
    }

    auto device = GetGraphicsDevice().GetD3D12Device();
    auto cmdQueue = GetGraphicsDevice().GetCmdQueue();
    {
        DXGI_SWAP_CHAIN_DESC1 sd{
            .Width = m_ScreenSize[0],
            .Height = m_ScreenSize[1],
            .Format = kBackBufferFormat,
            .Stereo = false,
            .SampleDesc{
                .Count = 1,
                .Quality = 0,
            },
            .BufferUsage = DXGI_USAGE_RENDER_TARGET_OUTPUT,
            .BufferCount = kNumBackBuffer,
            .Scaling = DXGI_SCALING_NONE,
            .SwapEffect = DXGI_SWAP_EFFECT_FLIP_DISCARD,
            .AlphaMode = DXGI_ALPHA_MODE_UNSPECIFIED,
            .Flags = DXGI_SWAP_CHAIN_FLAG_FRAME_LATENCY_WAITABLE_OBJECT,
        };

        IDXGIFactory4 *dxgiFactory = nullptr;
        IDXGISwapChain1 *swapChain1 = nullptr;
        CreateDXGIFactory1(IID_PPV_ARGS(&dxgiFactory));
        dxgiFactory->CreateSwapChainForHwnd(cmdQueue, m_Hwnd, &sd, nullptr, nullptr, &swapChain1);
        swapChain1->QueryInterface(IID_PPV_ARGS(&m_Swapchain));
        swapChain1->Release();
        dxgiFactory->Release();
        m_Swapchain->SetMaximumFrameLatency(kNumBackBuffer);
        m_SwapchainWaitableObject = m_Swapchain->GetFrameLatencyWaitableObject();
    }

    {
        D3D12_DESCRIPTOR_HEAP_DESC desc{
            .Type = D3D12_DESCRIPTOR_HEAP_TYPE_RTV,
            .NumDescriptors = kNumBackBuffer,
            .Flags = D3D12_DESCRIPTOR_HEAP_FLAG_NONE,
        };
        device->CreateDescriptorHeap(&desc, IID_PPV_ARGS(&m_RtvDescHeap));
        CreateSwapchainBuffer();
    }

    {

        {
            D3D12_DESCRIPTOR_HEAP_DESC desc{.Type = D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV,
                                            .NumDescriptors = 1,
                                            .Flags = D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE};
            device->CreateDescriptorHeap(&desc, IID_PPV_ARGS(&m_ImGuiSrvDescHeap));
        }

        // Setup Dear ImGui context
        IMGUI_CHECKVERSION();
        ImGui::CreateContext();
        ImGuiIO &io = ImGui::GetIO();
        (void)io;
        io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard; // Enable Keyboard Controls
        io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;  // Enable Gamepad Controls
        io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;     // Enable Docking
        io.ConfigFlags |= ImGuiConfigFlags_ViewportsEnable;   // Enable Multi-Viewport / Platform Windows
        // io.ConfigViewportsNoAutoMerge = true;
        // io.ConfigViewportsNoTaskBarIcon = true;

        // Setup Dear ImGui style
        ImGui::StyleColorsDark();
        // ImGui::StyleColorsLight();

        // When viewports are enabled we tweak WindowRounding/WindowBg so platform windows can look identical to regular
        // ones.
        ImGuiStyle &style = ImGui::GetStyle();
        if (io.ConfigFlags & ImGuiConfigFlags_ViewportsEnable)
        {
            style.WindowRounding = 0.0f;
            style.Colors[ImGuiCol_WindowBg].w = 1.0f;
        }

        // Setup Platform/Renderer backends
        ImGui_ImplWin32_Init(m_Hwnd);
        ImGui_ImplDX12_Init(device, kNumFrameInFlights, kBackBufferFormat, m_ImGuiSrvDescHeap.Get(),
                            m_ImGuiSrvDescHeap->GetCPUDescriptorHandleForHeapStart(),
                            m_ImGuiSrvDescHeap->GetGPUDescriptorHandleForHeapStart());
    }
}

RenderManager::~RenderManager()
{
    if (m_SwapchainWaitableObject)
    {
        CloseHandle(m_SwapchainWaitableObject);
        m_SwapchainWaitableObject = nullptr;
    }

    ImGui_ImplDX12_Shutdown();
    ImGui_ImplWin32_Shutdown();
    ImGui::DestroyContext();
}

Math::uint4 RenderManager::GetWindowAbsoluteRect() const
{
    RECT clientRect{};
    GetClientRect(m_Hwnd, &clientRect);
    POINT topLeft = {clientRect.left, clientRect.top};
    POINT bottomRight = {clientRect.right,
                         clientRect.bottom}; 
    ClientToScreen(m_Hwnd, &topLeft);
    ClientToScreen(m_Hwnd, &bottomRight);
    return {topLeft.x, topLeft.y, bottomRight.x, bottomRight.y};
}

void RenderManager::MainLoop(int argc, const char *argv[])
{
    float aspect = static_cast<float>(m_ScreenSize[0]) / m_ScreenSize[1];
    m_MainCamera = Camera(Math::float3(0.0f, 0.0f, -1.0f), 90.0f, 0.0f, aspect);
    m_Renderer = std::make_unique<Renderer>();
    if (argc > 1)
    {
        m_Renderer->LoadScene(argv[1]);
    }

    ShowWindow(m_Hwnd, SW_SHOWDEFAULT);
    UpdateWindow(m_Hwnd);

    MSG msg = {};
    while (msg.message != WM_QUIT)
    {
        if (PeekMessage(&msg, nullptr, 0, 0, PM_REMOVE))
        {
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }
    }
    GetGraphicsDevice().WaitForLastSubmittedFrame(GetAndResetSwapchainWaitableObject());
}

void RenderManager::LoopBody()
{
    m_MainCamera.DefaultFrameControl();

    // Start the Dear ImGui frame
    ImGui_ImplDX12_NewFrame();
    ImGui_ImplWin32_NewFrame();
    ImGui::NewFrame();
    // ImGui::DockSpaceOverViewport();

    m_Renderer->RenderGui();

    ImGui::Render();

    GetGraphicsDevice().FrameWaitAndBegin(GetAndResetSwapchainWaitableObject());

    auto backBufferIndex = m_Swapchain->GetCurrentBackBufferIndex();

    D3D12_RESOURCE_BARRIER barrier = {.Type = D3D12_RESOURCE_BARRIER_TYPE_TRANSITION,
                                      .Flags = D3D12_RESOURCE_BARRIER_FLAG_NONE,
                                      .Transition = {
                                          .pResource = m_SwapchainBuffer[backBufferIndex].Get(),
                                          .Subresource = D3D12_RESOURCE_BARRIER_ALL_SUBRESOURCES,
                                          .StateBefore = D3D12_RESOURCE_STATE_PRESENT,
                                          .StateAfter = D3D12_RESOURCE_STATE_RENDER_TARGET,
                                      }};

    auto cmd = GetGraphicsDevice().CurCmdList();

    cmd->ResourceBarrier(barrier);
    cmd->OMSetRenderTargets(1, &m_SwapchainBufferHandle[backBufferIndex], true, nullptr);
    cmd->ClearRenderTargetView(m_SwapchainBufferHandle[backBufferIndex], DirectX::Colors::Black, 0, nullptr);

    m_Renderer->Render(cmd);

    {
        ID3D12DescriptorHeap *descHeaps[] = {m_ImGuiSrvDescHeap.Get()};
        cmd->SetDescriptorHeaps(static_cast<UINT>(std::size(descHeaps)), descHeaps);
        ImGui_ImplDX12_RenderDrawData(ImGui::GetDrawData(), cmd->GetD3D12CmdList());
    }

    std::swap(barrier.Transition.StateBefore, barrier.Transition.StateAfter);
    cmd->ResourceBarrier(barrier);

    GetGraphicsDevice().FrameExecute();

    // Update and Render additional Platform Windows
    if (ImGui::GetIO().ConfigFlags & ImGuiConfigFlags_ViewportsEnable)
    {
        ImGui::UpdatePlatformWindows();
        ImGui::RenderPlatformWindowsDefault(nullptr, (void *)cmd);
    }

    m_Swapchain->Present(0, 0);
    SetNeedWaitSwapchain();
    GetGraphicsDevice().FrameSignal();
}

void RenderManager::DropEvent(const WCHAR *path)
{
    m_Renderer->LoadResource(path);
}

void RenderManager::WMResize(UINT width, UINT height)
{
    if (m_Swapchain == nullptr)
        return;
    m_ScreenSize = {width, height};
    GetGraphicsDevice().WaitForLastSubmittedFrame(GetAndResetSwapchainWaitableObject());
    for (int i = 0; i < kNumBackBuffer; i++)
        m_SwapchainBuffer[i].Reset();
    m_Swapchain->ResizeBuffers(0, width, height, DXGI_FORMAT_UNKNOWN,
                               DXGI_SWAP_CHAIN_FLAG_FRAME_LATENCY_WAITABLE_OBJECT);
    CreateSwapchainBuffer();
}

void RenderManager::CreateSwapchainBuffer()
{
    auto device = GetGraphicsDevice().GetD3D12Device();
    auto rtvDescriptorSize = device->GetDescriptorHandleIncrementSize(D3D12_DESCRIPTOR_HEAP_TYPE_RTV);
    auto rtvHandle = m_RtvDescHeap->GetCPUDescriptorHandleForHeapStart();
    for (int i = 0; i < kNumBackBuffer; i++)
    {
        m_Swapchain->GetBuffer(i, IID_PPV_ARGS(&m_SwapchainBuffer[i]));
        m_SwapchainBufferHandle[i] = rtvHandle;
        device->CreateRenderTargetView(m_SwapchainBuffer[i].Get(), nullptr, rtvHandle);
        rtvHandle.ptr += rtvDescriptorSize;
    }
}

