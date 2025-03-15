#pragma once

#include <dxgi1_6.h>

#include "Core/Common.h"
#include "Core/Math.h"
#include "Renderer.h"
#include "Camera.h"

static constexpr DXGI_FORMAT kBackBufferFormat = DXGI_FORMAT_R8G8B8A8_UNORM;

class RenderManager &GetRenderManager();

class RenderManager
{
    RenderManager();
    friend RenderManager &GetRenderManager();

  public:
    static constexpr int kNumBackBuffer = 3;

    RenderManager(const RenderManager &) = delete;
    RenderManager &operator=(const RenderManager &) = delete;

    ~RenderManager();

    Math::uint2 GetScreenSize() const
    {
        return m_ScreenSize;
    }

    Math::uint4 GetWindowAbsoluteRect() const;

    void MainLoop(int argc, const char *argv[]);

    void LoopBody();

    void DropEvent(const WCHAR* path);

    void WMResize(UINT width, UINT height);

    HWND GetHwnd() const
    {
        return m_Hwnd;
    }

    Camera &MainCamera()
    {
        return m_MainCamera;
    }

  protected:
    void CreateSwapchainBuffer();

    HANDLE GetAndResetSwapchainWaitableObject()
    {
        HANDLE res = nullptr;
        if (m_SwapchainNeedWait)
        {
            res = m_SwapchainWaitableObject;
            m_SwapchainNeedWait = false;
        }
        return res;
    }

    void SetNeedWaitSwapchain()
    {
        m_SwapchainNeedWait = true;
    }

    HWND m_Hwnd;
    ComPtr<ID3D12DescriptorHeap> m_ImGuiSrvDescHeap;
    ComPtr<ID3D12DescriptorHeap> m_RtvDescHeap;
    HANDLE m_SwapchainWaitableObject = nullptr;
    bool m_SwapchainNeedWait = false;
    ComPtr<IDXGISwapChain3> m_Swapchain;
    ComPtr<ID3D12Resource> m_SwapchainBuffer[kNumBackBuffer];
    D3D12_CPU_DESCRIPTOR_HANDLE m_SwapchainBufferHandle[kNumBackBuffer];

    Math::uint2 m_ScreenSize{1280, 800};
    std::unique_ptr<Renderer> m_Renderer;

    Camera m_MainCamera;
};
