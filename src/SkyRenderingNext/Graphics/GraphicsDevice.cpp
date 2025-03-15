#include "GraphicsDevice.h"

#include <fstream>
#include <sstream>
#include <string_view>

#include "Texture2D.h"

#define UNWRAP(exp)                                                                                                    \
    do                                                                                                                 \
    {                                                                                                                  \
        auto ret = exp;                                                                                                \
        ASSERT(#exp &&SUCCEEDED(ret));                                                                                 \
    } while (0)

GraphicsDevice &GetGraphicsDevice()
{
    static GraphicsDevice device;
    return device;
}

GraphicsDevice::GraphicsDevice()
{
    // Init DX12
#ifdef _DEBUG
    ID3D12Debug *dx12Debug = nullptr;
    UNWRAP(D3D12GetDebugInterface(IID_PPV_ARGS(&dx12Debug)));
    dx12Debug->EnableDebugLayer();
#endif

    UNWRAP(D3D12CreateDevice(nullptr, D3D_FEATURE_LEVEL_12_1, IID_PPV_ARGS(&m_Device)));

#ifdef _DEBUG
    ID3D12InfoQueue *dx12InfoQueue = nullptr;
    m_Device->QueryInterface(IID_PPV_ARGS(&dx12InfoQueue));
    dx12InfoQueue->SetBreakOnSeverity(D3D12_MESSAGE_SEVERITY_CORRUPTION, true);
    dx12InfoQueue->SetBreakOnSeverity(D3D12_MESSAGE_SEVERITY_ERROR, true);

    // disable warning break for DLSS
    // dx12InfoQueue->SetBreakOnSeverity(D3D12_MESSAGE_SEVERITY_WARNING, true);
    dx12InfoQueue->Release();
    dx12Debug->Release();
#endif

    {
        D3D12_COMMAND_QUEUE_DESC desc = {
            .Type = D3D12_COMMAND_LIST_TYPE_DIRECT, .Flags = D3D12_COMMAND_QUEUE_FLAG_NONE, .NodeMask = 1};

        m_Device->CreateCommandQueue(&desc, IID_PPV_ARGS(&m_CmdQueue));
    }

    {
        m_DescHeap.Init(m_Device.Get());
    }

    {

        for (int i = 0; i < kNumFrameInFlights; ++i)
        {
            m_Device->CreateCommandAllocator(D3D12_COMMAND_LIST_TYPE_DIRECT, IID_PPV_ARGS(&m_FrameContext[i].cmdAlloc));
        }
        m_TransientBuffer.Init(m_Device.Get());

        m_CmdList = std::make_unique<GraphicsCommandList>(m_Device.Get(), CurFrameCtx().cmdAlloc.Get());
        m_CmdList->Close();
    }

    {
        m_Device->CreateFence(0, D3D12_FENCE_FLAG_NONE, IID_PPV_ARGS(&m_Fence));

        m_FenceEvent = CreateEvent(nullptr, FALSE, FALSE, nullptr);
        ASSERT(m_FenceEvent != nullptr);
    }
}

GraphicsDevice::~GraphicsDevice()
{
    if (m_FenceEvent)
    {
        CloseHandle(m_FenceEvent);
        m_FenceEvent = nullptr;
    }
}

void GraphicsDevice::FrameWaitAndBegin(HANDLE additionalWaitHandle)
{
    WaitFrameInterval(kNumFrameInFlights - 1, additionalWaitHandle);

    m_FrameIndex++;

    CurFrameCtx().cmdAlloc->Reset();
    m_CmdList->GetD3D12CmdList()->Reset(CurFrameCtx().cmdAlloc.Get(), nullptr);

    m_DescHeap.SetFrame(m_CmdList.get());

    m_TransientBuffer.FrameUpdate(m_FrameIndex - 1);
    for (const auto &handle : CurFrameCtx().toBeFreedHandles)
    {
        m_DescHeap.Free(handle);
    }
    CurFrameCtx().toBeFreedHandles.clear();

    for (const auto &resource : CurFrameCtx().toBeFreedResources)
    {
        resource->Release();
    }
    CurFrameCtx().toBeFreedResources.clear();

    if (!m_BuildinTextures)
    {
        m_BuildinTextures = std::make_unique<BuildinTextures>(m_CmdList.get());
    }
}

void GraphicsDevice::FrameExecute()
{
    auto cmd = CurCmdList();

    cmd->Close();

    ID3D12CommandList *cmdsLists[] = {cmd->GetD3D12CmdList()};
    m_CmdQueue->ExecuteCommandLists(static_cast<UINT>(std::size(cmdsLists)), cmdsLists);
}

void GraphicsDevice::FrameSignal()
{
    m_CmdQueue->Signal(m_Fence.Get(), m_FrameIndex);
}

void GraphicsDevice::WaitForLastSubmittedFrame(HANDLE additionalWaitHandle)
{
    WaitFrameInterval(0, additionalWaitHandle);
}

void GraphicsDevice::WaitFrameInterval(int interval, HANDLE additionalWaitHandle)
{
    DWORD numWaitableObjects = 0;
    HANDLE waitableObjects[] = {nullptr, nullptr};
    if (additionalWaitHandle)
    {
        waitableObjects[numWaitableObjects++] = additionalWaitHandle;
    }
    if (m_Fence->GetCompletedValue() + interval < m_FrameIndex)
    {
        m_Fence->SetEventOnCompletion(m_FrameIndex - interval, m_FenceEvent);
        waitableObjects[numWaitableObjects++] = m_FenceEvent;
    }
    WaitForMultipleObjects(numWaitableObjects, waitableObjects, TRUE, INFINITE);
}

DescriptorHandle GraphicsDevice::AllocTransientBufferDescriptor(UINT offset, UINT size)
{
    auto descriptorIndex = AllocBindlessDescriptor();
    D3D12_CONSTANT_BUFFER_VIEW_DESC desc{
        .BufferLocation = m_TransientBuffer.GetGpuAddr(offset),
        .SizeInBytes = size,
    };
    m_Device->CreateConstantBufferView(&desc, GetBindlessCpuHandle(descriptorIndex));

    return descriptorIndex;
}

DescriptorHandle::~DescriptorHandle()
{
    if (!handle.Empty())
    {
        GetGraphicsDevice().CurFrameCtx().toBeFreedHandles.emplace_back(handle);
        handle = {};
    }
}

void DelayRelease(IUnknown *ptr)
{
    if (ptr)
    {
#if 1
        GetGraphicsDevice().CurFrameCtx().toBeFreedResources.emplace_back(ptr);
#else
        ptr->Release();
#endif
    }
}
