#pragma once

#include "BindlessDescriptorHeap.h"
#include "TransientBuffer.h"
#include "CommandList.h"

static constexpr int kNumFrameInFlights = 3;

class GraphicsDevice;
class BuildinTextures;

class DescriptorHandle
{
  public:
    DescriptorHandle() = default;

    DescriptorHandle(BindlessDescriptorHeap::Handle handle) : handle(handle)
    {
    }

    ~DescriptorHandle();

    DescriptorHandle(const DescriptorHandle &) = delete;

    DescriptorHandle(DescriptorHandle &&rhs) noexcept : DescriptorHandle()
    {
        swap(rhs);
    }

    DescriptorHandle &operator=(DescriptorHandle rhs) noexcept
    {
        swap(rhs);
        return *this;
    }

    UINT index() const noexcept
    {
        return handle.index;
    }

    void swap(DescriptorHandle &rhs) noexcept
    {
        std::swap(handle, rhs.handle);
    }

    bool Empty() const noexcept
    {
        return handle.Empty();
    }

  private:
    friend GraphicsDevice;

    BindlessDescriptorHeap::Handle handle;
};

void DelayRelease(IUnknown *ptr);

template <class T> class DelayReleasePtr
{
  public:
    DelayReleasePtr() = default;

    ~DelayReleasePtr()
    {
        DelayRelease(ptr);
        ptr = {};
    }

    DelayReleasePtr(const DelayReleasePtr &) = delete;

    DelayReleasePtr(DelayReleasePtr &&rhs) noexcept : DelayReleasePtr()
    {
        swap(rhs);
    }

    DelayReleasePtr &operator=(DelayReleasePtr rhs) noexcept
    {
        swap(rhs);
        return *this;
    }

    T *Get() noexcept
    {
        return ptr;
    }

    T **operator&()
    {
        *this = {};
        return &ptr;
    }

    T *operator->() noexcept
    {
        return ptr;
    }

    operator bool() const noexcept
    {
        return ptr;
    }

    void swap(DelayReleasePtr &rhs) noexcept
    {
        std::swap(ptr, rhs.ptr);
    }

  private:
    T* ptr{};
};

template<class T>
using ResPtr = DelayReleasePtr<T>;

class GraphicsDevice
{
  public:
    GraphicsDevice();

    ~GraphicsDevice();

    GraphicsDevice(const GraphicsDevice &) = delete;
    GraphicsDevice &operator=(const GraphicsDevice &) = delete;

    ID3D12Device5 *GetD3D12Device() const
    {
        return m_Device.Get();
    }

    ID3D12CommandQueue *GetCmdQueue() const
    {
        return m_CmdQueue.Get();
    }

    GraphicsCommandList *CurCmdList() const
    {
        return m_CmdList.get();
    }

    ID3D12RootSignature *GetBindlessRootSignature() const
    {
        return m_DescHeap.GetRootSignature();
    }

    DescriptorHandle AllocBindlessDescriptor()
    {
        return m_DescHeap.Alloc();
    }

    D3D12_CPU_DESCRIPTOR_HANDLE GetBindlessCpuHandle(const DescriptorHandle &handle)
    {
        return m_DescHeap.GetCpuHandle(handle.handle);
    }

    template <class T> DescriptorHandle AllocTransientBuffer(const T &data)
    {
        auto range = m_TransientBuffer.AllocAndWrite(data);
        return AllocTransientBufferDescriptor(range.first, range.second);
    }

    template <class T> UINT64 AllocTransientBufferAddress(const T &data)
    {
        auto range = m_TransientBuffer.AllocAndWrite(data);
        return m_TransientBuffer.GetGpuAddr(range.first);
    }

    const BuildinTextures& GetBuildinTextures() const
    {
        return *m_BuildinTextures;
    }

    UINT64 GetFrameIndex() const
    {
        return m_FrameIndex;
    }

    UINT64 GetCompleteFrameIndex() const
    {
        return m_Fence->GetCompletedValue();
    }

    void FrameWaitAndBegin(HANDLE additionalWaitHandle);

    void FrameExecute();

    void FrameSignal();

    void WaitForLastSubmittedFrame(HANDLE additionalWaitHandle);

  protected:
    friend DescriptorHandle;
    friend void DelayRelease(IUnknown *ptr);

    void WaitFrameInterval(int interval, HANDLE additionalWaitHandle);

    struct FrameContext
    {
        ComPtr<ID3D12CommandAllocator> cmdAlloc;
        std::vector<BindlessDescriptorHeap::Handle> toBeFreedHandles;
        std::vector<IUnknown *> toBeFreedResources;
    };

    DescriptorHandle AllocTransientBufferDescriptor(UINT offset, UINT size);

    FrameContext &CurFrameCtx()
    {
        return m_FrameContext[m_FrameIndex % kNumFrameInFlights];
    }

    ComPtr<ID3D12Device5> m_Device;
    ComPtr<ID3D12CommandQueue> m_CmdQueue;
    std::unique_ptr<GraphicsCommandList> m_CmdList;
    FrameContext m_FrameContext[kNumFrameInFlights];
    TransientBuffer m_TransientBuffer;

    BindlessDescriptorHeap m_DescHeap;
    std::unique_ptr<BuildinTextures> m_BuildinTextures;

    ComPtr<ID3D12Fence> m_Fence;
    HANDLE m_FenceEvent = nullptr;
    UINT64 m_FrameIndex = 0;
};

GraphicsDevice &GetGraphicsDevice();
