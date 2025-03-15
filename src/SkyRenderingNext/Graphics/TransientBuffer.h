#pragma once

#include <deque>

#include "Core/Common.h"

class TransientBuffer
{
  public:
    void Init(ID3D12Device *device);

    void FrameUpdate(UINT64 lastSubmittedFrameIndex);

    std::pair<UINT, UINT> Alloc(UINT size);

    void *GetCpuPtr(UINT offset) const
    {
        return reinterpret_cast<UINT8 *>(m_BasePtr) + offset;
    }

    UINT64 GetGpuAddr(UINT offset) const
    {
        return m_GpuBaseAddr + offset;
    }

    template<class T> 
    std::pair<UINT, UINT> AllocAndWrite(const T &data)
    {
        auto rawSize = static_cast<UINT>(sizeof(data));
        auto res = Alloc(rawSize);
        memcpy(GetCpuPtr(res.first), &data, rawSize);
        return res;
    }

  private:
    static UINT GetOffset(UINT rawOffset)
    {
        return rawOffset & (kBufferSize - 1);
    }

    static constexpr UINT kBufferSize = 128u * 1024 * 1024;
    static_assert((kBufferSize & (kBufferSize - 1)) == 0);

    ComPtr<ID3D12Resource> m_Buffer{};
    void *m_BasePtr{};
    UINT64 m_GpuBaseAddr{};
    UINT m_OffsetHead{};
    UINT m_OffsetTail{};

    struct FrameOffset
    {
        UINT64 frameIndex;
        UINT32 offset;
    };
    std::deque<FrameOffset> m_ToBeFreed;
};
